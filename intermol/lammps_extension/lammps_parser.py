"""LAMMPS extension for InterMol.

.. module:: lammps_parser
.. moduleauthor:: Christoph Klein <christoph.t.klein@me.com>
"""
import os
import pdb
import logging
import numpy as np
from collections import OrderedDict

import intermol.unit as units
from intermol.atom import Atom
from intermol.molecule import Molecule
from intermol.system import System
from intermol.types import *
from intermol.forces import *

logger = logging.getLogger('InterMolLog')

class LammpsParser(object):
    """A class containing methods to read and write LAMMPS files."""

    def __init__(self):
        """
        """

        self.box_vector = np.zeros(shape=(3, 3), dtype=float)

    def read_system(self, input_file):
        """Reads a LAMMPS input file and a data file specified within.

        Args:
            input_file (str): Name of LAMMPS input file to read in.
        """
        self.read_input(input_file)
        if self.data_file:
            self.read_data(self.data_file)
        else:
            raise Exception("No data file found in input script")

    def read_input(self, input_file):
        """Reads a LAMMPS input file.

        Args:
            input_file (str): Name of LAMMPS input file to read in.
        """
        self.basepath = os.path.dirname(os.path.realpath(input_file))
        parsable_keywords = {'units': self.parse_units,
                'atom_style': self.parse_atom_style,
                'dimension': self.parse_dimension,
                'boundary': self.parse_boundary,
                'pair_style': self.parse_pair_style,
                'kspace_style': self.parse_kspace_style,
                'pair_modify': self.parse_pair_modify,
                'bond_style': self.parse_bond_style,
                'angle_style': self.parse_angle_style,
                'dihedral_style': self.parse_dihedral_style,
                'improper_style': self.parse_improper_style,
                'special_bonds': self.parse_special_bonds,
                'read_data': self.parse_read_data,
                'fix': self.parse_fix}

        with open(input_file, 'r') as input_lines:
            for line in input_lines:
                if line.strip():
                    keyword = line.split()[0]
                    if keyword in parsable_keywords:
                        parsable_keywords[keyword](line.split())

    def read_data(self, data_file):
        """Reads a LAMMPS data file.

        Args:
            data_file (str): name of LAMMPS data file to read in.
        """
        # Read box, masses and forcefield info from data file.
        parsable_keywords = {'Masses': self.parse_masses,
                'Pair Coeffs': self.parse_pair_coeffs,
                'Bond Coeffs': self.parse_bond_coeffs,
                'Angle Coeffs': self.parse_angle_coeffs,
                'Dihedral Coeffs': self.parse_dihedral_coeffs,
                'Improper Coeffs': self.parse_improper_coeffs}

        with open(data_file, 'r') as data_lines:
            for line in data_lines:
                if line.strip():
                    line = line.partition('#')[0] # Remove trailing comment

                    fields = line.split()
                    if len(fields) == 2 and fields[1] == 'atoms':
                        self.natoms = int(fields[0])

                    # catch all box dimensions
                    if (('xlo' in line) and
                         ('xhi' in line)):
                        self.parse_box(line.split(), 0)
                    elif (('ylo' in line) and
                         ('yhi' in line)):
                        self.parse_box(line.split(), 1)
                    elif (('zlo' in line) and
                         ('zhi' in line)):
                        self.parse_box(line.split(), 2)
                    # other headers
                    else:
                        keyword = line.strip()
                        if keyword in parsable_keywords:
                            parsable_keywords[keyword](data_lines)

        # Read atoms and connectivity information from data file.
        parsable_keywords = {'Atoms': self.parse_atoms,
                'Bonds': self.parse_bonds,
                'Angles': self.parse_angles,
                'Dihedrals': self.parse_dihedrals,
                'Impropers': self.parse_impropers}

        with open(data_file, 'r') as data_lines:
            for line in data_lines:
                if line.strip():
                    keyword = line.partition('#')[0].strip()
                    if keyword in parsable_keywords:
                        parsable_keywords[keyword](data_lines)

        # SETTLE constraints (for now limited to rigid 3-point light water)
        for mol_type in System._sys._molecules.itervalues():
            if (len(mol_type.angleForceSet) == 1 and 
                len(mol_type.bondForceSet) == 2):
                mol_mass = units.sum([a._mass[0] for a in mol_type.moleculeSet[0].getAtoms()])
                if np.round(mol_mass.in_units_of(units.amu)._value) != 18.0:
                    continue

                is_rigid_water = True
                for angle in mol_type.angleForceSet.itervalues():
                    if isinstance(angle, Angle) and angle.c:
                        thetaHOH = angle.theta
                    else:
                        is_rigid_water = False
                for bond in mol_type.bondForceSet.itervalues():
                    if isinstance(bond, Bond) and bond.c:
                        dOH = bond.length
                    else:
                        is_rigid_water = False
                
                if is_rigid_water:
                    dHH = units.sqrt(2 * dOH**2 * (1 - units.cos(thetaHOH)))
                    for i, atom in enumerate(mol_type.moleculeSet[0]._atoms, start=1):
                        if atom._mass[0].in_units_of(units.amu)._value > 15:
                            iOW = i
                            break
                    mol_type.settles = Settles(iOW, dOH, dHH)
                    mol_type.nrexcl = 1
                    mol_type.exclusions.add(Exclusions([2, 3]))
                    mol_type.exclusions.add(Exclusions([1, 3]))
                    mol_type.exclusions.add(Exclusions([1, 2]))
        
        # Indentify 1-3 and 1-4 neighbors
        onethree = [[] for i in range(self.natoms + 1)]
        onefour  = [[] for i in range(self.natoms + 1)]
        for i in range(1, self.natoms):
            # 1-3 neighbors
            for j in self.onetwo[i]:
                for k in self.onetwo[j]:
                    if not ((k == i) or (k in self.onetwo[i])):
                        onethree[i].append(k)
            
            # 1-4 neighbors
            for j in onethree[i]:
                for k in self.onetwo[j]:
                    if not ((k == i) or (k in self.onetwo[i]) or (k in onethree[i])):
                        onefour[i].append(k)
                                    

    def parse_units(self, line):
        """ """
        assert(len(line) == 2), "Invalid units specified in input file."
        self.unit_set = line[1]

        self.RAD = units.radians
        self.DEGREE = units.degrees
        if self.unit_set == 'real':
            self.DIST = units.angstroms
            self.VEL = units.angstroms / units.femtosecond
            self.ENERGY = units.kilocalorie / units.mole
            self.MASS = units.grams / units.mole
            self.CHARGE = units.elementary_charge
            self.MOLE = units.mole
        else:
            raise Exception("Unsupported unit set specified in input file: "
                    "{0}".format(self.unit_set))

    def parse_atom_style(self, line):
        """
        """
        if len(line) > 2:
            logger.warn("Unsupported atom_style in input file.")
        self.atom_style = line[1]

        # Currently only support atom_style 'full'
        if self.atom_style != 'full':
            raise Exception("Unsupported atom_style in input file: {0}".format(self.atom_style))        

    def parse_dimension(self, line):
        """ """
        self.dimension = int(line[1])
        if self.dimension not in [2, 3]:
            raise ValueError("Invalid dimension specified in input file"
                    " (must be 2 or 3).")

    def parse_boundary(self, line):
        """ """
        self.boundaries = [line[1], line[2], line[3]]
        if len(self.boundaries) != self.dimension:
            raise ValueError("Boundaries do not match specified dimension "
                    "in input file")

    def parse_pair_style(self, line):
        """ """
        self.pair_style = []
        if line[1] == 'hybrid':
            logger.warn("Hybrid pair styles not yet implemented.")
        elif line[1] == 'lj/cut/coul/long':
            self.pair_style.append(line[1])
            System._sys.nonbonded_function = 1

    def parse_kspace_style(self, line):
        """
        Note:
            Currently ignored.
        """
        if line[1] == 'pppm':
            pass

    def parse_pair_modify(self, line):
        """
        """
        if line[1] == 'mix':
            if line[2] == 'geometric':
                System._sys.combination_rule = 3
            elif line[2] == 'arithmetic':
                System._sys.combination_rule = 2
            else:
                logger.warn("Unsupported pair_modify mix argument in input file!")
        else:
            logger.warn("Unsupported pair_modify style in input file!")

    def parse_bond_style(self, line):
        """ """
        self.bond_style = []
        self.hybrid_bond_style = False
        if len(line) == 2:
            self.bond_style.append(line[1])
        elif len(line) > 2 and line[1] == 'hybrid':
            self.hybrid_bond_style = True
            for style in line[2:]:
                self.bond_style.append(style)
        else:
            raise ValueError("Invalid bond_style in input file!")

    def parse_angle_style(self, line):
        """ """
        self.angle_style = []
        self.hybrid_angle_style = False
        if len(line) == 2:
            self.angle_style.append(line[1])
        elif len(line) > 2 and line[1] == 'hybrid':
            self.hybrid_angle_style = True
            for style in line[2:]:
                self.angle_style.append(style)
        else:
            raise ValueError("Invalid angle_style in input file!")

    def parse_dihedral_style(self, line):
        """ """
        self.dihedral_style = []
        self.hybrid_dihedral_style = False
        if len(line) == 2:
            self.dihedral_style.append(line[1])
            # TODO: correctly determine gen-pairs state
            if self.dihedral_style == 'opls':
                System._sys.genpairs = 'yes'
        elif len(line) > 2 and line[1] == 'hybrid':
            self.hybrid_dihedral_style = True
            for style in line[2:]:
                self.dihedral_style.append(style)
        else:
            raise ValueError("Invalid dihedral_style in input file!")

    def parse_improper_style(self, line):
        """ """
        self.improper_style = []
        self.hybrid_improper_style = False
        if len(line) == 2:
            self.improper_style.append(line[1])
        elif len(line) > 2 and line[1] == 'hybrid':
            self.hybrid_improper_style = True
            for style in line[2:]:
                self.improper_style.append(style)
        else:
            raise ValueError("Invalid improper_style in input file!")

    def parse_special_bonds(self, line):
        """ """
        if 'lj/coul' in line:
            System._sys.lj_correction = float(line[line.index('lj/coul') + 3])
            System._sys.coulomb_correction = float(line[line.index('lj/coul') + 3])
        elif 'lj' in line and 'coul' in line:
            System._sys.lj_correction = float(line[line.index('lj') + 3])
            System._sys.coulomb_correction = float(line[line.index('coul') + 3])
        elif 'lj' in line:
            System._sys.lj_correction = float(line[line.index('lj') + 3])
        elif 'coul' in line:
            System._sys.coulomb_correction = float(line[line.index('coul') + 3])
        else:
            logger.warn("Unsupported special_bonds in input file.")

    def parse_read_data(self, line):
        """ """
        if len(line) == 2:
            self.data_file = os.path.join(self.basepath, line[1])
        else:
            logger.warn("Unsupported read_data arguments in input file.")

    def parse_fix(self, line):
        """ """
        if len(line) > 3 and line[3] == 'shake':
            self.parse_fix_shake(line)

    def parse_fix_shake(self, line):
        """ """
        if line[2] != 'all':
            logger.warn("Unsupported group-ID in fix shake command")
            return            
        if 'mol' in line:
            logger.warn("Unsupported keyword 'mol' in fix shake command")
            return
        
        self.shake_bond_types_i = set()
        self.shake_angle_types_i = set()
        self.shake_masses = set()

        line = line[7:]
        for field in line:
            if field == 't':
                logger.warn("SHAKE 't' (atom type) constraints not yet supported: fix shake ignored.")
                return
            elif field == 'b':
                container = self.shake_bond_types_i
                str2num = int
            elif field == 'a':
                container = self.shake_angle_types_i
                str2num = int
            elif field == 'm':
                container = self.shake_masses
                str2num = float
            else:
                container.add(str2num(field))
            
    def parse_box(self, line, dim):
        """Read box information from data file.

        Args:
            line (str): Current line in input file.
            dim (int): Dimension specified in line.
        """
        fields = [float(field) for field in line[:2]]
        box_length = fields[1] - fields[0]
        if box_length > 0:
            self.box_vector[dim, dim] = box_length
        else:
            raise ValueError("Negative box length specified in data file.")
        System._sys.box_vector = self.box_vector * self.DIST

    def parse_masses(self, data_lines):
        """Read masses from data file."""
        next(data_lines)  # toss out blank line
        self.mass_dict = dict()
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = line.partition('#')[0].split()
            self.mass_dict[int(fields[0])] = float(fields[1]) * self.MASS

    def parse_pair_coeffs(self, data_lines):
        """Read pair coefficients from data file."""
        next(data_lines)  # toss out blank line
        self.nb_types = dict()
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = [float(field) for field in line.partition('#')[0].split()]
            if len(self.pair_style) == 1:
                # TODO: lookup of type of pairstyle to determine format
                if System._sys.nonbonded_function == 1:
                    self.nb_types[int(fields[0])] = [fields[1] * self.ENERGY,
                                                     fields[2] * self.DIST]
                else:
                    logger.warn("Unsupported pair coeff formatting in data file!")
            else:
                logger.warn("Unsupported pair coeff formatting in data file!")

    def parse_bond_coeffs(self, data_lines):
        """Read bond coefficients from data file."""
        next(data_lines)  # toss out blank line
        self.bond_types = dict()
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = line.partition('#')[0].split()
            if self.hybrid_bond_style:
                style = fields[1]
                if style not in self.bond_style:
                    raise Exception("Bond type found in Bond Coeffs that "
                            "was not specified in bond_style: {0}".format(style))
                coeffs = fields[2:]
            else:
                style = self.bond_style[0]
                coeffs = fields[1:]

            if style == 'harmonic':
                self.bond_types[int(fields[0])] = [
                    style,
                    2 * float(coeffs[0]) * self.ENERGY / (self.DIST*self.DIST),
                    float(coeffs[1]) * self.DIST]
            elif style == 'morse':
                self.bond_types[int(fields[0])] = [
                    style,
                    float(coeffs[0]) * self.ENERGY,
                    float(coeffs[1]) * self.DIST**(-1),
                    float(coeffs[2]) * self.DIST]
            else:
                logger.warn("Unsupported bond style: {0}".format(style))

    def parse_angle_coeffs(self, data_lines):
        """Read angle coefficients from data file."""
        next(data_lines)  # toss out blank line
        self.angle_types = dict()
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = line.partition('#')[0].split()
            if self.hybrid_angle_style:
                style = fields[1]
                if style not in self.angle_style:
                    raise Exception("Angle type found in Angle Coeffs that "
                            "was not specified in angle_style: {0}".format(style))
                coeffs = fields[2:]
            else:
                style = self.angle_style[0]
                coeffs = fields[1:]

            if style == 'harmonic':
                self.angle_types[int(fields[0])] = [
                    style,
                    2 * float(coeffs[0]) * self.ENERGY / self.RAD**2,
                    float(coeffs[1]) * self.DEGREE]
            else:
                logger.warn("Unsupported angle style: {0}".format(style))

    def parse_dihedral_coeffs(self, data_lines):
        """Read dihedral coefficients from data file."""
        next(data_lines)  # toss out blank line
        self.dihedral_types = dict()
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = line.partition('#')[0].split()
            if self.hybrid_dihedral_style:
                style = fields[1]
                if style not in self.dihedral_style:
                    raise Exception("Dihedral type found in Dihedral Coeffs that "
                            "was not specified in dihedral_style: {0}".format(style))
                coeffs = fields[2:]
            else:
                style = self.dihedral_style[0]
                coeffs = fields[1:]

            if style == 'opls':
                self.dihedral_types[int(fields[0])] = [
                    style,
                    float(coeffs[0]) * self.ENERGY,
                    float(coeffs[1]) * self.ENERGY,
                    float(coeffs[2]) * self.ENERGY,
                    float(coeffs[3]) * self.ENERGY]
            elif style == 'multi/harmonic':
                self.dihedral_types[int(fields[0])] = [
                    style,
                    float(coeffs[0]) * self.ENERGY,
                    float(coeffs[1]) * self.ENERGY,
                    float(coeffs[2]) * self.ENERGY,
                    float(coeffs[3]) * self.ENERGY,
                    float(coeffs[4]) * self.ENERGY]
            else:
                logger.warn("Unsupported dihedral style: {0}".format(style))

    def parse_improper_coeffs(self, data_lines):
        """Read improper coefficients from data file."""
        next(data_lines)  # toss out blank line
        self.improper_types = dict()
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = line.partition('#')[0].split()
            if self.hybrid_improper_style:
                style = fields[1]
                if style not in self.improper_style:
                    raise Exception("Improper type found in Improper Coeffs that "
                            "was not specified in improper_style: {0}".format(style))
                coeffs = fields[2:]
            else:
                style = self.improper_style[0]
                coeffs = fields[1:]

            if style == 'harmonic':
                self.improper_types[int(fields[0])] = [
                    style,
                    float(coeffs[0]) * self.ENERGY / self.RAD**2,
                    float(coeffs[1]) * self.DEGREE]
            else:
                logger.warn("Unsupported improper style: {0}".format(style))

    def parse_atoms(self, data_lines):
        """Read atoms from data file."""
        molecules = dict()
        next(data_lines)  # toss out blank line
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = line.partition('#')[0].split()

            if len(fields) in [7, 10]:
                if len(fields) == 10:
                    # TODO: store image flags?
                    pass
                new_atom_type = None
                bondtype = atomtype = 'lmp_{:03d}'.format(int(fields[2]))
                if System._sys.combination_rule == 1:
                    logger.warn("Combination rule '1' not yet implemented")
                elif System._sys.combination_rule in [2, 3]:
                    new_atom_type = AtomCR23Type(atomtype, bondtype,
                        -1,                               # atomic_number
                        self.mass_dict[int(fields[2])],
                        float(fields[3]) * self.CHARGE,
                        'A',                              # ptype
                        self.nb_types[int(fields[2])][1], # sigma
                        self.nb_types[int(fields[2])][0]) # epsilon

                System._sys._atomtypes.add(new_atom_type)

                atom = Atom(int(fields[0]), 
                            '',             # atom name (set below)
                            1)              # residue index
                atom.setAtomType(0, atomtype)
                atom.cgnr = 1
                # TODO: Set cg nr 'automatically' in gromacs_topology_parser.
                # See, e.g., topogromacs.tcl in VMD topotools.
                atom.setCharge(0, float(fields[3]) * self.CHARGE)
                atom.setMass(0, self.mass_dict[int(fields[2])])
                atom.setPosition(float(fields[4]) * self.DIST,
                        float(fields[5]) * self.DIST,
                        float(fields[6]) * self.DIST)

                for ab_state, atom_type in enumerate(atom._atomtype):
                    # Searching for a matching atom_type
                    temp = AbstractAtomType(atom._atomtype[ab_state])
                    atom_type = System._sys._atomtypes.get(temp)
                    if atom_type:
                        atom.setSigma(ab_state, atom_type.sigma)
                        atom.setEpsilon(ab_state, atom_type.epsilon)
                        atom.bondtype = atom_type.bondtype
                    else:
                        logger.warn("Corresponding AtomType was not found. "
                                "Insert missing values yourself.")

                if int(fields[1]) not in molecules:
                    molecules[int(fields[1])] = Molecule()
                molecules[int(fields[1])].addAtom(atom)

        # Add molecules to system
        atomtype_list_old = []
        moleculetype_i = 0
        mol_name = None
        self.nr = dict() # atom index => index within moleculetype (i.e. nr)
        self.mol_type = dict() # atom index => MoleculeType 
        for molecule in [molecules[i] for i in sorted(molecules.keys())]:
            molecule.getAtoms().list.sort()

            atomtype_list = [atom.getAtomType(0) for atom in molecule.getAtoms()]
            if atomtype_list != atomtype_list_old: # new moleculetype
                moleculetype_i += 1
                mol_name = 'moleculetype{:02d}'.format(moleculetype_i)
                atomtype_list_old = atomtype_list

            molecule.name = mol_name
            System._sys.add_molecule(molecule)

            # TODO: Move this elsewhere and determine nrexcl from special_bonds
            System._sys._molecules[molecule.name].nrexcl = 3

            for i, atom in enumerate(molecule.getAtoms(), start=1):
                self.nr[atom.index] = i
                self.mol_type[atom.index] = System._sys._molecules[molecule.name]
                atom.residue_name = 'R{:02d}'.format(moleculetype_i)
                atom.name = 'X{:d}'.format(i)

    def parse_bonds(self, data_lines):
        """Read bonds from data file."""
        self.onetwo = [[] for i in range(self.natoms + 1)] # 1-2 neighbors
        next(data_lines)  # toss out blank line
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = [int(field) for field in line.partition('#')[0].split()]

            new_bond_force = None
            coeff_num = fields[1]
            ai = self.nr[fields[2]]
            aj = self.nr[fields[3]]
            if coeff_num in self.shake_bond_types_i:
                constrained = True
            else:
                constrained = False

            # Bond
            if self.bond_types[coeff_num][0] == 'harmonic':
                r = self.bond_types[coeff_num][2]
                k = self.bond_types[coeff_num][1]
                new_bond_force = Bond(ai, aj, r, k, c=constrained)
            # Morse
            elif self.bond_types[coeff_num][0] == 'morse':
                r = self.bond_types[coeff_num][3]
                D = self.bond_types[coeff_num][1]
                beta = self.bond_types[coeff_num][2]
                new_bond_force = MorseBond(ai, aj, r, D, beta)
            
            self.current_mol_type = self.mol_type[fields[2]]
            if (self.mol_type[fields[3]] is not self.current_mol_type):
                raise Exception("Attempted to define bond between atoms {0:d}, {1:d} across different molecules.".format(fields[2], fields[3]))
            self.current_mol_type.bondForceSet.add(new_bond_force)

            # Keep track of 1st neighbors
            self.onetwo[fields[2]].append(fields[3])
            self.onetwo[fields[3]].append(fields[2])

    def parse_angles(self, data_lines):
        """Read angles from data file."""
        next(data_lines)  # toss out blank line
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = [int(field) for field in line.partition('#')[0].split()]

            new_angle_force = None
            coeff_num = fields[1]
            ai = self.nr[fields[2]]
            aj = self.nr[fields[3]]
            ak = self.nr[fields[4]]
            if coeff_num in self.shake_angle_types_i:
                constrained = True
            else:
                constrained = False            

            # Angle
            if self.angle_types[coeff_num][0] == 'harmonic':
                theta = self.angle_types[coeff_num][2]
                k = self.angle_types[coeff_num][1]
                new_angle_force = Angle(ai, aj, ak, theta, k, c=constrained)

            self.current_mol_type = self.mol_type[fields[2]]
            if ((self.mol_type[fields[3]] is not self.current_mol_type) or 
                (self.mol_type[fields[4]] is not self.current_mol_type)):
                raise Exception("Attempted to define angle between atoms {0:d}, {1:d}, {2:d} across different molecules.".format(fields[2], fields[3], fields[4]))
            self.current_mol_type.angleForceSet.add(new_angle_force)

    def parse_dihedrals(self, data_lines):
        """Read dihedrals from data file."""
        next(data_lines)  # toss out blank line
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = [int(field) for field in line.partition('#')[0].split()]

            new_dihed_force = None
            coeff_num = fields[1]
            ai = self.nr[fields[2]]
            aj = self.nr[fields[3]]
            ak = self.nr[fields[4]]
            al = self.nr[fields[5]]

            if  self.dihedral_types[coeff_num][0] == 'opls':
                fc0, fc1, fc2, fc3, fc4, fc5, fc6 = ConvertDihedralFromFourierToDihedralTrig(
                    self.dihedral_types[coeff_num][1],
                    self.dihedral_types[coeff_num][2],
                    self.dihedral_types[coeff_num][3],
                    self.dihedral_types[coeff_num][4])
            elif self.dihedral_types[coeff_num][0] == 'multi/harmonic':
                fc0, fc1, fc2, fc3, fc4, fc5, fc6 = ConvertDihedralFromRBToDihedralTrig(
                    self.dihedral_types[coeff_num][1],
                    -self.dihedral_types[coeff_num][2],
                    self.dihedral_types[coeff_num][3],
                    -self.dihedral_types[coeff_num][4],
                    self.dihedral_types[coeff_num][5],
                    0 * self.ENERGY,
                    0 * self.ENERGY)

            new_dihed_force = DihedralTrigDihedral(
                ai, aj, ak, al,
                0 * self.DEGREE, fc0, fc1, fc2, fc3, fc4, fc5, fc6)

            self.current_mol_type = self.mol_type[fields[2]]
            if ((self.mol_type[fields[3]] is not self.current_mol_type) or 
                (self.mol_type[fields[4]] is not self.current_mol_type) or 
                (self.mol_type[fields[5]] is not self.current_mol_type)):
                raise Exception("Attempted to define dihedral between atoms {0:d}, {1:d}, {2:d}, {3:d} across different molecules.".format(fields[2], fields[3], fields[4], fields[5]))
            self.current_mol_type.dihedralForceSet.add(new_dihed_force)

    def parse_impropers(self, data_lines):
        """Read impropers from data file."""
        next(data_lines)  # toss out blank line
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = [int(field) for field in line.partition('#')[0].split()]

            new_dihed_force = None
            coeff_num = fields[1]
            ai = self.nr[fields[2]]
            aj = self.nr[fields[3]]
            ak = self.nr[fields[4]]
            al = self.nr[fields[5]]

            if  self.improper_types[coeff_num][0] == 'harmonic':
                k = self.improper_types[coeff_num][1]
                xi = self.improper_types[coeff_num][2]
                new_dihed_force = ImproperHarmonicDihedral(ai, aj, ak, al, xi, k)

            self.current_mol_type = self.mol_type[fields[2]]
            if ((self.mol_type[fields[3]] is not self.current_mol_type) or 
                (self.mol_type[fields[4]] is not self.current_mol_type) or 
                (self.mol_type[fields[5]] is not self.current_mol_type)):
                raise Exception("Attempted to define improper dihedral between atoms {0:d}, {1:d}, {2:d}, {3:d} across different molecules.".format(fields[2], fields[3], fields[4], fields[5]))
            self.current_mol_type.dihedralForceSet.add(new_dihed_force)

    def write(self, data_file, unit_set='real', verbose=False):
        """Writes a LAMMPS data and corresponding input file.

        Args:
            data_file (str): Name of LAMMPS data file to write to.
            unit_set (str): LAMMPS unit set for output file.
        """
        self.RAD = units.radians
        self.DEGREE = units.degrees
        if unit_set == 'real':
            self.DIST = units.angstroms
            self.VEL = units.angstroms / units.femtosecond
            self.ENERGY = units.kilocalorie / units.mole
            self.MASS = units.grams / units.mole
            self.CHARGE = units.elementary_charge
            self.MOLE = units.mole
        else:
            raise Exception("Unsupported unit set specified: {0}".format(unit_set))

        # Containers for lines which are ultimately written to output files.
        mass_list = list()
        mass_list.append('\n')
        mass_list.append('Masses\n')
        mass_list.append('\n')

        pair_coeff_list = list()
        pair_coeff_list.append('\n')
        pair_coeff_list.append('Pair Coeffs\n')
        pair_coeff_list.append('\n')

        bond_coeffs = list()
        bond_coeffs.append('\n')
        bond_coeffs.append('Bond Coeffs\n')
        bond_coeffs.append('\n')

        angle_coeffs = list()
        angle_coeffs.append('\n')
        angle_coeffs.append('Angle Coeffs\n')
        angle_coeffs.append('\n')

        dihedral_coeffs = list()
        dihedral_coeffs.append('\n')
        dihedral_coeffs.append('Dihedral Coeffs\n')
        dihedral_coeffs.append('\n')

        improper_coeffs = list()
        improper_coeffs.append('\n')
        improper_coeffs.append('Improper Coeffs\n')
        improper_coeffs.append('\n')

        atom_list = list()
        atom_list.append('\n')
        atom_list.append('Atoms\n')
        atom_list.append('\n')

        vel_list = list()
        vel_list.append('\n')
        vel_list.append('Velocities\n')
        vel_list.append('\n')

        bond_list = list()
        bond_list.append('\n')
        bond_list.append('Bonds\n')
        bond_list.append('\n')

        angle_list = list()
        angle_list.append('\n')
        angle_list.append('Angles\n')
        angle_list.append('\n')

        dihedral_list = list()
        dihedral_list.append('\n')
        dihedral_list.append('Dihedrals\n')
        dihedral_list.append('\n')

        improper_list = list()
        improper_list.append('\n')
        improper_list.append('Impropers\n')
        improper_list.append('\n')

        # dicts for type information
        atom_type_dict = dict()  # str_type:int_type
        a_type_i = 1  # counter for atom types

        bond_style = set()
        bond_type_dict = dict()  # typeObject:int_type
        b_type_i = 1  # counter for bond types

        angle_style = set()
        angle_type_dict = dict()  # typeObject:int_type
        ang_type_i = 1

        dihedral_style = set()
        dihedral_type_dict = dict()  # typeObject:int_type
        dih_type_i = 1

        improper_style = set()
        improper_type_dict = dict()  # typeObject:int_type
        imp_type_i = 1

        # read all atom specific and FF information
        offset = 0
        bond_i = 1
        angle_i = 1
        dihedral_i = 1
        improper_i = 1
        shake_bond_types = set()
        shake_angle_types = set()
        x_min = y_min = z_min = np.inf
        for mol_type in System._sys._molecules.itervalues():
            logger.debug("    Writing moleculetype {0}...".format(mol_type.name))

            # Settles (for rigid water in GROMACS)
            # We'll convert these to SHAKE constraints in the input script.
            if mol_type.settles:
                for bond in mol_type.bondForceSet.itervalues():
                    shake_bond_types.add(BondType(
                            mol_type.moleculeSet[0]._atoms[bond.atom1 - 1].bondtype,
                            mol_type.moleculeSet[0]._atoms[bond.atom2 - 1].bondtype,
                            bond.length,
                            bond.k))
                for angle in mol_type.angleForceSet.itervalues():
                    shake_angle_types.add(AngleType(
                            mol_type.moleculeSet[0]._atoms[angle.atom1 - 1].bondtype,
                            mol_type.moleculeSet[0]._atoms[angle.atom2 - 1].bondtype,
                            mol_type.moleculeSet[0]._atoms[angle.atom3 - 1].bondtype,
                            angle.theta,
                            angle.k))

#            molecule = mol_type.moleculeSet[0]
#            atoms = molecule._atoms

#            for i, offset in enumerate(offsets):
            for molecule in mol_type.moleculeSet:

                # bonds
                logger.debug("        Writing bonds...")
                for bond in mol_type.bondForceSet.itervalues():
                    atomtype1 = molecule._atoms[bond.atom1 - 1].bondtype
                    atomtype2 = molecule._atoms[bond.atom2 - 1].bondtype

                    if isinstance(bond, Bond):
                        style = 'harmonic'
                        temp = BondType(atomtype1, atomtype2,
                                bond.length, bond.k)
                        # NOTE: k includes the factor of 0.5 for harmonic in LAMMPS
                        if temp not in bond_type_dict:
                            bond_type_dict[temp] = b_type_i
                            bond_coeffs.append('{0:d} {1} {2:18.8f} {3:18.8f} # {4:2s}-{5:2s}\n'.format(
                                    b_type_i,
                                    style,
                                    0.5 * bond.k.in_units_of(self.ENERGY / (self.DIST*self.DIST))._value,
                                    bond.length.in_units_of(self.DIST)._value,
                                    atomtype1, atomtype2))
                            b_type_i += 1
                    elif isinstance(bond, MorseBond):
                        style = 'morse'
                        temp = MorseBondType(atomtype1, atomtype2,
                                bond.length, bond.D, bond.beta)
                        if temp not in bond_type_dict:
                            bond_type_dict[temp] = b_type_i
                            bond_coeffs.append('{0:d} {1} {2:18.8f} {3:18.8f} {4:18.8f}\n'.format(
                                    b_type_i,
                                    style,
                                    bond.D.in_units_of(self.ENERGY)._value,
                                    bond.beta.in_units_of(self.DIST**(-1))._value,
                                    bond.length.in_units_of(self.DIST)._value))
                            b_type_i += 1
                    else:
                        logger.warn("Found unimplemented bond type for LAMMPS!")
                        continue

                    bond_list.append('{0:-6d} {1:6d} {2:6d} {3:6d}\n'.format(
                            bond_i,
                            bond_type_dict[temp],
                            bond.atom1 + offset,
                            bond.atom2 + offset))
                    bond_i += 1
                    bond_style.add(style)
                if len(bond_style) > 1:
                    logger.warn("More than one bond style found!")

                # angles
                logger.debug("        Writing angles...")
                for angle in mol_type.angleForceSet.itervalues():
                    atomtype1 = molecule._atoms[angle.atom1 - 1].bondtype
                    atomtype2 = molecule._atoms[angle.atom2 - 1].bondtype
                    atomtype3 = molecule._atoms[angle.atom3 - 1].bondtype

                    if isinstance(angle, Angle):
                        style = 'harmonic'
                        temp = AngleType(atomtype1, atomtype2, atomtype3,
                                angle.theta, angle.k)
                        # NOTE: k includes the factor of 0.5 for harmonic in LAMMPS
                        if temp not in angle_type_dict:
                            angle_type_dict[temp] = ang_type_i
                            angle_coeffs.append('{0:d} {1} {2:18.8f} {3:18.8f} # {4:2s}-{5:2s}-{6:2s}\n'.format(
                                    ang_type_i,
                                    style,
                                    0.5 * angle.k.in_units_of(self.ENERGY / self.RAD**2)._value,
                                    angle.theta.in_units_of(self.DEGREE)._value,
                                    atomtype1, atomtype2, atomtype3))
                            ang_type_i += 1
                    elif isinstance(angle, UreyBradleyAngle):
                        style = 'charmm'
                        temp = UreyBradleyAngleType(atomtype1, atomtype2, atomtype3,
                                angle.theta, angle.k, angle.r, angle.kUB)
                        # NOTE: k includes the factor of 0.5 for harmonic in LAMMPS
                        if temp not in angle_type_dict:
                            angle_type_dict[temp] = ang_type_i
                            angle_coeffs.append('{0:d} {1} {2:18.8f} {3:18.8f} {4:18.8f} {5:18.8f}\n'.format(
                                    ang_type_i,
                                    style,
                                    0.5 * angle.k.in_units_of(self.ENERGY / self.RAD**2)._value,
                                    angle.theta.in_units_of(self.DEGREE)._value,
                                    0.5 * angle.kUB.in_units_of(self.ENERGY / self.DIST**2)._value,
                                    angle.r.in_units_of(self.DIST)._value))
                            ang_type_i += 1
                    elif isinstance(angle, G96Angle):
                        style = 'cosine/squared'
                        temp = G96AngleType(atomtype1, atomtype2, atomtype3,
                                angle.theta, angle.k)
                        # NOTE: k includes the factor of 0.5 for harmonic in LAMMPS
                        if temp not in angle_type_dict:
                            angle_type_dict[temp] = ang_type_i
                            angle_coeffs.append('{0:d} {1} {2:18.8f} {3:18.8f}\n'.format(
                                    ang_type_i,
                                    style,
                                    0.5 * angle.k.in_units_of(self.ENERGY)._value,
                                    angle.theta.in_units_of(self.DEGREE)._value))
                            ang_type_i += 1
                    else:
                        logger.warn("Found unimplemented angle type for LAMMPS!")
                        continue

                    angle_list.append('{0:-6d} {1:6d} {2:6d} {3:6d} {4:6d}\n'.format(
                            angle_i,
                            angle_type_dict[temp],
                            angle.atom1 + offset,
                            angle.atom2 + offset,
                            angle.atom3 + offset))
                    angle_i += 1
                    angle_style.add(style)
                if len(angle_style) > 1:
                    logger.warn("More than one angle style found!")

                # dihedrals
                logger.debug("        Writing dihedrals...")

                for dihedral in mol_type.dihedralForceSet.itervalues():
                    atomtype1 = molecule._atoms[dihedral.atom1 - 1].bondtype
                    atomtype2 = molecule._atoms[dihedral.atom2 - 1].bondtype
                    atomtype3 = molecule._atoms[dihedral.atom3 - 1].bondtype
                    atomtype4 = molecule._atoms[dihedral.atom4 - 1].bondtype

                    if isinstance(dihedral, DihedralTrigDihedral):
                        coefficients = [dihedral.fc1, dihedral.fc2, dihedral.fc3,
                                dihedral.fc4, dihedral.fc5, dihedral.fc6]
                        if dihedral.improper:
                            found_nonzero = False
                            for n, coeff in enumerate(coefficients):
                                if coeff._value != 0.0:
                                    if found_nonzero == False:
                                        found_nonzero = True
                                    else:
                                        raise ValueError("Found more than one nonzero "
                                                "coefficient in improper trigonal dihedral!")
                                    style = 'charmm'
                                    temp = ProperPeriodicDihedralType(atomtype1, atomtype2,
                                            atomtype3, atomtype4,
                                            dihedral.phi, coeff, n + 1)
                                    if temp not in dihedral_type_dict:
                                        dihedral_type_dict[temp] = dih_type_i
                                        # NOTE: weighting factor assumed to be 0.0
                                        # May need to add some additional checks here in the future
                                        dihedral_coeffs.append('{0:d} {1} {2:18.8f} {3:18d} '
                                                    '{4:18d} {5:18.4f}\n'.format(
                                                dih_type_i, style,
                                                coeff.in_units_of(self.ENERGY)._value,
                                                n + 1,
                                                int(dihedral.phi.in_units_of(units.degrees)._value),
                                                0.0))
                                        dih_type_i += 1
                                    dihedral_list.append('{0:-6d} {1:6d} {2:6d} {3:6d} {4:6d} {5:6d}\n'.format(
                                            dihedral_i,
                                            dihedral_type_dict[temp],
                                            dihedral.atom1 + offset,
                                            dihedral.atom2 + offset,
                                            dihedral.atom3 + offset,
                                            dihedral.atom4 + offset))
                                    dihedral_i += 1
                                    dihedral_style.add(style)

                        else:
                            # NOTE: the following logic could instead default to printing
                            # out a series of charmm style dihedrals instead of attempting
                            # to write a multi/harmonic. I presume one multi/harmonic vs.
                            # multiple charmm dihedrals may be slightly (but probably
                            # negligibly) faster but if anyone has a better reason to do
                            # one or the other, please chime in!
                            rb_coeffs = ConvertDihedralFromDihedralTrigToRB(
                                    np.cos(dihedral.phi.in_units_of(units.radians)._value),
                                    dihedral.phi, dihedral.fc0, *coefficients)

                            # LAMMPS only supports multi/harmonic (Ryckaert-Bellemans)
                            # up to 5 coefficients.
                            if (dihedral.phi in [0*units.degrees, 180*units.degrees] and
                                    rb_coeffs[5]._value == rb_coeffs[6]._value == 0.0):
                                style = 'multi/harmonic'
                                temp = RBDihedralType(atomtype1, atomtype2,
                                        atomtype3, atomtype4,
                                        rb_coeffs[0],
                                        rb_coeffs[1],
                                        rb_coeffs[2],
                                        rb_coeffs[3],
                                        rb_coeffs[4],
                                        0.0*units.kilojoules_per_mole,
                                        0.0*units.kilojoules_per_mole)
                                if temp not in dihedral_type_dict:
                                    dihedral_type_dict[temp] = dih_type_i
                                    # multiple alternating powers by -1 for sign convention
                                    dihedral_coeffs.append('{0:d} {1} {2:18.8f} {3:18.8f} '
                                                '{4:18.8f} {5:18.8f} {6:18.8f}\n'.format(
                                            dih_type_i, style,
                                            rb_coeffs[0].in_units_of(self.ENERGY)._value,
                                            -rb_coeffs[1].in_units_of(self.ENERGY)._value,
                                            rb_coeffs[2].in_units_of(self.ENERGY)._value,
                                            -rb_coeffs[3].in_units_of(self.ENERGY)._value,
                                            rb_coeffs[4].in_units_of(self.ENERGY)._value))
                                    dih_type_i += 1
                                dihedral_list.append('{0:-6d} {1:6d} {2:6d} {3:6d} {4:6d} {5:6d}\n'.format(
                                        dihedral_i,
                                        dihedral_type_dict[temp],
                                        dihedral.atom1 + offset,
                                        dihedral.atom2 + offset,
                                        dihedral.atom3 + offset,
                                        dihedral.atom4 + offset))
                                dihedral_i += 1
                                dihedral_style.add(style)

                            # If the 6th and/or 7th coefficients are non-zero, we decompose
                            # the dihedral into multiple CHARMM style dihedrals.
                            else:
                                logger.warn("Found unsupported dihedral style.")
                                continue
                                """
                                for n, coeff in enumerate(coefficients):
                                    style = 'charmm'
                                    temp = ProperPeriodicDihedralType(atomtype1, atomtype2,
                                            atomtype3, atomtype4,
                                            dihedral.phi, coeff, n + 1)
                                    if temp not in dihedral_type_dict:
                                        dihedral_type_dict[temp] = dih_type_i
                                        # NOTE: weighting factor assumed to be 0.0
                                        # May need to add some additional checks here in the future
                                        dihedral_coeffs.append('{0:d} {1} {2:18.8f} {3:18d} '
                                                    '{4:18d} {5:18.4f}\n'.format(
                                                dih_type_i, style,
                                                coeff.in_units_of(self.ENERGY)._value,
                                                n + 1,
                                                int(dihedral.phi.in_units_of(units.degrees)._value),
                                                0.0))
                                        dih_type_i += 1
                                    dihedral_list.append('{0:-6d} {1:6d} {2:6d} {3:6d} {4:6d} {5:6d}\n'.format(
                                            i + j + 1,
                                            dihedral_type_dict[temp],
                                            dihedral.atom1 + offset,
                                            dihedral.atom2 + offset,
                                            dihedral.atom3 + offset,
                                            dihedral.atom4 + offset))
                                    dihedral_style.add(style)
                                """

                    elif isinstance(dihedral, ImproperHarmonicDihedral):
                        stlye = 'harmonic'
                        temp = ImproperHarmonicDihedralType(atomtype1, atomtype2,
                                atomtype3, atomtype4, dihedral.xi, dihedral.k)
                        if temp not in improper_type_dict:
                            improper_type_dict[temp] = imp_type_i
                            # NOTE: k includes the factor of 0.5 for harmonic in LAMMPS
                            improper_coeffs.append('{0:d} {1} {2:18.8f} {3:18.8f}\n'.format(
                                    imp_type_i, style,
                                    0.5 * dihedral.k.in_units_of(self.ENERGY / self.RAD**2)._value,
                                    dihedral.xi.in_units_of(self.DEGREE)._value))
                            imp_type_i += 1
                        improper_list.append('{0:-6d} {1:6d} {2:6d} {3:6d} {4:6d} {5:6d}\n'.format(
                                improper_i,
                                improper_type_dict[temp],
                                dihedral.atom1 + offset,
                                dihedral.atom2 + offset,
                                dihedral.atom3 + offset,
                                dihedral.atom4 + offset))
                        improper_i += 1
                        improper_style.add(style)
                    else:
                        raise Exception("InterMol expects all internally stored"
                                " dihedrals to be of types ImproperHarmonic"
                                " or DihedralTrig.")
                if len(dihedral_style) > 1:
                    logger.warn("More than one dihedral style found!")
                if len(improper_style) > 1:
                    logger.warn("More than one improper style found!")

                # atom specific information
                logger.debug("    Writing atoms...")

                for atom in molecule._atoms:
                    # type, mass and pair coeffs
                    if atom._atomtype[0] not in atom_type_dict:
                        atom_type_dict[atom._atomtype[0]] = a_type_i
                        mass_list.append('%d %8.4f # %s\n'
                                    % (a_type_i,
                                       atom._mass[0].in_units_of(self.MASS)._value,
                                       atom.bondtype))
                        pair_coeff_list.append('{0:d} {1:10.6f} {2:10.6f} # {3:s}\n'.format(
                                       a_type_i,
                                       atom._epsilon[0].in_units_of(self.ENERGY)._value,
                                       atom._sigma[0].in_units_of(self.DIST)._value,
                                       atom.bondtype))
                        a_type_i += 1

                    # box minima
                    x_coord = atom._position[0].in_units_of(self.DIST)._value
                    y_coord = atom._position[1].in_units_of(self.DIST)._value
                    z_coord = atom._position[2].in_units_of(self.DIST)._value
                    if x_coord < x_min:
                        x_min = x_coord
                    if y_coord < y_min:
                        y_min = y_coord
                    if z_coord < z_min:
                        z_min = z_coord

                    # atom
                    atom_list.append('{0:-6d} {1:-6d} {2:-6d} {3:5.8f} {4:8.5f} {5:8.5f} {6:8.5f}\n'.format(
                            atom.index + offset,
                            atom.residue_index,
                            atom_type_dict[atom._atomtype[0]],
                            atom._charge[0].in_units_of(self.CHARGE)._value,
                            x_coord,
                            y_coord,
                            z_coord))
                    # velocity
                    vel_list.append('{0:-6d} {1:8.4f} {2:8.4f} {3:8.4f}\n'.format(
                            atom.index + offset,
                            atom._velocity[0].in_units_of(self.VEL)._value,
                            atom._velocity[1].in_units_of(self.VEL)._value,
                            atom._velocity[2].in_units_of(self.VEL)._value))

                offset += len(molecule._atoms)

        # Write the actual data file.
        with open(data_file, 'w') as f:
            # front matter
            f.write(System._sys._name + '\n')
            f.write('\n')

            n_atoms = len(atom_list) - 3
            n_bonds = len(bond_list) - 3
            n_angles = len(angle_list) - 3
            n_dihedrals = len(dihedral_list) - 3
            n_impropers = len(improper_list) - 3

            n_atom_types = len(pair_coeff_list) - 3
            n_bond_types = len(bond_coeffs) - 3
            n_angle_types = len(angle_coeffs) - 3
            n_dihedral_types = len(dihedral_coeffs) - 3
            n_improper_types = len(improper_coeffs) - 3

            f.write('{0} atoms\n'.format(n_atoms))
            f.write('{0} bonds\n'.format(n_bonds))
            f.write('{0} angles\n'.format(n_angles))
            f.write('{0} dihedrals\n'.format(n_dihedrals))
            f.write('{0} impropers\n'.format(n_impropers))
            f.write('\n')

            f.write('{0} atom types\n'.format(n_atom_types))
            if n_bond_types > 0:
                f.write('{0} bond types\n'.format(n_bond_types))
            if n_angle_types > 0:
                f.write('{0} angle types\n'.format(n_angle_types))
            if n_dihedral_types > 0:
                f.write('{0} dihedral types\n'.format(n_dihedral_types))
            if n_improper_types > 0:
                f.write('{0} improper types\n'.format(n_improper_types))
            f.write('\n')

            # shifting of box dimensions
            f.write('{0:10.6f} {1:10.6f} xlo xhi\n'.format(
                    x_min,
                    x_min + System._sys.box_vector[0][0].in_units_of(self.DIST)._value))
            f.write('{0:10.6f} {1:10.6f} ylo yhi\n'.format(
                    y_min,
                    y_min + System._sys.box_vector[1][1].in_units_of(self.DIST)._value))
            f.write('{0:10.6f} {1:10.6f} zlo zhi\n'.format(
                    z_min,
                    z_min + System._sys.box_vector[2][2].in_units_of(self.DIST)._value))

            # masses
            for mass in mass_list:
                f.write(mass)

            # forcefield coefficients
            if len(pair_coeff_list) > 3:
                for pair in pair_coeff_list:
                    f.write(pair)
            if len(bond_coeffs) > 3:
                for bond in bond_coeffs:
                    f.write(bond)
            if len(angle_coeffs) > 3:
                for angle in angle_coeffs:
                    f.write(angle)
            if len(dihedral_coeffs) > 3:
                for dihedral in dihedral_coeffs:
                    f.write(dihedral)
            if len(improper_coeffs) > 3:
                for improper in improper_coeffs:
                    f.write(improper)

            # atoms and velocities
            for atom in atom_list:
                f.write(atom)
            for vel in vel_list:
                f.write(vel)

            # topology
            if len(bond_list) > 3:
                for bond in bond_list:
                    f.write(bond)
            if len(angle_list) > 3:
                for angle in angle_list:
                    f.write(angle)
            if len(dihedral_list) > 3:
                for dihedral in dihedral_list:
                    f.write(dihedral)
            if len(improper_list) > 3:
                for improper in improper_list:
                    f.write(improper)

        # Write the corresponding input file.
        basename = os.path.splitext(data_file)[0]
        input_filename = '{0}.input'.format(basename)
        with open(input_filename, 'w') as f:
            f.write('units {0}\n'.format(unit_set))
            f.write('atom_style full\n')  # TODO
            f.write('\n')

            f.write('dimension 3\n')  # TODO
            f.write('boundary p p p\n')  # TODO
            f.write('\n')

            # non-bonded
            f.write('pair_style lj/cut/coul/long 10.0 10.0\n')  # TODO: match mdp
            if System._sys.combination_rule == 3:
                f.write('pair_modify mix geometric\n')
            elif System._sys.combination_rule == 2:
                f.write('pair_modify mix arithmetic\n')
            else:
                logger.warn("Unsupported combination rule: {0}".format(
                        System._sys.combination_rule))
            f.write('kspace_style ewald 1.0e-5\n')  # TODO: match mdp
            f.write('\n')

            # bonded
            if len(bond_coeffs) > 3:
                f.write('bond_style hybrid {0}\n'.format(
                        " ".join(bond_style)))
            if len(angle_coeffs) > 3:
                f.write('angle_style hybrid {0}\n'.format(
                        " ".join(angle_style)))
            if len(dihedral_coeffs) > 3:
                f.write('dihedral_style hybrid {0}\n'.format(
                        " ".join(dihedral_style)))
            if len(improper_coeffs) > 3:
                f.write('improper_style hybrid {0}\n'.format(
                        " ".join(improper_style)))

            f.write('special_bonds lj {0} {1} {2} coul {3} {4} {5}\n'.format(
                    0.0,
                    0.0,
                    System._sys.lj_correction,
                    0.0,
                    0.0,
                    System._sys.coulomb_correction))
            f.write('\n')

            # read data
            f.write('read_data {0}\n'.format(os.path.basename(data_file)))
            f.write('\n')

            # output energies
            energy_terms = " ".join(['ebond',
                                     'eangle',
                                     'edihed',
                                     'eimp',
                                     'epair',
                                     'evdwl',
                                     'ecoul',
                                     'elong',
                                     'etail',
                                     'pe'])

            f.write('thermo_style custom {0}\n'.format(energy_terms))
            f.write('\n')

            # SHAKE constraints
            if len(shake_bond_types) > 0:
                f.write('fix fSHAKE all shake 1.0e-4 20 10 b')
                for btype in shake_bond_types:
                    f.write(' {0}'.format(bond_type_dict[btype]))
                if len(shake_angle_types) > 0:
                    f.write(' a')
                    for atype in shake_angle_types:
                        f.write(' {0}'.format(angle_type_dict[atype]))
                f.write('\n')

            f.write('run 0\n')
