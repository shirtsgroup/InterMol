"""LAMMPS extension for InterMol.

.. module:: lammps_parser
.. moduleauthor:: Christoph Klein <christoph.t.klein@me.com>
"""
import time
import os
from warnings import warn
import pdb

import numpy as np

import intermol.unit as units
from intermol.Atom import Atom
from intermol.Molecule import Molecule
from intermol.System import System
from intermol.Types import *
from intermol.Force import *


class LammpsParser(object):
    """A class containing methods to read and write LAMMPS files."""

    def __init__(self):
        """
        """

        self.box_vector = np.zeros(shape=(3, 3), dtype=float)

    def read_input(self, input_file):
        """Reads a LAMMPS input file.

        Args:
            input_file (str): Name of LAMMPS input file to read in.
        """
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
                'special_bonds': self.parse_special_bonds}

        with open(input_file, 'r') as input_lines:
            for line in input_lines:
                if line.strip():
                    keyword = line.split()[0]
                    if keyword in parsable_keywords:
                        parsable_keywords[keyword](line.split())

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
                'Dihedral Coeffs': self.parse_dihedral_coeffs}

        with open(data_file, 'r') as data_lines:
            self.molecule_name = next(data_lines).strip()
            # Currently only reading a single molecule/moleculeType
            # per LAMMPS file.
            self.current_mol = Molecule(self.molecule_name)
            System._sys.addMolecule(self.current_mol)
            self.current_mol_type = System._sys._molecules[self.molecule_name]
            self.current_mol_type.nrexcl = 3  # TODO: automate determination

            for line in data_lines:
                if line.strip():
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

        # Read atoms, velocities and connectivity information from data file.
        parsable_keywords = {'Atoms': self.parse_atoms,
                'Bonds': self.parse_bonds,
                'Angles': self.parse_angles,
                'Dihedrals': self.parse_dihedrals}

        with open(data_file, 'r') as data_lines:
            for line in data_lines:
                if line.strip():
                    keyword = line.strip()
                    if keyword in parsable_keywords:
                        parsable_keywords[keyword](data_lines)

    def parse_units(self, line):
        """ """
        assert(len(line) == 2), "Invalid units specified in input file."
        self.unit_set = line[1]

    def parse_atom_style(self, line):
        """
        Note:
            Assuming 'full' as default for everything else.
        """
        self.atom_style = line[1]
        if len(line) > 2:
            warn("Unsupported atom_style in input file.")

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
            warn("Hybrid pair styles not yet implemented.")
        elif line[1] == 'lj/cut/coul/long':
            self.pair_style.append(line[1])
            System._sys._nbFunc = 1

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
                System._sys._combinationRule = 3
            elif line[2] == 'arithmetic':
                System._sys._combinationRule = 2
            else:
                warn("Unsupported pair_modify mix argument in input file!")
        else:
            warn("Unsupported pair_modify style in input file!")

    def parse_bond_style(self, line):
        """ """
        self.bond_style = set()
        if len(line) == 2:
            self.bond_style.add(line[1])
        elif line[1] == 'hybrid':
            for style in line[2:]:
                self.bond_style.add(style)
        else:
            raise ValueError("Invalid bond_style in input file!")

    def parse_angle_style(self, line):
        """ """
        self.angle_style = []
        if len(line) == 2:
            self.angle_style.append(line[1])
        elif line[1] == 'hybrid':
            for style in line[2:]:
                self.angle_style.append(style)
        else:
            raise ValueError("Invalid angle_style in input file!")

    def parse_dihedral_style(self, line):
        """ """
        self.dihedral_style = []
        if len(line) == 2:
            self.dihedral_style.append(line[1])
            # TODO: correctly determine gen-pairs state
            if self.dihedral_style == 'opls':
                System._sys._genpairs = 'yes'
        elif line[1] == 'hybrid':
            self.dihedral_style = []
            for style in line[2:]:
                self.dihedral_style.append(style)
        else:
            raise ValueError("Invalid dihedral_style in input file!")

    def parse_improper_style(self, line):
        """ """
        self.improper_style = []
        if len(line) == 2:
            self.improper_style.append(line[1])
        elif line[1] == 'hybrid':
            self.improper_style = []
            for style in line[2:]:
                self.improper_style.append(style)
        else:
            raise ValueError("Invalid improper_style in input file!")

    def parse_special_bonds(self, line):
        """ """
        if 'lj/coul' in line:
            System._sys._ljCorrection = float(line[line.index('lj/coul') + 3])
            System._sys._coulombCorrection = float(line[line.index('lj/coul') + 3])
        elif 'lj' in line and 'coul' in line:
            System._sys._ljCorrection = float(line[line.index('lj') + 3])
            System._sys._coulombCorrection = float(line[line.index('coul') + 3])
        elif 'lj' in line:
            System._sys._ljCorrection = float(line[line.index('lj') + 3])
        elif 'coul' in line:
            System._sys._coulombCorrection = float(line[line.index('coul') + 3])
        else:
            warn("Unsupported special_bonds in input file.")

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
        System._sys.setBoxVector(self.box_vector * self.DIST)

    def parse_masses(self, data_lines):
        """Read masses from data file."""
        next(data_lines)  # toss out blank line
        self.mass_dict = dict()
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = line.split()
            self.mass_dict[int(fields[0])] = float(fields[1]) * self.MASS

    def parse_pair_coeffs(self, data_lines):
        """Read pair coefficients from data file."""
        next(data_lines)  # toss out blank line
        self.nb_types = dict()
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = [float(field) for field in line.split()]
            if len(self.pair_style) == 1:
                # TODO: lookup of type of pairstyle to determine format
                if System._sys._nbFunc == 1:
                    self.nb_types[int(fields[0])] = [fields[1] * self.ENERGY,
                                                     fields[2] * self.DIST]
                else:
                    warn("Unsupported pair coeff formatting in data file!")
            else:
                warn("Unsupported pair coeff formatting in data file!")

    def parse_bond_coeffs(self, data_lines):
        """Read bond coefficients from data file."""
        next(data_lines)  # toss out blank line
        self.bond_types = dict()
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = line.split()
            if len(self.bond_style) == 1:
                if 'harmonic' in self.bond_style:
                    self.bond_types[int(fields[0])] = [
                            'harmonic',
                            2 * float(fields[1]) * self.ENERGY / (self.DIST*self.DIST),
                            float(fields[2]) * self.DIST]
                if 'morse' in self.bond_style:
                    self.bond_types[int(fields[0])] = [
                            'morse',
                            float(fields[1]) * self.ENERGY,
                            float(fields[2]) * self.DIST**(-1),
                            float(fields[3]) * self.DIST]
            elif len(self.bond_style) > 1:
                style = fields[1]
                if style not in self.bond_style:
                    raise Exception("Bond type found in Bond Coeffs that "
                            "was not specified in bond_style: {0}".format(style))
                if style == 'harmonic':
                    self.bond_types[int(fields[0])] = [
                            style,
                            2 * float(fields[2]) * self.ENERGY / (self.DIST*self.DIST),
                            float(fields[3]) * self.DIST]
                if style == 'morse':
                    self.bond_types[int(fields[0])] = [
                            style,
                            float(fields[2]) * self.ENERGY,
                            float(fields[3]) * self.DIST**(-1),
                            float(fields[4]) * self.DIST]
            else:
                warn("Unsupported bond coeff formatting in data file!")

    def parse_angle_coeffs(self, data_lines):
        """Read angle coefficients from data file."""
        next(data_lines)  # toss out blank line
        self.angle_types = dict()
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = [float(field) for field in line.split()]
            if len(self.angle_style) == 1:
                if 'harmonic' in self.angle_style:
                    self.angle_types[int(fields[0])] = [
                            2 * fields[1] * self.ENERGY / (self.RAD*self.RAD),
                            fields[2] * self.DEGREE]
                else:
                    warn("Unsupported angle coeff formatting in data file!")
            else:
                warn("Unsupported angle coeff formatting in data file!")

    def parse_dihedral_coeffs(self, data_lines):
        """Read dihedral coefficients from data file."""
        next(data_lines)  # toss out blank line
        self.dihedral_types = dict()
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = [float(field) for field in line.split()]
            if len(self.dihedral_style) == 1:
                if 'opls' in self.dihedral_style:
                    self.dihedral_types[int(fields[0])] = [
                            fields[1] * self.ENERGY,
                            fields[2] * self.ENERGY,
                            fields[3] * self.ENERGY,
                            fields[4] * self.ENERGY]
                else:
                    warn("Unsupported dihedral coeff formatting in data file!")
            else:
                warn("Unsupported dihedral coeff formatting in data file!")

    def parse_atoms(self, data_lines):
        """Read atoms from data file."""
        next(data_lines)  # toss out blank line
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = line.split()

            if len(fields) in [7, 10]:
                if len(fields) == 10:
                    # TODO: store image flags?
                    pass
                new_atom_type = None
                if System._sys._combinationRule == 1:
                    warn("Combination rule '1' not yet implemented")
                elif System._sys._combinationRule in [2, 3]:
                    new_atom_type = AtomCR23Type(fields[2],    # atomtype
                            fields[2],                         # bondtype
                            -1,                                # Z
                            self.mass_dict[int(fields[2])],    # mass
                            float(fields[3]) * self.CHARGE,    # charge
                            'A',                               # ptype
                            self.nb_types[int(fields[2])][1],  # sigma
                            self.nb_types[int(fields[2])][0])  # epsilon

                System._sys._atomtypes.add(new_atom_type)

                atom = Atom(int(fields[0]),  # AtomNum
                        fields[2],           # atomName
                        int(fields[1]),      # resNum (molNum)
                        fields[1])           # resName (molNum)
                atom.setAtomType(0, fields[2])  # atomNum for LAMMPS
                atom.cgnr = 0  # TODO: look into alternatives
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
                        warn("Corresponding AtomType was not found. "
                                "Insert missing values yourself.")
            self.current_mol.addAtom(atom)

    def parse_bonds(self, data_lines):
        """Read bonds from data file."""
        next(data_lines)  # toss out blank line
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = [int(field) for field in line.split()]

            new_bond_force = None
            coeff_num = int(fields[1])
            if self.bond_types[coeff_num][0] == 'harmonic':
                r = self.bond_types[coeff_num][2]
                k = self.bond_types[coeff_num][1]
                new_bond_force = Bond(
                        fields[2], fields[3],
                        r, k)
            if self.bond_types[coeff_num][0] == 'morse':
                r = self.bond_types[coeff_num][3]
                D = self.bond_types[coeff_num][1]
                beta = self.bond_types[coeff_num][2]
                new_bond_force = Morse(
                        fields[2], fields[3],
                        r, D, beta)
            self.current_mol_type.bondForceSet.add(new_bond_force)

    def parse_angles(self, data_lines):
        """Read angles from data file."""
        next(data_lines)  # toss out blank line
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = [int(field) for field in line.split()]

            new_angle_force = None
            if len(self.angle_style) == 1:
                if self.angle_style[0] == 'harmonic':
                    theta = self.angle_types[int(fields[1])][1]
                    k = self.angle_types[int(fields[1])][0]
                    new_angle_force = Angle(
                            fields[2], fields[3], fields[4],
                            theta, k)
            else:
                warn("Hybrid angle styles not yet implemented")
            self.current_mol_type.angleForceSet.add(new_angle_force)

    def parse_dihedrals(self, data_lines):
        """Read dihedrals from data file."""
        next(data_lines)  # toss out blank line
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = [int(field) for field in line.split()]

            new_dihed_force = None
            if len(self.dihedral_style) == 1:
                if self.dihedral_style[0] == 'opls':
                    cs = [self.dihedral_types[fields[1]][0],
                          self.dihedral_types[fields[1]][1],
                          self.dihedral_types[fields[1]][2],
                          self.dihedral_types[fields[1]][3]]
                    new_dihed_force = FourierDihedral(
                            fields[2], fields[3], fields[4], fields[5],
                            cs[0], cs[1], cs[2], cs[3])

            else:
                warn("Hybrid dihedral styles not yet implemented")
            self.current_mol_type.dihedralForceSet.add(new_dihed_force)

    def write(self, data_file, unit_set='real', verbose=True):
        """Reads a LAMMPS data file.

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
        angle_type_dict = dict()
        ang_type_i = 1

        dihedral_style = set()
        dihedral_type_dict = dict()
        dih_type_i = 1

        improper_style = set()
        improper_type_dict = dict()
        imp_type_i = 1

        # read all atom specific and FF information
        for mol_type in System._sys._molecules.itervalues():
            if verbose:
                print "    Writing moleculetype {0}...".format(mol_type.name)
            # atom index offsets from 1 for each molecule
            offsets = list()
            for molecule in mol_type.moleculeSet:
                offsets.append(molecule._atoms[0].atomIndex - 1)

            molecule = mol_type.moleculeSet[0]
            atoms = molecule._atoms
            for i, offset in enumerate(offsets):
                if verbose:
                    start = time.time()
                    print "    Writing bonds..."
                    #import cProfile, pstats, StringIO
                    #pr = cProfile.Profile()
                    #pr.enable()

                for j, bond in enumerate(mol_type.bondForceSet.itervalues()):
                    #atom1 = mol_type.moleculeSet[0]._atoms[bond.atom1 - 1]
                    atom1 = molecule._atoms[bond.atom1 - 1]
                    atomtype1 = atom1.bondtype
                    #atom2 = mol_type.moleculeSet[0]._atoms[bond.atom2 - 1]
                    atom2 = molecule._atoms[bond.atom2 - 1]
                    atomtype2 = atom2.bondtype

                    if isinstance(bond, Bond):
                        style = 'harmonic'
                        temp = BondType(atomtype1, atomtype2,
                                bond.length, bond.k)
                        # NOTE: k includes the factor of 0.5 for harmonic in LAMMPS
                        if temp not in bond_type_dict:
                            bond_type_dict[temp] = b_type_i
                            bond_coeffs.append('{0:d} {1} {2:18.8f} {3:18.8f}\n'.format(
                                    b_type_i,
                                    style,
                                    0.5 * bond.k.in_units_of(self.ENERGY / (self.DIST*self.DIST))._value,
                                    bond.length.in_units_of(self.DIST)._value))
                            b_type_i += 1
                    elif isinstance(bond, Morse):
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
                        warn("Found unimplemented bond type for LAMMPS!")
                        continue

                    bond_list.append('{0:-6d} {1:6d} {2:6d} {3:6d}\n'.format(
                            i + j + 1,
                            bond_type_dict[temp],
                            bond.atom1 + offset,
                            bond.atom2 + offset))
                    bond_style.add(style)
                    if len(bond_style) > 1:
                        warn("More than one bond style found!")
                if verbose:
                    print "    Done. ({0:.2f} s)".format(time.time() - start)
                    #pr.disable()
                    #s = StringIO.StringIO()
                    #sortby = 'cumulative'
                    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
                    #ps.print_stats()
                    #print s.getvalue()

                # angle types
                if verbose:
                    start = time.time()
                    print "    Writing angles..."
                for j, angle in enumerate(mol_type.angleForceSet.itervalues()):
                    #atom1 = mol_type.moleculeSet[0]._atoms[angle.atom1 - 1]
                    atom1 = molecule._atoms[angle.atom1 - 1]
                    atomtype1 = atom1.bondtype
                    #atom2 = mol_type.moleculeSet[0]._atoms[angle.atom2 - 1]
                    atom2 = molecule._atoms[angle.atom2 - 1]
                    atomtype2 = atom2.bondtype
                    #atom3 = mol_type.moleculeSet[0]._atoms[angle.atom3 - 1]
                    atom3 = molecule._atoms[angle.atom3 - 1]
                    atomtype3 = atom3.bondtype

                    if isinstance(angle, Angle):
                        style = 'harmonic'
                        temp = AngleType(atomtype1, atomtype2, atomtype3,
                                angle.theta, angle.k)
                        # NOTE: k includes the factor of 0.5 for harmonic in LAMMPS
                        if temp not in angle_type_dict:
                            angle_type_dict[temp] = ang_type_i
                            angle_coeffs.append('{0:d} {1} {2:18.8f} {3:18.8f}\n'.format(
                                    ang_type_i,
                                    style,
                                    0.5 * angle.k.in_units_of(self.ENERGY / (self.RAD*self.RAD))._value,
                                    angle.theta.in_units_of(self.DEGREE)._value))
                            ang_type_i += 1
                    else:
                        warn("Found unimplemented angle type for LAMMPS!")
                        continue

                    angle_list.append('{0:-6d} {1:6d} {2:6d} {3:6d} {4:6d}\n'.format(
                            i + j + 1,
                            angle_type_dict[temp],
                            angle.atom1 + offset,
                            angle.atom2 + offset,
                            angle.atom3 + offset))

                    angle_style.add(style)
                    if len(angle_style) > 1:
                        warn("More than one angle style found!")
                if verbose:
                    print "    Done. ({0:.2f} s)".format(time.time() - start) 

                # dihedrals
                if verbose:
                    start = time.time()
                    print "    Writing dihedrals..."
                for j, dihedral in enumerate(mol_type.dihedralForceSet.itervalues()):
                    #atom1 = mol_type.moleculeSet[0]._atoms[dihedral.atom1 - 1]
                    atom1 = molecule._atoms[dihedral.atom1 - 1]
                    atomtype1 = atom1.bondtype
                    #atom2 = mol_type.moleculeSet[0]._atoms[dihedral.atom2 - 1]
                    atom2 = molecule._atoms[dihedral.atom2 - 1]
                    atomtype2 = atom2.bondtype
                    #atom3 = mol_type.moleculeSet[0]._atoms[dihedral.atom3 - 1]
                    atom3 = molecule._atoms[dihedral.atom3 - 1]
                    atomtype3 = atom3.bondtype
                    #atom4 = mol_type.moleculeSet[0]._atoms[dihedral.atom4 - 1]
                    atom4 = molecule._atoms[dihedral.atom4 - 1]
                    atomtype4 = atom4.bondtype

                    if isinstance(dihedral, FourierDihedral):
                        style = 'opls'
                        temp = FourierDihedralType(atomtype1, atomtype2,
                                atomtype3, atomtype4,
                                dihedral.c1, dihedral.c1,
                                dihedral.c2, dihedral.c3)
                        if temp not in dihedral_type_dict:
                            dihedral_type_dict[temp] = dih_type_i
                            dihedral_coeffs.append('{0:d} {1} {2:18.8f} {3:18.8f} {4:18.8f} {5:18.8f}\n'.format(
                                    dih_type_i,
                                    style,
                                    dihedral.c1.in_units_of(self.ENERGY)._value,
                                    dihedral.c2.in_units_of(self.ENERGY)._value,
                                    dihedral.c3.in_units_of(self.ENERGY)._value,
                                    dihedral.c4.in_units_of(self.ENERGY)._value))
                            dih_type_i += 1
                    else:
                        warn("Found unimplemented dihedral type for LAMMPS!")
                        continue

                    dihedral_list.append('{0:-6d} {1:6d} {2:6d} {3:6d} {4:6d} {5:6d}\n'.format(
                            i + j + 1,
                            dihedral_type_dict[temp],
                            dihedral.atom1 + offset,
                            dihedral.atom2 + offset,
                            dihedral.atom3 + offset,
                            dihedral.atom4 + offset))
                    dihedral_style.add(style)
                    if len(dihedral_style) > 1:
                        warn("More than one dihedral style found!")
                if verbose:
                    print "    Done. ({0:.2f} s)".format(time.time() - start)

            # atom specific information
            x_min = y_min = z_min = np.inf
            if verbose:
                start = time.time()
                print "    Writing atoms..."
            for molecule in mol_type.moleculeSet:
                for atom in molecule._atoms:
                    # type, mass and pair coeffs
                    if atom._atomtype[0] not in atom_type_dict:
                        atom_type_dict[atom._atomtype[0]] = a_type_i
                        mass_list.append('%d %8.4f\n'
                                    % (a_type_i,
                                       atom._mass[0].in_units_of(self.MASS)._value))
                        pair_coeff_list.append('{0:d} {1:8.4f} {2:8.4f}\n'.format(
                                       a_type_i,
                                       atom._epsilon[0].in_units_of(self.ENERGY)._value,
                                       atom._sigma[0].in_units_of(self.DIST)._value))
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
                            atom.atomIndex,
                            atom.residueIndex,
                            atom_type_dict[atom._atomtype[0]],
                            atom._charge[0].in_units_of(self.CHARGE)._value,
                            x_coord,
                            y_coord,
                            z_coord))
                    # velocity
                    vel_list.append('{0:-6d} {1:8.4f} {2:8.4f} {3:8.4f}\n'.format(
                            atom.atomIndex,
                            atom._velocity[0].in_units_of(self.VEL)._value,
                            atom._velocity[1].in_units_of(self.VEL)._value,
                            atom._velocity[2].in_units_of(self.VEL)._value))
            if verbose:
                print "    Done. ({0:.2f} s)".format(time.time() - start)

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
                    x_min + System._sys._boxVector[0][0].in_units_of(self.DIST)._value))
            f.write('{0:10.6f} {1:10.6f} ylo yhi\n'.format(
                    y_min,
                    y_min + System._sys._boxVector[1][1].in_units_of(self.DIST)._value))
            f.write('{0:10.6f} {1:10.6f} zlo zhi\n'.format(
                    z_min,
                    z_min + System._sys._boxVector[2][2].in_units_of(self.DIST)._value))

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
            f.write('pair_style lj/cut/coul/long 9.0 10.0\n')  # TODO: match mdp
            f.write('pair_modify mix geometric\n')  # TODO: match defaults
            f.write('kspace_style pppm 1.0e-5\n')  # TODO: match mdp
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
                    System._sys._ljCorrection,
                    0.0,
                    0.0,
                    System._sys._coulombCorrection))
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

            f.write('run 0\n')
