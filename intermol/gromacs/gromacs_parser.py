from collections import OrderedDict
import logging
import os
from warnings import warn

import simtk.unit as units
from intermol.atom import Atom

from intermol.forces import *
from intermol.forces.forcefunctions import (create_kwds_from_entries,
                                            build_paramlist, build_unitvars)
from intermol.molecule import Molecule
from intermol.moleculetype import MoleculeType
from intermol.system import System
from grofile_parser import GromacsGroParser


logger = logging.getLogger('InterMolLog')


# TODO: FIGURE OUT ATOM AND MOLECULE NUMBERING


def load_gromacs(top_file, gro_file, include_dir=None, defines=None):
    """Load a set of GROMACS input files into a `System`.

    Args:
        top_file:
        gro_file:
        include_dir:
        defines:
    Returns:
        system:
    """
    parser = GromacsParser(top_file, gro_file, include_dir=None, defines=None)
    return parser.read()


def write_gromacs(top_file, gro_file, system):
    """Load a set of GROMACS input files into a `System`.

    Args:
        top_file:
        gro_file:
        include_dir:
        defines:
    Returns:
        system:
    """
    parser = GromacsParser(top_file, gro_file, system)
    return parser.write()




def _default_gromacs_include_dir():
    """Find the location where gromacs #include files are referenced from, by
    searching for (1) gromacs environment variables, (2) just using the default
    gromacs install location, /usr/local/gromacs/share/gromacs/top. """
    if 'GMXDATA' in os.environ:
        return os.path.join(os.environ['GMXDATA'], 'top')
    if 'GMXBIN' in os.environ:
        return os.path.abspath(os.path.join(
            os.environ['GMXBIN'], '..', 'share', 'gromacs', 'top'))

    return '/usr/local/gromacs/share/gromacs/top'


class GromacsParser(object):
    """
    A class containing methods required to read in a Gromacs(4.5.4) Topology File
    """


    # 'lookup_*' is the inverse dictionary typically used for writing
    gromacs_combination_rules = {
        '1': 'Multiply-C6C12',
        '2': 'Lorentz-Berthelot',
        '3': 'Multiply-Sigeps'
        }
    lookup_gromacs_combination_rules = dict(
        (v, k) for k, v in gromacs_combination_rules.items())

    gromacs_bonds = {
        '1': HarmonicBond,
        '2': G96Bond,
        '3': MorseBond,
        '4': CubicBond,
        '5': ConnectionBond,
        '6': HarmonicPotentialBond,
        '7': FeneBond
        }
    lookup_gromacs_bonds = dict((v, k) for k, v in gromacs_bonds.items())

    # create a BondType dictonary
    gromacs_bond_types = dict(
        (k, eval(v.__name__ + 'Type')) for k, v in gromacs_bonds.items())

    def canonical_bond(self, params, bond, direction='into'):
        """
        Args:
            params:
            bond:
            direction:
        Returns:
        """
        if direction == 'into':
            return bond, params
        else:  # currently, no bonds need to be de-canonicalized
            b_type = self.lookup_gromacs_bonds[bond.__class__]
            if b_type:
                return b_type, params
            else:
                warn("WriteError: found unsupported bond type {0}".format(
                        bond.__class__.__name__))

    gromacs_angles = {
        '1': HarmonicAngle,
        '2': CosineSquaredAngle,
        '3': CrossBondBondAngle,
        '4': CrossBondAngleAngle,
        '5': UreyBradleyAngle,
        '6': QuarticAngle
        }

    lookup_gromacs_angles = dict((v, k) for k, v in gromacs_angles.items())

    gromacs_angle_types = dict(
        (k, eval(v.__name__ + 'Type')) for k, v in gromacs_angles.items())

    def canonical_angle(self, params, angle, direction = 'into'):
        """
        Args:
            params:
            angle:
            direction:
        Returns:
        """
        if direction == 'into':
            return angle, params
        else:  # currently, no angles need to be de-canonicalized
            a_type = self.lookup_gromacs_angles[angle.__class__]
            if a_type:
                return a_type, params
            else:
                warn("WriteError: found unsupported angle type {0}".format(
                        angle.__class__.__name__))

    gromacs_dihedrals = {
        # TrigDihedrals are actually used for 1, 4, and 9.  Can't use lists as keys!
        '1': ProperPeriodicDihedral,
        '2': ImproperHarmonicDihedral,
        '3': RbDihedral,
        '4': ProperPeriodicDihedral,
        '5': FourierDihedral,
        '9': ProperPeriodicDihedral,
        'Trig': TrigDihedral,
        }

    lookup_gromacs_dihedrals = {
        TrigDihedral: 'Trig',
        ImproperHarmonicDihedral: '2',
        RbDihedral: '3',
        FourierDihedral: '5'
        }

    gromacs_dihedral_types = dict(
        (k, eval(v.__name__ + 'Type')) for k, v in gromacs_dihedrals.items())

    def canonical_dihedral(self, params, dihedral, direction='into'):
        """

        We can fit everything into two types of dihedrals - dihedral_trig, and
        improper harmonic dihedral trig is of the form fc0 + sum_i=1^6 fci
        (cos(nx-phi) proper dihedrals can be stored easily in this form, since
        they have only 1 n improper dihedrals can as well (flag as improper) RB
        can be stored as well, assuming phi = 0 or 180 Fourier can also be
        stored.  a full dihedral trig can be decomposied in to multiple proper
        dihedrals.

        Will need to handle multiple dihedrals little differently in that we
        will need to add multiple 9 dihedrals together into a single
        dihedral_trig, as long as they have the same phi angle (seems to be
        always the case).

        Args:
            params:
            dihedral:
            direction:
        Returns:
        """
        if direction == 'into':
            # if we are converting a type
            if 'Type' in dihedral.__name__:
                # Convert the dihedral parameters to the form we want to
                # actually store.
                converted_dihedral = dihedral  # Default.
                if dihedral == ProperPeriodicDihedralType:  # Proper dihedral.
                    convertfunc = convert_dihedral_from_proper_to_trig
                    converted_dihedral = TrigDihedralType
                elif dihedral == ImproperHarmonicDihedralType:
                    convertfunc = convert_nothing
                elif dihedral == RbDihedralType:
                    convertfunc = convert_dihedral_from_RB_to_trig
                    converted_dihedral = TrigDihedralType
                elif dihedral == FourierDihedralType:
                    convertfunc = convert_dihedral_from_fourier_to_trig
                    converted_dihedral = TrigDihedralType
                else:
                    # need some kind of exception here
                    pass
                # Now actually convert the dihedral.
                params = convertfunc(params)
            else:
                # Most will be TrigDihedral.
                converted_dihedral = TrigDihedral
                if dihedral == TrigDihedral:  # Already converted!
                    convertfunc = convert_nothing
                # Proper dihedral, still need to convert
                elif dihedral == ProperPeriodicDihedral:
                    convertfunc = convert_dihedral_from_proper_to_trig
                elif dihedral == ImproperHarmonicDihedral:
                    convertfunc = convert_nothing
                    converted_dihedral = dihedral
                elif dihedral == RbDihedral:
                    convertfunc = convert_dihedral_from_RB_to_trig
                elif dihedral == FourierDihedral:
                    convertfunc = convert_dihedral_from_fourier_to_trig
                else:
                    # need some kind of exception here
                    pass
                params = convertfunc(params)
            return converted_dihedral, params

        else:
            d_type = self.lookup_gromacs_dihedrals[dihedral.__class__]
            if not d_type:
                warn("WriteError: found unsupported dihedral type {0}".format(
                    dihedral.__class__.__name__))

            # Translate the dihedrals back to write them out.
            if isinstance(dihedral, TrigDihedral):
                # TODO: Exceptions for cases where tmpparams and paramlist don't
                # get assigned below.
                if dihedral.improper:
                    # Must be a improper dihedral. Print these out as improper
                    # dihedrals (d_type = 4).
                    d_type = '4'
                    paramlist = convert_dihedral_from_trig_to_proper(params)
                else:
                    # Print as a RB dihedral if the phi is 0.
                    if params['phi'].value_in_unit(units.degrees) in [0, 180]:
                        tmpparams = convert_dihedral_from_trig_to_RB(params)
                        if tmpparams['C6']._value == 0:
                            d_type = '3'
                        else:
                            # If C6 is not zero, then we have to print it out as
                            # multiple propers.
                            d_type = '9'
                    if d_type in ['9', 'Trig']:
                        # Print as proper dihedral. If one nonzero term, as a
                        # type 1, if multiple, type 9.
                        paramlist = convert_dihedral_from_trig_to_proper(params)
                        if len(paramlist) == 1:
                            d_type = '1'
                        else:
                            d_type = '9'
                    elif d_type == '3':
                        paramlist = [tmpparams]
            else:
                paramlist = [params]
            return d_type, paramlist

    paramlist = build_paramlist('gromacs')
    unitvars = build_unitvars('gromacs', paramlist)

    class _TopMoleculeType(object):
        """Inner class to store information about a molecule type."""
        def __init__(self):
            self.nrexcl = -1
            self.atoms = []
            self.bonds = []
            self.angles = []
            self.dihedrals = []
            self.exclusions = []
            self.pairs = []
            self.cmaps = []

    def __init__(self, top_file, gro_file, system=None, include_dir=None, defines=None):
        """
        Initializes a GromacsTopologyParse object which serves to read in a Gromacs
        topology into the abstract representation.

        Args:
            defines: Sets of default defines to use while parsing.
        """
        self.top_file = top_file
        self.gro_file = gro_file
        if not system:
            system = System()
        self.system = system

        if not include_dir is None:
            include_dir = _default_gromacs_include_dir()
        self._include_dirs = (os.path.dirname(top_file), include_dir)
        # Most of the gromacs water itp files for different forcefields,
        # unless the preprocessor #define FLEXIBLE is given, don't define
        # bonds between the water hydrogen and oxygens, but only give the
        # constraint distances and exclusions.
        self._defines = {'FLEXIBLE': True}
        if defines is not None:
            self._defines.update(defines)

    def read(self):
        """

        Return:
            system
        """
        self._current_directive = None
        self._if_stack = list()
        self._else_stack = list()
        self._molecule_types = OrderedDict()
        self._molecules = list()
        self._current_molecule_type = None
        self._current_molecule = None
        self._atomtypes = dict()
        self._bondtypes = dict()
        self._angletypes = dict()
        self._dihedraltypes = dict()
        self._implicittypes = dict()
        self._pairtypes = dict()
        self._cmaptypes = dict()

        # Parse the top_file into a set of plain text, intermediate
        # _TopMoleculeType objects.
        self._process_file(self.top_file)

        # Open the corresponding gro file and push all the information to the
        # InterMol system.
        self.gro = GromacsGroParser(self.gro_file)
        self.gro.read()
        self.system.box_vector = self.gro.box_vector
        self.system.n_atoms = self.gro.positions.shape[0]

        self.n_atoms_added = 0
        for mol_name, mol_count in self._molecules:
            if mol_name not in self._molecule_types:
                e = ValueError("Unknown molecule type: {0}".format(mol_name))
                logger.exception(e)
            # Grab the relevent plain text molecule type.
            top_moltype = self._molecule_types[mol_name]

            self._create_moleculetype(top_moltype, mol_name, mol_count)

        return self.system

    def write(self):
        """Write this topology in GROMACS file format.

        Args:
            filename: the name of the file to write out to
        """
        gro = GromacsGroParser(self.gro_file)
        gro.write(self.system)

        with open(self.top_file, 'w') as top:
            self._write_defaults(top)
            self._write_atomtypes(top)
            if self.system.nonbonded_types:
                self._write_nonbonded_types(top)

            self._write_moleculetypes(top)

            self._write_system(top)
            self._write_molecules(top)

    def _write_defaults(self, top):
        top.write('[ defaults ]\n')
        top.write('; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n')
        top.write('{0:6d} {1:6s} {2:6s} {3:8.4f} {4:8.4f}\n\n'.format(
                   self.system.nonbonded_function,
                   self.lookup_gromacs_combination_rules[self.system.combination_rule],
                   self.system.genpairs,
                   self.system.lj_correction,
                   self.system.coulomb_correction))

    def _write_atomtypes(self, top):
        top.write('[ atomtypes ]\n')
        top.write(';type, bondtype, atomic_number, mass, charge, ptype, sigma, epsilon\n')
        for atomtype in sorted(self.system.atomtypes.itervalues(), key=lambda x: x.atomtype):
            if atomtype.atomtype.isdigit():
                atomtype.atomtype = "LMP_{0}".format(atomtype.atomtype)
            if atomtype.bondtype.isdigit():
                atomtype.bondtype = "LMP_{0}".format(atomtype.bondtype)

            top.write('{0:<11s} {1:5s} {2:6d} {3:18.8f} {4:18.8f} {5:5s}'.format(
                    atomtype.atomtype,
                    atomtype.bondtype,
                    int(atomtype.atomic_number),
                    atomtype.mass.value_in_unit(units.atomic_mass_unit),
                    atomtype.charge.value_in_unit(units.elementary_charge),
                    atomtype.ptype))

            if self.system.combination_rule == 'Multiply-C6C12':
                top.write('{0:18.8e} {1:18.8e}\n'.format(
                    atomtype.sigma.value_in_unit(units.kilojoules_per_mole * units.nanometers**(6)),
                    atomtype.epsilon.value_in_unit(units.kilojoules_per_mole * units.nanometers**(12))))
            elif self.system.combination_rule in ['Lorentz-Berthelot','Multiply-Sigeps']:
                top.write('{0:18.8e} {1:18.8e}\n'.format(
                    atomtype.sigma.value_in_unit(units.nanometers),
                    atomtype.epsilon.value_in_unit(units.kilojoules_per_mole)))
        top.write('\n')

    def _write_nonbonded_types(self, top):
        top.write('[ nonbond_params ]\n')
        for nbtype in sorted(self.system.nonbonded_types.itervalues(), key=lambda x: (x.atom1, x.atom2)):
            top.write('{0:6s} {1:6s} {2:3d}'.format(
                    nbtype.atom1, nbtype.atom2, nbtype.type))
            if self.system.combination_rule == 'Multiply-C6C12':
                top.write('{3:18.8e} {4:18.8e}\n'.format(
                    nbtype.sigma.value_in_unit(units.kilojoules_per_mole * units.nanometers**(6)),
                    nbtype.epsilon.value_in_unit(units.kilojoules_per_mole * units.nanometers**(12))))
            elif self.system.combination_rule in ['Lorentz-Berthelot', 'Multiply-Sigeps']:
                top.write('{3:18.8e} {4:18.8e}\n'.format(
                    nbtype.sigma.value_in_unit(units.nanometers),
                    nbtype.epsilon.value_in_unit(units.kilojoules_per_mole)))
        top.write('\n')

    def _write_moleculetypes(self, top):
        for mol_name, mol_type in self.system.molecule_types.iteritems():
            self._current_molecule_type = mol_type
            top.write('[ moleculetype ]\n')
            # Gromacs can't handle spaces in the molecule name.
            printname = mol_name
            printname = printname.replace(' ', '_')
            printname = printname.replace('"', '')
            top.write('{0:s} {1:10d}\n\n'.format(printname, mol_type.nrexcl))

            self._write_atoms(top)

            # if moleculeType.pairForceSet:
            #     lines += self.write_pairs(moleculeType.pairForceSet)
            #
            # if moleculeType.bondForceSet and not moleculeType.settles:  # if settles, no bonds
            #     lines += self.write_bonds(moleculeType.bondForceSet)
            #
            # if moleculeType.angleForceSet and not moleculeType.settles: # if settles, no angles
            #     lines += self.write_angles(moleculeType.angleForceSet)
            #
            # if moleculeType.dihedralForceSet:
            #     lines += self.write_dihedrals(moleculeType.dihedralForceSet)
            #
            # if moleculeType.virtualForceSet:
            #     lines += self.write_virtuals(moleculeType.virtualForceSet)
            #
            # if moleculeType.settles:
            #     lines += self.write_settles(moleculeType.settles)
            #
            # if moleculeType.exclusions:
            #     # [ exclusions ]
            #     lines += self.write_exclusions(moleculeType.exclusions)

    def _write_system(self, top):
        top.write('[ system ]\n')
        top.write('{0}\n\n'.format(self.system.name))

    def _write_molecules(self, top):
        top.write('[ molecules ]\n')
        top.write('; Compound        nmols\n')
        for mol_name, mol_type in self.system.molecule_types.iteritems():
            n_molecules = len(mol_type.molecules)
            # The following lines are more 'chemical'.
            printname = mol_name
            printname = printname.replace(' ', '_')
            printname = printname.replace('"', '')
            top.write('{0:<15s} {1:8d}\n'.format(printname, n_molecules))

    def _write_atoms(self, top):
        top.write('[ atoms ]\n')
        top.write(';num, type, resnum, resname, atomname, cgnr, q, m\n')

        # Start iterating the set to get the first entry (somewhat kludgy...)
        for i, atom in enumerate(next(iter(self._current_molecule_type.molecules)).atoms):
            if atom.name.isdigit():  # LAMMPS atom names can have digits
                atom.name = "LMP_{0}".format(atom.name)
            if atom.atomtype[0].isdigit():
                atom.atomtype[0] = "LMP_{0}".format(atom.atomtype[0])

            top.write('{0:6d} {1:18s} {2:6d} {3:8s} {4:8s} {5:6d} '
                      '{6:18.8f} {7:18.8f}'.format(
                        i + 1,
                        atom.atomtype[0],
                        atom.residue_index,
                        atom.residue_name,
                        atom.name,
                        atom.cgnr,
                        atom.charge[0].value_in_unit(units.elementary_charge),
                        atom.mass[0].value_in_unit(units.atomic_mass_unit)))

            # Alternate states -- only one for now.
            if atom.atomtype.get(1):
                top.write('{0:18s} {1:18.8f} {2:18.8f}'.format(
                        atom.atomtype[1],
                        atom.charge[1].value_in_unit(units.elementary_charge),
                        atom.mass[1].value_in_unit(units.atomic_mass_unit)))
            top.write('\n')
        top.write('\n')

    def _create_moleculetype(self, top_moltype, mol_name, mol_count):
        # Create an intermol moleculetype.
        moltype = MoleculeType(mol_name)
        moltype.nrexcl = top_moltype.nrexcl
        self.system.add_molecule_type(moltype)
        self._current_molecule_type = moltype

        # Create all the intermol molecules of the current type.
        for n_mol in range(mol_count):
            self._create_molecule(top_moltype, mol_name)

    def _create_molecule(self, top_moltype, mol_name):
        molecule = Molecule(mol_name)
        self.system.add_molecule(molecule)
        self._current_molecule = molecule
        for atom in top_moltype.atoms:
            self._create_atom(atom)
        for bond in top_moltype.bonds:
            self._create_bond(bond)

    def _create_atom(self, temp_atom):
        index = self.n_atoms_added + 1
        atomtype = temp_atom[1]
        res_id = int(temp_atom[2])
        res_name = temp_atom[3]
        atom_name = temp_atom[4]
        cgnr = int(temp_atom[5])
        charge = float(temp_atom[6]) * units.elementary_charge
        mass = float(temp_atom[7]) * units.amu

        atom = Atom(index, atom_name, res_id, res_name)
        atom.cgnr = cgnr

        atom.atomtype = (0, atomtype)
        atom.charge = (0, charge)
        atom.mass = (0, mass)
        if len(temp_atom) == 11:
            atomtype = temp_atom[8]
            charge = float(temp_atom[9]) * units.elementary_charge
            mass = float(temp_atom[10]) * units.amu
            atom.atomtype = (1, atomtype)
            atom.charge = (1, charge)
            atom.mass = (1, mass)

        atom.position = self.gro.positions[self.n_atoms_added]

        for state, atomtype in atom.atomtype.iteritems():
            intermol_atomtype = self.system.atomtypes.get(atomtype)
            if not intermol_atomtype:
                logger.warn('A corresponding AtomType for {0} was not'
                            ' found.'.format(atom))
                continue
            atom.atomic_number = intermol_atomtype.atomic_number
            if not atom.bondtype:
                if intermol_atomtype.bondtype:
                    atom.bondtype = intermol_atomtype.bondtype
                else:
                    logger.warn("Suspicious bondtype parameter found for atom "
                                "{0}. Visually inspect before using.".format(atom))
            if not atom.mass.get(state):
                if intermol_atomtype.mass._value >= 0:
                    atom.mass = (state, intermol_atomtype.mass)
                else:
                    logger.warn("Suspicious mass parameter found for atom "
                                "{0}. Visually inspect before using.".format(atom))
            atom.sigma = (state, intermol_atomtype.sigma)
            atom.epsilon = (state, intermol_atomtype.epsilon)

        self._current_molecule.add_atom(atom)
        self.n_atoms_added += 1

    def _create_bond(self, bond):
        if len(bond) == 3:



    def _process_file(self, top_file):
        append = ''
        for line in open(top_file):
            if line.strip().endswith('\\'):
                append = '{0} {1}'.format(append, line[:line.rfind('\\')])
            else:
                self._process_line(top_file, '{0} {1}'.format(append, line))
                append = ''

    def _process_line(self, top_file, line):
        """Process one line from a file."""
        if ';' in line:
            line = line[:line.index(';')]
        stripped = line.strip()
        ignore = not all(self._if_stack)
        if stripped.startswith('*') or len(stripped) == 0:
            # A comment or empty line.
            return

        elif stripped.startswith('[') and not ignore:
            # The start of a category.
            if not stripped.endswith(']'):
                e =  ValueError('Illegal line in .top file: '+line)
                logger.exception(e)
            self._current_directive = stripped[1:-1].strip()
            logger.debug("Parsing {0}...".format(self._current_directive))

        elif stripped.startswith('#'):
            # A preprocessor command.
            fields = stripped.split()
            command = fields[0]
            if len(self._if_stack) != len(self._else_stack):
                e = RuntimeError('#if/#else stack out of sync')
                logger.exception(e)

            if command == '#include' and not ignore:
                # Locate the file to include
                name = stripped[len(command):].strip(' \t"<>')
                search_dirs = self._include_dirs+(os.path.dirname(top_file),)
                for sub_dir in search_dirs:
                    top_file = os.path.join(sub_dir, name)
                    if os.path.isfile(top_file):
                        # We found the file, so process it.
                        self._process_file(top_file)
                        break
                else:
                    e = ValueError('Could not locate #include file: '+name)
                    logger.exception(e)
            elif command == '#define' and not ignore:
                # Add a value to our list of defines.
                if len(fields) < 2:
                    e = ValueError('Illegal line in .top file: '+line)
                    logger.exception(e)
                name = fields[1]
                value_start = stripped.find(name, len(command))+len(name)+1
                value = line[value_start:].strip()
                self._defines[name] = value
            elif command == '#ifdef':
                # See whether this block should be ignored.
                if len(fields) < 2:
                    e = ValueError('Illegal line in .top file: '+line)
                    logger.exception(e)
                name = fields[1]
                self._if_stack.append(name in self._defines)
                self._else_stack.append(False)
            elif command == '#ifndef':
                # See whether this block should be ignored.
                if len(fields) < 2:
                    e = ValueError('Illegal line in .top file: '+line)
                    logger.exception(e)
                name = fields[1]
                self._if_stack.append(name not in self._defines)
                self._else_stack.append(False)
            elif command == '#endif':
                # Pop an entry off the if stack
                if len(self._if_stack) == 0:
                    e = ValueError('Unexpected line in .top file: '+line)
                    logger.exception(e)
                del(self._if_stack[-1])
                del(self._else_stack[-1])
            elif command == '#else':
                # Reverse the last entry on the if stack
                if len(self._if_stack) == 0:
                    e = ValueError('Unexpected line in .top file: '+line)
                    logger.exception(e)
                if self._else_stack[-1]:
                    e = ValueError('Unexpected line in .top file: '
                                   '#else has already been used ' + line)
                    logger.exception(e)
                self._if_stack[-1] = (not self._if_stack[-1])
                self._else_stack[-1] = True

        elif not ignore:
            # A line of data for the current category
            if self._current_directive is None:
                e = ValueError('Unexpected line in .top file: {0}'.format(line))
                logger.exception(e)
            if self._current_directive == 'defaults':
                self._process_defaults(line)
            elif self._current_directive == 'moleculetype':
                self._process_moleculetype(line)
            elif self._current_directive == 'molecules':
                self._process_molecule(line)
            elif self._current_directive == 'atoms':
                self._process_atom(line)
            elif self._current_directive == 'bonds':
                self._process_bond(line)
            elif self._current_directive == 'angles':
                self._process_angle(line)
            elif self._current_directive == 'dihedrals':
                self._process_dihedral(line)
            elif self._current_directive == 'exclusions':
                self._process_exclusion(line)
            elif self._current_directive == 'pairs':
                self._process_pair(line)
            elif self._current_directive == 'cmap':
                self._process_cmap(line)
            elif self._current_directive == 'atomtypes':
                self._process_atomtype(line)
            elif self._current_directive == 'bondtypes':
                self._process_bondtype(line)
            elif self._current_directive == 'angletypes':
                self._process_angletype(line)
            elif self._current_directive == 'dihedraltypes':
                self._process_dihedraltype(line)
            elif self._current_directive == 'implicit_genborn_params':
                self._process_implicittype(line)
            elif self._current_directive == 'pairtypes':
                self._process_pairtype(line)
            elif self._current_directive == 'cmaptypes':
                self._process_cmaptype(line)

    def _process_defaults(self, line):
        """Process the [ defaults ] line."""
        fields = line.split()
        if len(fields) < 4:
            self.too_few_fields(line)
        self.system.nonbonded_function = int(fields[0])
        self.system.combination_rule = self.gromacs_combination_rules[fields[1]]
        self.system.genpairs = fields[2]
        self.system.lj_correction = float(fields[3])
        self.system.coulomb_correction = float(fields[4])

    def _process_moleculetype(self, line):
        """Process a line in the [ moleculetypes ] category."""
        fields = line.split()
        if len(fields) < 1:
            self.too_few_fields(line)
        mol_type = self._TopMoleculeType()
        mol_type.nrexcl = int(fields[1])
        self._molecule_types[fields[0]] = mol_type
        self._current_molecule_type = mol_type

    def _process_molecule(self, line):
        """Process a line in the [ molecules ] category."""
        fields = line.split()
        if len(fields) < 2:
            self.too_few_fields(line)
        self._molecules.append((fields[0], int(fields[1])))

    def _process_atom(self, line):
        """Process a line in the [ atoms ] category."""
        if self._current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 5:
            self.too_few_fields(line)
        if len(fields) not in [8, 11]:
            self.invalid_line(line)
        self._current_molecule_type.atoms.append(fields)

    def _process_bond(self, line):
        """Process a line in the [ bonds ] category."""
        if self._current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 3:
            self.too_few_fields(line)
        self._current_molecule_type.bonds.append(fields)

    def _process_angle(self, line):
        """Process a line in the [ angles ] category."""
        if self._current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 4:
            self.too_few_fields(line)
        self._current_molecule_type.angles.append(fields)

    def _process_dihedral(self, line):
        """Process a line in the [ dihedrals ] category."""
        if self._current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 5:
            self.too_few_fields(line)
        self._current_molecule_type.dihedrals.append(fields)

    def _process_exclusion(self, line):
        """Process a line in the [ exclusions ] category."""
        if self._current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 2:
            self.too_few_fields(line)
        self._current_molecule_type.exclusions.append(fields)

    def _process_pair(self, line):
        """Process a line in the [ pairs ] category."""
        if self._current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 3:
            self.too_few_fields(line)
        self._current_molecule_type.pairs.append(fields)

    def _process_cmap(self, line):
        """Process a line in the [ cmaps ] category."""
        if self._current_molecule_type is None:
            self.directive_before_moleculetype('cmap')
        fields = line.split()
        if len(fields) < 6:
            self.too_few_fields(line)
        self._current_molecule_type.cmaps.append(fields)

    def _process_atomtype(self, line):
        """Process a line in the [ atomtypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            self.too_few_fields(line)
        if len(fields[3]) == 1:
            # Bonded type and atomic number are both missing.
            fields.insert(1, None)
            fields.insert(1, None)
        elif len(fields[4]) == 1 and len(fields[5]) > 1:
            if fields[1][0].isalpha():
                # Atomic number is missing.
                fields.insert(2, None)
            else:
                # Bonded type is missing.
                fields.insert(1, None)

        atomtype = fields[0]
        bondingtype = fields[1]
        if fields[2]:
            atomic_number = fields[2]
        else:
            atomic_number = -1
        mass = float(fields[3]) * units.amu
        charge = float(fields[4]) * units.elementary_charge
        ptype = fields[5]
        # Add correct units to the LJ parameters.
        if self.system.combination_rule == "Multiply-C6C12":
            lj_param1 = (float(fields[6]) *
                         units.kilojoules_per_mole * units.nanometers**(6))
            lj_param2 = (float(fields[7]) *
                         units.kilojoules_per_mole * units.nanometers**(12))
            AtomtypeClass = AtomCType
        elif self.system.combination_rule in ['Multiply-Sigeps', 'Lorentz-Berthelot']:
            lj_param1 = float(fields[6]) * units.nanometers           # sigma
            lj_param2 = float(fields[7]) * units.kilojoules_per_mole  # epsilon
            AtomtypeClass = AtomSigepsType
        else:
            e = ValueError("Unknown combination rule: {0}".format(
                self.system.combination_rule))
            logger.exception(e)
        new_atom_type = AtomtypeClass(atomtype, bondingtype, atomic_number,
                                      mass, charge, ptype, lj_param1, lj_param2)
        self.system.add_atomtype(new_atom_type)

    def _process_forcetype(self, forcename, line, n_atoms, gromacs_force_types,
                           canonical_force):
        """ """
        fields = line.split()

        numeric_forcetype = fields[n_atoms]
        gromacs_force_type = gromacs_force_types[numeric_forcetype]
        forcevars = fields[0:n_atoms]  # the bonding types

        kwds = create_kwds_from_entries(self.unitvars, self.paramlist, fields,
                                        gromacs_force_type, offset=n_atoms+1)
        CanonicalForceType, kwds = canonical_force(
            kwds, gromacs_force_type, direction='into')
        force_type = CanonicalForceType(*forcevars, **kwds)

        if not force_type:
            warn("{0} is not a supported {1} type".format(fields[2], forcename))
            return None
        else:
            return force_type

    def _process_bondtype(self, line):
        """Process a line in the [ bondtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            self.too_few_fields(line)

        bond_type = self._process_forcetype('bond', line, 2,
                self.gromacs_bond_types, self.canonical_bond)
        self._bondtypes[tuple(fields[:2])] = bond_type

    def _process_angletype(self, line):
        """Process a line in the [ angletypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            self.too_few_fields(line)
        angle_type = self._process_forcetype('angle', line, 3,
                self.gromacs_angle_types, self.canonical_angle)
        self._angletypes[tuple(fields[:3])] = angle_type

    def _process_dihedraltype(self, line):
        """Process a line in the [ dihedraltypes ] category."""
        fields = line.split()
        if len(fields) < 7:
            self.too_few_fields(line)


        # Some gromacs parameters don't include sufficient numbers of types.
        # Add some zeros (bit of a kludge).
        line += ' 0.0 0.0 0.0'
        fields = line.split()

        # Check whether they are using 2 or 4 atom types
        if fields[2].isdigit():
            d = 2
        elif fields[4].isdigit():
            d = 4

        dihedral_type = self._process_forcetype('dihedral', line, d,
                self.gromacs_dihedral_types, self.canonical_dihedral)

        # Still need a bit more information
        numeric_dihedraltype = fields[d]
        dihedral_type.improper = numeric_dihedraltype in ['2', '4']

        key = tuple(fields[:4])
        if numeric_dihedraltype == '9' and key in self._dihedraltypes:
            # There are multiple dihedrals defined for these atom types.
            self._dihedraltypes[key].append(dihedral_type)
        else:
            self._dihedraltypes[key] = [dihedral_type]

    def _process_implicittype(self, line):
        """Process a line in the [ implicit_genborn_params ] category."""
        fields = line.split()
        if len(fields) < 6:
            self.too_few_fields(line)
        self._implicittypes[fields[0]] = fields

    def _process_pairtype(self, line):
        """Process a line in the [ pairtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            self.too_few_fields(line)
        self._pairtypes[tuple(fields[:2])] = fields

    def _process_cmaptype(self, line):
        """Process a line in the [ cmaptypes ] category."""
        fields = line.split()
        if len(fields) < 8 or len(fields) < 8+int(fields[6])*int(fields[7]):
            self.too_few_fields(line)
        self._cmaptypes[tuple(fields[:5])] = fields

    def too_few_fields(self, line):
        e = ValueError('Too few fields in [ {0} ] line: {1}'.format(
                self._current_directive, line))
        logger.exception(e)

    def invalid_line(self, line):
        e = ValueError('Invalid format in [ {0} ] line: {1}'.format(
            self._current_directivedirective, line))
        logger.exception(e)

    def directive_before_moleculetype(self):
        e = ValueError('Found [ {0} ] directive before [ moleculetype ]'.format(
            self._current_directive))
        logger.exception(e)

if __name__ == "__main__":
    import pdb
    import intermol.tests

    tests_path = os.path.dirname(intermol.tests.__file__)
    top_file = os.path.join(tests_path, 'gromacs/unit_tests/dihedral9/dihedral9.top')
    gro_file = os.path.join(tests_path, 'gromacs/unit_tests/dihedral9/dihedral9.gro')
    gmx_system = load_gromacs(top_file, gro_file)

    print(gmx_system.n_atoms)
    top_file = 'converted_dihedral9.top'
    gro_file = 'converted_dihedral9.gro'
    write_gromacs(top_file, gro_file, gmx_system)

