import os
import logging
import re

import parmed.unit as units
import numpy as np

from intermol.forces import *
import intermol.forces.forcefunctions as ff
from intermol.exceptions import (UnimplementedFunctional, UnsupportedFunctional,
                                 UnimplementedSetting, UnsupportedSetting,
                                 LammpsError, InterMolError)
from intermol.atom import Atom
from intermol.molecule import Molecule
from intermol.moleculetype import MoleculeType
from intermol.system import System

logger = logging.getLogger('InterMolLog')

ENGINE = 'lammps'


def load(in_file):
    """Load a LAMMPS input file into a `System`.

    Args:
        in_file:
        include_dir:
        defines:
    Returns:
        system:
    """
    parser = LammpsParser(in_file)
    return parser.read()


def save(in_file, system, unit_set='real',
                 nonbonded_style='pair_style lj/cut/coul/long 9.0 9.0\nkspace_style pppm 1e-6\n'):

    """write a `System` into LAMMPS input file.

    Args:
        in_file:
        system:
        unit_set:
        nonbonded_style: default is for a periodic system

    Returns:
        system:
    """
    parser = LammpsParser(in_file, system)
    return parser.write(unit_set=unit_set, nonbonded_style=nonbonded_style)


class LammpsParser(object):
    """A class containing methods to read and write LAMMPS files. """

    SCALE_INTO = 2.0
    SCALE_FROM = 0.5

    lammps_bonds = {
        'harmonic': HarmonicBond,
        'morse': MorseBond,
        'class2': QuarticBond,  # TODO: coefficients need special handling.
        'fene': FeneExpandableBond,  # TODO: need special handling for LJ terms
        'fene/expand': FeneExpandableBond,  # TODO: need special handling for LJ terms
        'quartic': QuarticBreakableBond,
        'nonlinear': NonlinearBond
        }
    lookup_lammps_bonds = {v: k for k, v in lammps_bonds.items()}
    # Add some non 1-to-1 mappings.
    lookup_lammps_bonds[HarmonicPotentialBond] = 'harmonic'
    lammps_bond_types = dict(
        (k, eval(v.__name__ + 'Type')) for k, v in lammps_bonds.items())

    def canonical_bond(self, params, bond, direction='into'):
        """Convert to/from the canonical form of this interaction. """
        # TODO: Gromacs says harmonic potential bonds do not have constraints or
        #       exclusions. Check that this logic is supported.
        if direction == 'into':
            canonical_force_scale = self.SCALE_INTO
        else:
            try:
                typename = self.lookup_lammps_bonds[bond.__class__]
            except KeyError:
                if bond.__class__.__name__ in ['FeneBond', 'ConnectionBond']:
                    raise UnimplementedFunctional(bond, ENGINE)
                else:
                    raise UnsupportedFunctional(bond, ENGINE)
            canonical_force_scale = self.SCALE_FROM

        if bond.__class__ in [HarmonicBond, HarmonicPotentialBond]:
            params['k'] *= canonical_force_scale

        if bond.__class__ == HarmonicPotentialBond:
            typename = 'harmonic'

        if direction == 'into':
            return bond, params
        else:
            return typename, [params]  # we expect a list

    lammps_angles = {
        'harmonic': HarmonicAngle,
        'cosine': CosineAngle,
        'cosine/squared': CosineSquaredAngle,
        'charmm': UreyBradleyAngle
        }
    lookup_lammps_angles = dict((v, k) for k, v in lammps_angles.items())
    lammps_angle_types = dict(
        (k, eval(v.__name__ + 'Type')) for k, v in lammps_angles.items())

    def canonical_angle(self, params, angle, direction):
        """Convert from the canonical form of this interaction. """
        if direction == 'into':
            canonical_force_scale = self.SCALE_INTO
            angletest = angle
        else:
            try:
                typename = self.lookup_lammps_angles[angle.__class__]
            except KeyError:
                raise UnsupportedFunctional(angle, ENGINE)
            angletest = angle.__class__
            canonical_force_scale = self.SCALE_FROM

        if angletest in [HarmonicAngle, CosineSquaredAngle, UreyBradleyAngle]:
            params['k'] *= canonical_force_scale

        if angletest == UreyBradleyAngle:
            params['kUB'] *= canonical_force_scale

        if direction == 'into':
            return angle, params
        else:
            return typename, [params]  # We expect a list

    lammps_dihedrals = {
        'opls': FourierDihedral,
        'multi/harmonic': RbDihedral,
        'charmm': ProperPeriodicDihedral,
        # not quite canonical form, but easily interconvertible
        }
    # Have to manually reverse dihedrals -- not unique.
    lookup_lammps_dihedrals = {
        TrigDihedral: 'Trig',
        RbDihedral: 'multi/harmonic',
        FourierDihedral: 'opls',
        ProperPeriodicDihedral: 'charmm'
        # not quite canonical form, but easily interconvertible
        }
    lammps_dihedral_types = dict(
        (k, eval(v.__name__ + 'Type')) for k, v in lammps_dihedrals.items())

    lammps_impropers = {
        'harmonic': ImproperHarmonicDihedral,
        'cvff': TrigDihedral,
        }
    lookup_lammps_impropers = dict((v, k) for k, v in lammps_impropers.items())
    lammps_improper_types = dict(
        (k, eval(v.__name__ + 'Type')) for k, v in lammps_impropers.items())

    def canonical_dihedral(self, params, dihedral, direction='into'):
        """Convert from the canonical form of this interaction. """

        if direction == 'into':
            converted_dihedral = dihedral  # Default
            if dihedral == ProperPeriodicDihedral:  # Proper dihedral
                convertfunc = convert_dihedral_from_proper_to_trig
                converted_dihedral = TrigDihedral
            elif dihedral == ImproperHarmonicDihedral:
                convertfunc = convert_nothing
            elif dihedral == RbDihedral:
                convertfunc = convert_dihedral_from_RB_to_trig
                converted_dihedral = TrigDihedral
            elif dihedral == FourierDihedral:
                convertfunc = convert_dihedral_from_fourier_to_trig
                converted_dihedral = TrigDihedral
                # Now actually convert the dihedral.
            params = convertfunc(params)

            # Adjust scaling conventions.
            canonical_force_scale = self.SCALE_INTO
            if converted_dihedral == ImproperHarmonicDihedralType:
                params['k'] *= canonical_force_scale

            return converted_dihedral, params

        else:
            dihedral_class = dihedral.__class__
            canonical_force_scale = self.SCALE_FROM
            if dihedral_class == TrigDihedral:
                if False: #dihedral.improper:
                    # TODO
                    d_type = '4'
                    paramlist = convert_dihedral_from_trig_to_proper(params)
                else:
                    if (params['phi'].value_in_unit(units.degrees) in [0, 180] and
                                params['fc5']._value == 0 and
                                params['fc6']._value == 0):
                        typename = 'multi/harmonic'
                        parameters = convert_dihedral_from_trig_to_RB(params)
                        paramlist = [parameters]
                    else:
                        # Print as proper dihedral. If one nonzero term, as a
                        # type 1, if multiple, type 9.
                        typename = 'charmm'
                        paramlist = convert_dihedral_from_trig_to_proper(params)

            elif dihedral_class == ImproperHarmonicDihedral:
                params['k'] *= canonical_force_scale
                typename = 'harmonic'
                paramlist = [params]
            else:
                raise UnsupportedFunctional(dihedral, ENGINE)
            return typename, paramlist

    def create_kwds_from_entries(self, entries, force_class, offset=0):
        return ff.create_kwds_from_entries(self.unitvars, self.paramlist,
                entries, force_class, offset=offset)

    def get_parameter_list_from_force(self, force):
        return ff.get_parameter_list_from_force(force, self.paramlist)

    def get_parameter_kwds_from_force(self, force):
        return ff.get_parameter_kwds_from_force(
                force, self.get_parameter_list_from_force, self.paramlist)

    def __init__(self, in_file, system=None):
        """
        """
        self.in_file = in_file
        if not system:
            system = System()
        self.system = system
        self.data_file = None

    def set_units(self, unit_set):
        """Set what unit set to use. """

        self.RAD = units.radians
        self.DEGREE = units.degrees
        self.MOLE = units.mole
        self.TEMP = units.kelvin
        if unit_set == 'real':
            self.DIST = units.angstroms
            self.VEL = units.angstroms / units.femtosecond
            self.ENERGY = units.kilocalorie / units.mole
            self.MASS = units.grams / units.mole
            self.CHARGE = units.elementary_charge
        elif unit_set == 'metal':
            self.DIST = units.angstroms
            self.VEL = units.angstroms / units.picosecond
            self.ENERGY = units.joule / units.coulomb * units.elementary_charge
            self.MASS = units.grams / units.mole
            self.CHARGE = units.elementary_charge
        elif unit_set == 'si':
            self.DIST = units.meters
            self.VEL = units.meters / units.second
            self.ENERGY = units.joules
            self.MASS = units.kilograms
            self.CHARGE = units.coulomb
        elif unit_set == 'cgs':
            self.DIST = units.centimeter
            self.VEL = units.centimeter / units.second
            self.ENERGY = units.erg
            self.MASS = units.grams
            self.CHARGE = np.sqrt(units.erg * units.centimeter)
        elif unit_set == 'micro':
            self.DIST = units.micrometers
            self.VEL = units.nanometers / units.nanosecond
            self.ENERGY = units.picogram * (
                units.micrometer / units.microsecond) ^ 2
            self.MASS = units.attograms
            self.CHARGE = units.elementary_charge
        elif unit_set == 'nano':
            self.DIST = units.nanometers
            self.VEL = units.nanometer / units.nanosecond
            self.ENERGY = units.attogram * (
                units.nanometer / units.nanosecond) ^ 2
            self.MASS = units.attograms
            self.CHARGE = units.elementary_charge
        elif unit_set == 'lj':
            self.DIST = units.dimensionless
            self.VEL = units.dimensionless
            self.ENERGY = units.dimensionless
            self.MASS = units.dimensionless
            self.CHARGE = units.dimensionless
            logger.warning("Using unit type lj: All values are dimensionless. "
                           "This is untested and will likely fail. "
                           "See LAMMPS doc for more.")
        elif unit_set == 'electron':
            self.DIST = units.bohr
            self.VEL = units.bohr / units.atu
            self.ENERGY = units.hartree
            self.MASS = units.amu
            self.CHARGE = units.elementary_charge
        else:
            raise LammpsError('Unsupported unit set specified: {0}'.format(unit_set))

        # Now create the dictionary of which units go in which order
        # for each command.  we need to pass 'self' so that we can
        # access the different unit sets, but the function unitvars is
        # not actually a member, so we have to do it in a nonstandard way.
        self.paramlist = ff.build_paramlist('lammps')
        self.unitvars = ff.build_unitvars('lammps', self.paramlist, dumself=self)

    def read(self):
        """Reads a LAMMPS input file and a data file specified within.

        Returns:
            system:
        """
        self.read_input()
        if self.data_file:
            self.read_data(self.data_file)
        else:
            raise LammpsError("No data file found in input script")
        return self.system

    def read_input(self):
        """Reads a LAMMPS input file.

        Args:
            input_file (str): Name of LAMMPS input file to read in.
        """
        self.input_dir = os.path.dirname(os.path.realpath(self.in_file))
        parsable_keywords = {
            'units': self.parse_units,
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
            'read_data': self.parse_read_data}

        defaults = [
            'units lj',
            'atom_style atomic',
            'dimension 3',
            'boundary p p p',
            'pair_style none',
            'kspace_style none',
            'pair_modify mix geometric shift no table 12 tabinner sqrt(2.0) tail no compute yes',
            'bond_style none',
            'angle_style none',
            'dihedral_style none',
            'improper_style none',
            'special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0 angle no dihedral no extra 0']

        keyword_defaults = {x.split()[0]: x for x in defaults}
        keyword_check = {x: False for x in keyword_defaults.keys()}

        with open(self.in_file, 'r') as input_lines:
            for line in input_lines:
                if line.strip():
                    keyword = line.split()[0]
                    if keyword in parsable_keywords:
                        parsable_keywords[keyword](line.split())
                        keyword_check[keyword] = True

        for key in keyword_check.keys():
            if not (keyword_check[key]):
                logger.warning('Keyword {0} not set, using LAMMPS default value {1}'.format(
                        key, " ".join(keyword_defaults[key].split()[1:])))
                parsable_keywords[key](keyword_defaults[key].split())

        self.set_units(self.unit_set)

        # Run some checks that cannot be evaluated until the full file has been parsed
        self.check_boundary_dimension_compatibility()

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
            self.molecule_name = next(data_lines).strip()
            # Currently only reading a single molecule/moleculeType
            # per LAMMPS file.
            self.current_mol_type = MoleculeType(self.molecule_name)
            self.current_mol_type.nrexcl = 3  # TODO: automate determination!
            # NOTE: nrexcl is a global in lammps and should probably be 
            # determined in parse_special_bonds
            self.system.add_molecule_type(self.current_mol_type)
            self.current_mol = Molecule(self.molecule_name)
            self.system.add_molecule(self.current_mol)

            for line in data_lines:
                if line.strip():
                    # Remove trailing comment
                    line = line.partition('#')[0]
                    # Catch all box dimensions.
                    if ('xlo' in line) and ('xhi' in line):
                        self.parse_box(line.split(), 0)
                    elif ('ylo' in line) and ('yhi' in line):
                        self.parse_box(line.split(), 1)
                    elif ('zlo' in line) and ('zhi' in line):
                        self.parse_box(line.split(), 2)
                    elif ('xy xz yz' in line):
                        self.parse_box(line.split(), 'tilt')
                    # Other headers.
                    else:
                        keyword = line.strip()
                        if keyword in parsable_keywords:
                            parsable_keywords[keyword](data_lines)

        # Read atoms, velocities and connectivity information from data file.
        parsable_keywords = {'Atoms': self.parse_atoms,
                             'Velocities': self.parse_velocities,
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

        # Indentify 1-2, 1-3, and 1-4 neighbors and create pair forces
        for mol_type in self.system.molecule_types.values():
            molecule = list(mol_type.molecules)[0]
            onetwo =   [set() for i in range(len(molecule.atoms) + 1)]
            onethree = [set() for i in range(len(molecule.atoms) + 1)]
            onefour =  [set() for i in range(len(molecule.atoms) + 1)]

            # 1-2 neighbors
            for bond in mol_type.bond_forces:
                onetwo[bond.atom1].add(bond.atom2)
                onetwo[bond.atom2].add(bond.atom1)

            for ai in [atom.index for atom in molecule.atoms]:
                # 1-3 neighbors
                for aj in onetwo[ai]:
                    for ak in onetwo[aj]:
                        if not ((ak == ai) or (ak in onetwo[ai])):
                            onethree[ai].add(ak)

                # 1-4 neighbors
                for aj in onethree[ai]:
                    for ak in onetwo[aj]:
                        if not ((ak == ai) or (ak in onetwo[ai]) or (ak in onethree[ai])):
                            onefour[ai].add(ak)

            # Generate 1-4 pairs (need to check nrexcl, lj/coulomb correction)
            for ai in [atom.index for atom in molecule.atoms]:
                for aj in onefour[ai]:
                    if aj >= ai:
                        mol_type.pair_forces.add(LjDefaultPair(ai, aj))

    def parse_units(self, line):
        """ """
        assert (len(line) == 2), "Invalid units specified in input file."
        self.unit_set = line[1]

    def parse_atom_style(self, line):
        """
        Note:
            Assuming 'full' as default for everything else.
        """
        self.atom_style = line[1]
        if len(line) > 2:
            raise UnimplementedSetting(line, ENGINE)

        # TODO: Add remaining atom_styles
        # http://lammps.sandia.gov/doc/atom_style.html
        # See issue:
        self.q_idx = None
        self.res_idx = None
        if self.atom_style == 'full':
            self.res_idx = 1
            self.type_idx = 2
            self.q_idx = 3
            self.pos_idx = (4, 5, 6)
        elif self.atom_style == 'molecular':
            self.res_idx = 1
            self.type_idx = 2
            self.pos_idx = (3, 4, 5)
        elif self.atom_style == 'charge':
            self.type_idx = 1
            self.q_idx = 2
            self.pos_idx = (3, 4, 5)
        elif self.atom_style in ('angle', 'atomic', 'bond'):
            self.type_idx = 1
            self.pos_idx = (2, 3, 4)

    def parse_dimension(self, line):
        """ """
        self.dimension = int(line[1])
        if self.dimension not in [2, 3]:
            raise LammpsError("Invalid dimension specified in input file "
                              "(must be 2 or 3).")

    def parse_boundary(self, line):
        """ """
        self.boundaries = [line[1], line[2], line[3]]

    def check_boundary_dimension_compatibility(self):
        """ """
        if len(self.boundaries) != self.dimension:
            raise LammpsError("Boundaries do not match specified dimension "
                              "in input file")

    def parse_pair_style(self, line):
        """ """
        self.pair_style = []
        if line[1] in ('lj/cut/coul/long', 'lj/cut', 'lj/cut/coul/cut'):
            self.pair_style.append(line[1])
            self.system.nonbonded_function = 1
        else:
            raise UnimplementedSetting(line, ENGINE)

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
                self.system.combination_rule = 'Multiply-Sigeps'
            elif line[2] == 'arithmetic':
                self.system.combination_rule = 'Lorentz-Berthelot'
            else:
                raise UnimplementedSetting(line, ENGINE)
        else:
            raise UnimplementedSetting(line, ENGINE)

    def parse_bonded_style(self, line):
        """ """
        style_set = set()
        if len(line) == 2:
            style_set.add(line[1])
        elif line[1] == 'hybrid':
            for style in line[2:]:
                style_set.add(style)
        else:
            raise LammpsError("Invalid style in input file: {}!".format(line))
        return style_set

    def parse_bond_style(self, line):
        """ """
        self.bond_style = self.parse_bonded_style(line)

    def parse_angle_style(self, line):
        """ """
        self.angle_style = self.parse_bonded_style(line)

    def parse_dihedral_style(self, line):
        """ """
        self.dihedral_style = self.parse_bonded_style(line)
        # TODO: correctly determine gen-pairs state
        if self.dihedral_style == 'opls':
            self.system.genpairs = 'yes'

    def parse_improper_style(self, line):
        """ """
        self.improper_style = self.parse_bonded_style(line)

    def parse_special_bonds(self, line):
        """ """
        if 'lj/coul' in line:
            self.system.lj_correction = float(line[line.index('lj/coul') + 3])
            self.system.coulomb_correction = float(
                line[line.index('lj/coul') + 3])
        elif 'lj' in line and 'coul' in line:
            self.system.lj_correction = float(line[line.index('lj') + 3])
            self.system.coulomb_correction = float(line[line.index('coul') + 3])
        elif 'lj' in line:
            self.system.lj_correction = float(line[line.index('lj') + 3])
        elif 'coul' in line:
            self.system.coulomb_correction = float(line[line.index('coul') + 3])
        else:
            raise UnimplementedSetting(line, ENGINE)

    def parse_read_data(self, line):
        """ """
        if len(line) == 2:
            self.data_file = os.path.join(self.input_dir, line[1])
        else:
            raise UnimplementedSetting(line, ENGINE)

    def parse_box(self, line, dim):
        """Read box information from data file.

        Args:
            line (str): Current line in input file.
            dim (int): Dimension specified in line.
        """
        if dim == 'tilt':
            fields = [(float(field) * self.DIST) for field in line[:3]]
            dirs = line[3:]
            for f,d in zip(fields,dirs):
                if d == 'xy':
                    self.system.box_vector[1,0] = f
                elif d == 'xz':
                    self.system.box_vector[2,0] = f
                elif d == 'yz':
                    self.system.box_vector[2,1] = f
        else:
            fields = [float(field) for field in line[:2]]
            box_length = fields[1] - fields[0]
            if box_length > 0:
                self.system.box_vector[dim, dim] = box_length * self.DIST
            else:
                raise LammpsError("Negative box length specified in data file.")

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
                if self.system.nonbonded_function == 1:
                    self.nb_types[int(fields[0])] = [fields[1] * self.ENERGY,
                                                     fields[2] * self.DIST]
                else:
                    raise UnimplementedSetting(line, ENGINE)
            else:
                raise UnimplementedFunctional(line, ENGINE)

    def parse_force_coeffs(self, data_lines, force_name, force_classes,
                           force_style, lammps_forces, canonical_force):
        """Read force coefficients from data file."""
        next(data_lines)  # toss out blank line

        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = line.partition('#')[0].split()

            warn = False
            if len(force_style) == 1:
                style = list(force_style)[0]  # awkward to have to translate to list to get the only member!
                if style == fields[1]:
                    field_offset = 2
                else:
                    if re.search('[a-zA-Z]+', fields[1]):
                        if style == 'none':
                            style = fields[1]
                            field_offset = 2
                        else:
                            warn = True
                    else:
                        field_offset = 1

            elif len(force_style) > 1:
                style = fields[1]
                field_offset = 2
                if style not in force_style:
                    warn = True
            else:
                raise LammpsError("No entries found in '%s_style'." % (force_name))

            if warn:
                logger.warning('{0} type found in {1} Coeffs that was not '
                               'specified in {2}_style: {3}'.format(force_name, force_name, force_name, style))

            # what internal force correspond to this style
            force_class = lammps_forces[style]

            # Get the parameters from the line and translate into keywords
            kwds = self.create_kwds_from_entries(fields, force_class,
                                                 offset=field_offset)
            # translate the force into canonical form
            force_class, kwds = canonical_force(kwds, force_class,
                                                direction='into')
            # add to the dictionary of this force term
            force_classes[int(fields[0])] = [force_class, kwds]

    def parse_bond_coeffs(self, data_lines):

        self.bond_classes = dict()
        self.parse_force_coeffs(data_lines, "Bond",
                                self.bond_classes, self.bond_style,
                                self.lammps_bonds, self.canonical_bond)

    def parse_angle_coeffs(self, data_lines):

        self.angle_classes = dict()
        self.parse_force_coeffs(data_lines, "Angle",
                                self.angle_classes, self.angle_style,
                                self.lammps_angles, self.canonical_angle)

    def parse_dihedral_coeffs(self, data_lines):

        self.dihedral_classes = dict()
        self.parse_force_coeffs(data_lines, "Dihedral",
                                self.dihedral_classes, self.dihedral_style,
                                self.lammps_dihedrals, self.canonical_dihedral)

    def parse_improper_coeffs(self, data_lines):

        self.improper_classes = dict()
        self.parse_force_coeffs(data_lines, "Improper",
                                self.improper_classes, self.improper_style,
                                self.lammps_impropers, self.canonical_dihedral)

    def parse_atoms(self, data_lines):
        """Read atoms from data file."""
        next(data_lines)  # toss out blank line

        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = line.partition('#')[0].split()

            if len(fields) == 10:
                # TODO: store image flags?
                pass
            new_atom_type = None
            if self.system.combination_rule == "Multiply-C6C12":
                logger.warning("Combination rule 'Multiply-C6C12' not yet implemented")
            elif self.system.combination_rule in ['Multiply-Sigeps',
                                                  'Lorentz-Berthelot']:
                atomtype = 'lmp_{:03d}'.format(int(fields[self.type_idx]))
                bondtype = atomtype
                new_atom_type = AtomSigepsType(
                    atomtype,  # atomtype
                    bondtype,  # bondtype
                    -1,  # atomic_number
                    self.mass_dict[int(fields[self.type_idx])],  # mass
                    0 * self.CHARGE,  # charge (0 for atomtype)
                    'A',  # ptype
                    self.nb_types[int(fields[self.type_idx])][1],  # sigma
                    self.nb_types[int(fields[self.type_idx])][0])  # epsilon

            self.system.add_atomtype(new_atom_type)

            atom = Atom(int(fields[0]),  # index
                        'A{:x}'.format(int(fields[0])),  # name (A + idx in hex)
                        1,  # residue_index (lammps doesn't track this)
                        'R{:02d}'.format(1))  # residue_name (R + moltype num)
            atom.atomtype = (0, atomtype)
            atom.atomic_number = 0  #TODO: this must be defined for Desmond output; can we get this from LAMMPS?
            atom.cgnr = 0  # TODO: look into alternatives
            atom.mass = (0, self.mass_dict[int(fields[self.type_idx])])
            atom.epsilon = (0, self.nb_types[int(fields[self.type_idx])][0])
            atom.sigma = (0, self.nb_types[int(fields[self.type_idx])][1])
            atom.bondingtype = bondtype

            if self.q_idx:
                atom.charge = (0, float(fields[self.q_idx]) * self.CHARGE)
            else:
                atom.charge = (0, 0.0 * self.CHARGE)

            atom.position = [float(fields[self.pos_idx[0]]) * self.DIST,
                             float(fields[self.pos_idx[1]]) * self.DIST,
                             float(fields[self.pos_idx[2]]) * self.DIST]

            self.current_mol.add_atom(atom)
            
    def parse_velocities(self, data_lines):
        """ """
        next(data_lines)
        atoms = self.current_mol.atoms
        vel_dict = dict()
        for line in data_lines:
            if not line.strip():
                break
            fields = [field for field in line.partition('#')[0].split()]
            vel_dict[int(fields[0])] = fields[1:4]
        for atom in atoms:
            atom._velocity = [float(vel) * self.VEL for vel in
                              vel_dict[atom.index]]

    def parse_force(self, data_lines, force_classes, forceSet, n=0):
        """Read bonds, angles, dihedrals, impropers from data file."""
        next(data_lines)  # toss out blank line
        for line in data_lines:
            if not line.strip():
                break  # found another blank line
            fields = [int(field) for field in line.partition('#')[0].split()]

            coeff_num = fields[1]
            atom_nums = fields[2:n + 2]
            paraminfo = force_classes[coeff_num]
            kwds = paraminfo[1]
            new_force = paraminfo[0](*atom_nums, **kwds)
            forceSet.add(new_force)

    def parse_bonds(self, data_lines):
        self.parse_force(data_lines, self.bond_classes,
                         self.current_mol_type.bond_forces, n=2)

    def parse_angles(self, data_lines):
        self.parse_force(data_lines, self.angle_classes,
                         self.current_mol_type.angle_forces, n=3)

    def parse_dihedrals(self, data_lines):
        self.parse_force(data_lines, self.dihedral_classes,
                         self.current_mol_type.dihedral_forces, n=4)

    def parse_impropers(self, data_lines):
        self.parse_force(data_lines, self.improper_classes,
                         self.current_mol_type.dihedral_forces, n=4)

    def get_force_atoms(self, force, forceclass):
        """Return the atoms involved in a force. """
        if forceclass in ['Bond', 'Pair']:
            return [force.atom1, force.atom2]
        elif forceclass in ['Angle']:
            return [force.atom1, force.atom2, force.atom3]
        elif forceclass in ['Dihedral', 'Improper']:
            return [force.atom1, force.atom2, force.atom3, force.atom4]
        else:
            logger.warning("No interaction type %s defined!" % (forceclass))

    def get_force_bondingtypes(self, force, forceclass):
        """Return the atoms involved in a force. """
        if forceclass in ['Bond', 'Pair']:
            return [force.bondingtype1, force.bondingtype2]
        elif forceclass in ['Angle']:
            return [force.bondingtype1, force.bondingtype2, force.bondingtype3]
        elif forceclass in ['Dihedral', 'Improper']:
            return [force.bondingtype1, force.bondingtype2, force.bondingtype3,
                    force.bondingtype4]
        else:
            logger.warning("No interaction type %s defined!" % (forceclass))

    def write_forces(self, forces, offset, force_name,
                     lookup_lammps_force, lammps_force_types, canonical_force):
        """The general force writing function.

        Currently supports bonds, angles, dihedrals, impropers.
        """
        logger.debug("        Writing {0:s}s...".format(force_name))
        force_list = self.force_dict[force_name]
        force_count = len(force_list)

        coeff_name = '{0} Coeffs'.format(force_name)
        coeff_list = self.force_dict[coeff_name]

        numeric_coeff = self.numeric_coeff[coeff_name]
        type_count = len(numeric_coeff) + 1

        style_set = self.style_dict[force_name]

        for force in forces:
            atom_indices = self.get_force_atoms(force, force_name)
            atom_bondingtypes = self.get_force_bondingtypes(force, force_name)
            try:
                lookup_lammps_force[force.__class__]
            except KeyError:
                logger.warning("Found unimplemented {0} type {1} for LAMMPS!".format(
                    force_name, force.__class__.__name__))

            # Get the parameters of the force.
            kwds = self.get_parameter_kwds_from_force(force)

            # Convert keywords from canonical form.
            style, kwdslist = canonical_force(kwds, force, direction='from')
            force_type = lammps_force_types[style]
            style_set.add(style)

            # A single force can produce multiple forces.
            for kwds in kwdslist:
                temp_force_type = force_type(*atom_bondingtypes, **kwds)

                # New type found. Write out the force coefficients.
                # TODO: Check by content would be nice but may lead to the
                # equality issues again.
                if temp_force_type not in numeric_coeff:
                    # Get the numerical type for this interaction.
                    numeric_coeff[temp_force_type] = type_count
                    line = '{0:d} {1}'.format(type_count, style)
                    type_count += 1

                    # Generate the list of parameters for this force in the
                    # order they appear in the file format.
                    params = self.get_parameter_list_from_force(temp_force_type)

                    # Generate the units for this force.
                    u = self.unitvars[force_type.__name__]
                    for i, p in enumerate(params):
                        if p.unit == units.dimensionless and isinstance(p._value, int):
                            # LAMMPS expects an integer.
                            line += "%10d" % (p.value_in_unit(u[i]))
                        elif style == 'charmm' and p.unit == units.degrees:
                            if force_name == 'Dihedral':
                                # LAMMPS loves enforcing unnecessary integers.
                                line += "%10d" % (p.value_in_unit(u[i]))
                            else:
                                line += "%18.8e" % (p.value_in_unit(u[i]))
                        else:
                            line += "%18.8e" % (p.value_in_unit(u[i]))
                    line += '\n'
                    coeff_list.append(line)

                # Write out the force entry.
                line = '{0:-6d} {1:6d}'.format(force_count, numeric_coeff[temp_force_type])
                for atom in atom_indices:
                    line += ' {0:d}'.format(atom + offset)
                line += '\n'
                force_list.append(line)
                force_count += 1

    def write_bonds(self, mol_type, offset):
        bonds = sorted(mol_type.bond_forces, key=lambda x: (x.atom1, x.atom2))
        self.write_forces(bonds, offset, "Bond",
                          self.lookup_lammps_bonds,
                          self.lammps_bond_types,
                          self.canonical_bond)


    def write_angles(self, mol_type, offset):
        angles = sorted(mol_type.angle_forces, key=lambda x: (x.atom1, x.atom2, x.atom3))
        return self.write_forces(angles, offset, "Angle",
                                 self.lookup_lammps_angles,
                                 self.lammps_angle_types,
                                 self.canonical_angle)

    def write_dihedrals(self, mol_type, offset):
        """Separate dihedrals from impropers. """
        dihedral_forces = {force for force in mol_type.dihedral_forces
                           if force.__class__ != ImproperHarmonicDihedral}
        dihedrals = sorted(dihedral_forces, key=lambda x: (x.atom1, x.atom2, x.atom3, x.atom4))
        return self.write_forces(dihedrals, offset, "Dihedral",
                                 self.lookup_lammps_dihedrals,
                                 self.lammps_dihedral_types,
                                 self.canonical_dihedral)

    def write_impropers(self, mol_type, offset):
        """Separate dihedrals from impropers. """
        improper_forces = {force for force in mol_type.dihedral_forces
                           if force.__class__ == ImproperHarmonicDihedral}
        impropers = sorted(improper_forces, key=lambda x: (x.atom1, x.atom2, x.atom3, x.atom4))
        return self.write_forces(impropers, offset, "Improper",
                                 self.lookup_lammps_impropers,
                                 self.lammps_improper_types,
                                 self.canonical_dihedral)

    def write_virtuals(self, mol_type, offset):
        if len(mol_type.virtual_forces) > 0:
            logger.warning('Virtuals not currently supported: will need to be '
                           'implemeneted from shake and rigid')

    def write(self, unit_set='real',nonbonded_style=None):
        """Writes a LAMMPS data and corresponding input file.

        Args:
            data_file (str): Name of LAMMPS data file to write to.
            unit_set (str): LAMMPS unit set for output file.
        """
        self.data_file = os.path.splitext(self.in_file)[0] + '.lmp'
        self.set_units(unit_set)

        # Containers for lines which are ultimately written to output files.
        mass_list = list()
        mass_list.append('\nMasses\n\n')

        pair_coeffs = list()

        atom_list = list()
        atom_list.append('\nAtoms\n\n')

        vel_list = list()
        vel_list.append('\nVelocities\n\n')

        # Dicts for type information.
        atom_type_dict = dict()  # str_type:int_type
        a_type_i = 1  # counter for atom types

        # Dicts to store the final outputs.
        self.force_dict = {'Bond': ['\nBonds\n\n'],
                           'Bond Coeffs': ['\nBond Coeffs\n\n'],
                           'Angle': ['\nAngles\n\n'],
                           'Angle Coeffs': ['\nAngle Coeffs\n\n'],
                           'Dihedral': ['\nDihedrals\n\n'],
                           'Dihedral Coeffs': ['\nDihedral Coeffs\n\n'],
                           'Improper': ['\nImpropers\n\n'],
                           'Improper Coeffs': ['\nImproper Coeffs\n\n']
                           }

        # Dicts to store the numeric values for each type (force:int).
        self.numeric_coeff = {'Bond Coeffs': {},
                              'Angle Coeffs': {},
                              'Dihedral Coeffs': {},
                              'Improper Coeffs': {}
                              }

        self.style_dict = {'Bond': set(),
                           'Angle': set(),
                           'Dihedral': set(),
                           'Improper': set()
                           }

        # Read all atom specific and FF information.
        offset = 0
        for mol_name, mol_type in self.system.molecule_types.items():
            logger.debug(
                "    Writing moleculetype {0}...".format(mol_name))

            # OrderedSet isn't indexable so get the first molecule by iterating.
            molecule = next(iter(mol_type.molecules))
            atoms_per_molecule = len(molecule.atoms)

            for molecule in mol_type.molecules:
                # Atom index offsets from 1 for each molecule.
                self.write_bonds(mol_type, offset)
                self.write_angles(mol_type, offset)
                self.write_dihedrals(mol_type, offset)
                self.write_impropers(mol_type, offset)
                # Only issues warning now.
                self.write_virtuals(mol_type, offset)
                offset += atoms_per_molecule

            # Atom specific information.
            x_min = y_min = z_min = np.inf
            logger.debug("    Writing atoms...")
            atom_charges = False
            for molecule in mol_type.molecules:
                for atom in molecule.atoms:
                    # Type, mass and pair coeffs.
                    if atom.atomtype[0] not in atom_type_dict:
                        atom_type_dict[atom.atomtype[0]] = a_type_i
                        mass_list.append('{0:d} {1:11.7f}\n'.format(
                                a_type_i,
                                atom.mass[0].value_in_unit(self.MASS)))
                        pair_coeffs.append('pair_coeff {0:d} {0:d} {1:11.7f} {2:11.7f}\n'.format(
                                a_type_i,
                                atom.epsilon[0].value_in_unit(self.ENERGY),
                                atom.sigma[0].value_in_unit(self.DIST)))
                        a_type_i += 1

                    # Box minima.
                    x_coord = atom.position[0].value_in_unit(self.DIST)
                    y_coord = atom.position[1].value_in_unit(self.DIST)
                    z_coord = atom.position[2].value_in_unit(self.DIST)
                    if x_coord < x_min:
                        x_min = x_coord
                    if y_coord < y_min:
                        y_min = y_coord
                    if z_coord < z_min:
                        z_min = z_coord

                    atom_list.append('{0:-6d} {1:-6d} {2:-6d} {3:5.8f} {4:12.7f} {5:12.7f} {6:12.7f}\n'.format(
                            atom.index,
                            atom.residue_index,
                            atom_type_dict[atom.atomtype[0]],
                            atom.charge[0].value_in_unit(self.CHARGE),
                            x_coord,
                            y_coord,
                            z_coord))

                    if atom.charge[0]._value != 0:
                        atom_charges = True
                    if np.any(atom.velocity):
                        vel_list.append(
                            '{0:-6d} {1:11.7f} {2:11.7f} {3:11.7f}\n'.format(
                                atom.index,
                                atom.velocity[0].value_in_unit(self.VEL),
                                atom.velocity[1].value_in_unit(self.VEL),
                                atom.velocity[2].value_in_unit(self.VEL)))
                    else:
                        vel_list.append(
                            '{0:-6d} {1:11.7f} {2:11.7f} {3:11.7f}\n'.format(
                                atom.index, 0, 0, 0))

            for pair in mol_type.pair_forces:
                if not isinstance(pair, (LjDefaultPairType, LjqDefaultPairType)):
                    atom1_type = int(atom_list[pair.atom1].split()[2])
                    atom2_type = int(atom_list[pair.atom2].split()[2])
                    if atom2_type < atom1_type:  # LAMMPS requires i < j
                        atom1_type, atom2_type = atom2_type, atom1_type
                    pair_coeffs.append('pair_coeff {0:d} {1:d} {2:11.7f} {3:11.7f}\n'.format(
                                atom1_type,
                                atom2_type,
                                pair.epsilon.value_in_unit(self.ENERGY),
                                pair.sigma.value_in_unit(self.DIST)))

        bond_list = self.force_dict['Bond']
        angle_list = self.force_dict['Angle']
        dihedral_list = self.force_dict['Dihedral']
        improper_list = self.force_dict['Improper']

        bond_coeffs = self.force_dict['Bond Coeffs']
        angle_coeffs = self.force_dict['Angle Coeffs']
        dihedral_coeffs = self.force_dict['Dihedral Coeffs']
        improper_coeffs = self.force_dict['Improper Coeffs']

        bond_styles = self.style_dict['Bond']
        angle_styles = self.style_dict['Angle']
        dihedral_styles = self.style_dict['Dihedral']
        improper_styles = self.style_dict['Improper']

        # Write the actual data file.
        with open(self.data_file, 'w') as f:
            # Front matter.
            f.write(self.system.name + '\n')
            f.write('\n')

            n_atoms = len(atom_list) - 1
            n_bonds = len(bond_list) - 1
            n_angles = len(angle_list) - 1
            n_dihedrals = len(dihedral_list) - 1
            n_impropers = len(improper_list) - 1

            n_atom_types = len(atom_type_dict)
            n_bond_types = len(bond_coeffs) - 1
            n_angle_types = len(angle_coeffs) - 1
            n_dihedral_types = len(dihedral_coeffs) - 1
            n_improper_types = len(improper_coeffs) - 1

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

            # Shifting of box dimensions.
            f.write('{0:11.7f} {1:11.7f} xlo xhi\n'.format(
                    x_min, x_min + self.system.box_vector[0][0].value_in_unit(
                            self.DIST)))
            f.write('{0:11.7f} {1:11.7f} ylo yhi\n'.format(
                    y_min, y_min + self.system.box_vector[1][1].value_in_unit(
                            self.DIST)))
            f.write('{0:11.7f} {1:11.7f} zlo zhi\n'.format(
                    z_min, z_min + self.system.box_vector[2][2].value_in_unit(
                            self.DIST)))
            f.write('{0:11.7f} {1:11.7f} {2:11.7f} xy xz yz\n'.format(
                    self.system.box_vector[1][0].value_in_unit(self.DIST),
                    self.system.box_vector[2][0].value_in_unit(self.DIST),
                    self.system.box_vector[2][1].value_in_unit(self.DIST)))

            for mass in mass_list:
                f.write(mass)

            # Forcefield coefficients.
            coeff_types = [bond_coeffs, angle_coeffs,
                           dihedral_coeffs, improper_coeffs]
            for coefficients in coeff_types:
                if len(coefficients) > 1:
                    for coeff in coefficients:
                        f.write(coeff)

            # Atoms and velocities.
            for atom in atom_list:
                f.write(atom)
            for vel in vel_list:
                f.write(vel)

            # Forces.
            force_lists = [bond_list, angle_list, dihedral_list, improper_list]
            for force_list in force_lists:
                if len(force_list) > 1:
                    for force in force_list:
                        f.write(force)

        # Write the corresponding input file.
        with open(self.in_file, 'w') as f:
            f.write('units {0}\n'.format(unit_set))
            f.write('atom_style full\n')  # TODO
            f.write('\n')

            f.write('dimension 3\n')  # TODO
            f.write('boundary p p p\n')  # TODO
            f.write('\n')

            # bonded
            if len(bond_coeffs) > 1:
                f.write('bond_style hybrid {0}\n'.format(
                    " ".join(bond_styles)))
            if len(angle_coeffs) > 1:
                f.write('angle_style hybrid {0}\n'.format(
                    " ".join(angle_styles)))
            if len(dihedral_coeffs) > 1:
                f.write('dihedral_style hybrid {0}\n'.format(
                    " ".join(dihedral_styles)))
            if len(improper_coeffs) > 1:
                f.write('improper_style hybrid {0}\n'.format(
                    " ".join(improper_styles)))

            f.write('special_bonds lj {0} {1} {2} coul {3} {4} {5}\n'.format(
                0.0, 0.0, self.system.lj_correction,
                0.0, 0.0, self.system.coulomb_correction))
            f.write('\n')

            # Specify the path to the corresponding data file that we just wrote.
            f.write('read_data {0}\n'.format(os.path.basename(self.data_file)))
            f.write('\n')

            # non-bonded: either defaults, or specified by user
            f.write(nonbonded_style)

            for line in pair_coeffs:
                f.write(line)
            f.write('\n')

            if self.system.combination_rule == 'Lorentz-Berthelot':
                f.write('pair_modify mix arithmetic\n')
            elif self.system.combination_rule == 'Multiply-Sigeps':
                f.write('pair_modify mix geometric\n')
            else:
                logger.warning("Unsupported pair combination rule on writing input file!")
            f.write('\n')

            if len(mol_type.rigidwaters) > 0:
                f.write('fix settle all shake 0.000001 100 0 t')

                for rigidwater in mol_type.rigidwaters:
                    molecules = list(mol_type.molecules)
                    a1 = atom_type_dict[molecules[0].atoms[rigidwater.atom1-1].atomtype[0]]
                    a2 = atom_type_dict[molecules[0].atoms[rigidwater.atom1].atomtype[0]]
                    # get the atom types of the first two molecules

                    # settles should be numbered from 0, not 1?

                    # first, write out all the atom types involved: should be the first two in the molecule.
                    f.write(' {0:d} {1:d} a'.format(a1,a2))
                    angle_i = 0
                    # add up the angles until the settle.  I think this is problematic because
                    # it doesn't take into account any transformation to the number of
                    # angles when lammps writes out.
                    for mol_name_j, mol_type_j in self.system.molecule_types.items():
                        if mol_name_j != mol_name:
                            angle_i += len(mol_type_j.angle_forces)*len(mol_type_j.molecules)
                        elif mol_name_j == mol_name:
                            break

                # only one angle per settle
                angle_range = np.arange(angle_i+1,angle_i+len(mol_type.molecules)+1)
                for a in angle_range:
                    f.write(' {0:d}'.format(a))
                    if (a-angle_i)%10 == 0:
                        f.write(' &\n')
                f.write('\n')

            # Specify the output energies that we are interested in.
            energy_terms = " ".join(['ebond', 'eangle', 'edihed', 'eimp',
                                     'epair', 'evdwl', 'ecoul', 'elong',
                                     'etail', 'pe'])

            f.write('thermo_style custom {0}\n'.format(energy_terms))
            f.write('\n')

            f.write('run 0\n')
