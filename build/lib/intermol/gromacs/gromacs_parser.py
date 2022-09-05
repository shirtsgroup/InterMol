from collections import OrderedDict, defaultdict
import logging
import os
import math

import parmed.unit as units
from intermol.atom import Atom

from intermol.forces import *
import intermol.forces.forcefunctions as ff
from intermol.exceptions import (UnimplementedFunctional, UnsupportedFunctional,
                                 UnimplementedSetting, UnsupportedSetting,
                                 GromacsError, InterMolError)
from intermol.molecule import Molecule
from intermol.moleculetype import MoleculeType
from intermol.system import System
from intermol.gromacs.grofile_parser import GromacsGroParser


logger = logging.getLogger('InterMolLog')

ENGINE = 'gromacs'


def load(top_file, gro_file, include_dir=None, defines=None):
    """Load a set of GROMACS input files into a `System`.

    Args:
        top_filename:
        gro_file:
        include_dir:
        defines:
    Returns:
        system:
    """
    parser = GromacsParser(top_file, gro_file,
                           include_dir=include_dir, defines=defines)
    return parser.read()


def save(top_file, gro_file, system):
    """Load a set of GROMACS input files into a `System`.

    Args:
        top_filename:
        gro_file:
        include_dir:
        defines:
    Returns:
        system:
    """
    parser = GromacsParser(top_file, gro_file, system)
    return parser.write()


def default_gromacs_include_dir():
    """Find the location where gromacs #include files are referenced from, by
    searching for (1) gromacs environment variables, (2) just using the default
    gromacs install location, /usr/local/gromacs/share/gromacs/top. """
    if 'GMXLIB' in os.environ:
        return os.environ['GMXLIB']
    if 'GMXDATA' in os.environ:
        return os.path.join(os.environ['GMXDATA'], 'top')
    if 'GMXBIN' in os.environ:
        return os.path.abspath(os.path.join(
            os.environ['GMXBIN'], '..', 'share', 'gromacs', 'top'))
    return '/usr/local/share/gromacs/top'


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

    gromacs_pairs = {
        # First three correspond to pairtype 1, last two pairtype 2.
        # Letter is arbitrary.
        '1A': LjCPair,
        '1B': LjSigepsPair,
        '1C': LjDefaultPair,
        '2A': LjqCPair,
        '2B': LjqSigepsPair,
        '2C': LjqDefaultPair
        }

    lookup_gromacs_pairs = dict((v, k) for k, v in gromacs_pairs.items())

    gromacs_pair_types = dict(
        (k, eval(v.__name__ + 'Type')) for k, v in gromacs_pairs.items())

    gromacs_virtuals = {
        '2-1': TwoVirtual,
        '3-1': ThreeLinearVirtual,
        '3-2': ThreeFdVirtual,
        '3-3': ThreeFadVirtual,
        '3-4': ThreeOutVirtual,
        '4-2': FourFdnVirtual
    }

    lookup_gromacs_virtuals = dict((v, k) for k, v in gromacs_virtuals.items())

    gromacs_virtual_types = dict(
        (k, eval(v.__name__ + 'Type')) for k, v in gromacs_virtuals.items())

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
            try:
                b_type = self.lookup_gromacs_bonds[bond.__class__]
            except KeyError:
                raise UnsupportedFunctional(bond, ENGINE)
            return b_type, params

    gromacs_angles = {
        '1': HarmonicAngle,
        '2': CosineSquaredAngle,
        '3': CrossBondBondAngle,
        '4': CrossBondAngleAngle,
        '5': UreyBradleyAngle,
        '6': QuarticAngle,
        '10': RestrictedBendingAngle
        }
    lookup_gromacs_angles = dict((v, k) for k, v in gromacs_angles.items())

    gromacs_angle_types = dict(
        (k, eval(v.__name__ + 'Type')) for k, v in gromacs_angles.items())

    def canonical_angle(self, params, angle, direction='into'):
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
            try:
                a_type = self.lookup_gromacs_angles[angle.__class__]
            except KeyError:
                raise UnsupportedFunctional(angle, ENGINE)
            return a_type, params

    gromacs_dihedrals = {
        # TrigDihedrals are actually used for 1, 4, and 9.  Can't use lists as keys!
        '1': ProperPeriodicDihedral,
        '2': ImproperHarmonicDihedral,
        '3': RbDihedral,
        '4': ProperPeriodicDihedral,
        '5': FourierDihedral,
        '9': ProperPeriodicDihedral,
        '10': RestrictedBendingDihedral,
        '11': BendingTorsionDihedral,
        'Trig': TrigDihedral
        }

    # have to invert manually because of canonical conversion collapse of types.
    lookup_gromacs_dihedrals = {
        TrigDihedral: 'Trig',
        ImproperHarmonicDihedral: '2',
        RbDihedral: '3',
        FourierDihedral: '5',
        RestrictedBendingDihedral: '10',
        BendingTorsionDihedral: '11'
        }

    gromacs_dihedral_types = dict(
        (k, eval(v.__name__ + 'Type')) for k, v in gromacs_dihedrals.items())

    def canonical_dihedral(self, params, dihedral, direction='into'):
        """

        We can fit everything into two types of dihedrals - dihedral_trig, and
        improper harmonic. Dihedral trig is of the form

            fc0 + sum_i=1^6 fci (cos(nx-phi)

        Proper dihedrals can be stored easily in this form, since they have
        only 1 n. Improper dihedrals can as well (flag as improper). RB can be
        stored as well, assuming phi = 0 or 180. Fourier can also be stored. A
        full dihedral trig can be decomposed into multiple proper dihedrals.

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
            if dihedral == ProperPeriodicDihedralType:
                convertfunc = convert_dihedral_from_proper_to_trig
                converted_dihedral = TrigDihedralType
            elif dihedral == ProperPeriodicDihedral:
                convertfunc = convert_dihedral_from_proper_to_trig
                converted_dihedral = TrigDihedral
            elif dihedral == RbDihedralType:
                convertfunc = convert_dihedral_from_RB_to_trig
                # Sign convention from psi to phi.
                params['C1'] *= -1
                params['C3'] *= -1
                params['C5'] *= -1
                converted_dihedral = TrigDihedralType
            elif dihedral == RbDihedral:
                convertfunc = convert_dihedral_from_RB_to_trig
                # Sign convention from psi to phi.
                params['C1'] *= -1
                params['C3'] *= -1
                params['C5'] *= -1
                converted_dihedral = TrigDihedral
            elif dihedral == FourierDihedralType:
                convertfunc = convert_dihedral_from_fourier_to_trig
                converted_dihedral = TrigDihedralType
            elif dihedral == FourierDihedral:
                convertfunc = convert_dihedral_from_fourier_to_trig
                converted_dihedral = TrigDihedral
            elif dihedral in (ImproperHarmonicDihedralType, ImproperHarmonicDihedral,
                              TrigDihedralType, TrigDihedral,
                              BendingTorsionDihedralType, BendingTorsionDihedral,
                              RestrictedBendingDihedralType, RestrictedBendingDihedral
                              ):
                convertfunc = convert_nothing
                converted_dihedral = dihedral
            else:
                raise GromacsError('Unable to convert dihedral: {0}'.format(dihedral))
            params = convertfunc(params)
            return converted_dihedral, params
        else:
            if isinstance(dihedral, TrigDihedral):
                if dihedral.improper:
                    d_type = '4'
                    paramlist = convert_dihedral_from_trig_to_proper(params)
                else:
                    if (params['phi'].value_in_unit(units.degrees) in [0, 180] and
                                params['fc6']._value == 0):
                        d_type = '3'
                        params = convert_dihedral_from_trig_to_RB(params)
                        # Sign convention from phi to psi.
                        params['C1'] *= -1
                        params['C3'] *= -1
                        params['C5'] *= -1
                        paramlist = [params]
                    else:
                        # Print as proper dihedral. If one nonzero term, as a
                        # type 1, if multiple, type 9.
                        paramlist = convert_dihedral_from_trig_to_proper(params)
                        if len(paramlist) == 1:
                            d_type = '1'
                        else:
                            d_type = '9'
            else:
                try:
                    d_type = self.lookup_gromacs_dihedrals[dihedral.__class__]
                except KeyError:
                    raise UnsupportedFunctional(dihedral, ENGINE)
                paramlist = [params]
            return d_type, paramlist

    def choose_parameter_kwds_from_forces(self, entries, n_atoms, force_type,
                                          gromacs_force):
        """Extract a force's parameters into a keyword dictionary.

        Args:
            entries (str): The `split()` line being parsed.
            n_atoms (int): The number of atoms in the force.
            force_type: The type of the force.
            gromacs_force: The
        Returns:
            kwds (dict): The force's parameters, e.g.
                {'length': Quantity(value=0.13, unit=nanometers),
                 'k': ...
                }
        """
        n_entries = len(entries)
        gromacs_force_type = gromacs_force.__base__  # what's the base class
        typename = gromacs_force_type.__name__
        u = self.unitvars[typename]
        params = self.paramlist[typename]
        kwds = dict()
        if n_entries > n_atoms + 2:
            for i, p in enumerate(params):
                kwds[p] = float(entries[n_atoms + 1 + i]) * u[i]
        elif n_entries in [n_atoms + 1 or n_atoms + 2]:
            # Check to see if the force is defined exists
            if isinstance(force_type, gromacs_force_type):
                force_type_params = self.get_parameter_list_from_force(force_type)
                # Note: for now, not passing the bonding variables.
                for i, p in enumerate(params):
                    kwds[p] = force_type_params[i]
            else:
                logger.warning("No forcetype defined for: {0}".format(entries))
        return kwds

    paramlist = ff.build_paramlist('gromacs')
    unitvars = ff.build_unitvars('gromacs', paramlist)

    def create_kwds_from_entries(self, entries, force_class, offset=0):
        return ff.create_kwds_from_entries(self.unitvars, self.paramlist,
                                           entries, force_class, offset=offset)

    def get_parameter_list_from_force(self, force):
        return ff.get_parameter_list_from_force(force, self.paramlist)

    def get_parameter_kwds_from_force(self, force):
        return ff.get_parameter_kwds_from_force(
                force, self.get_parameter_list_from_force, self.paramlist)

    class TopMoleculeType(object):
        """Inner class to store information about a molecule type."""
        def __init__(self):
            self.nrexcl = -1
            self.atoms = []
            self.bonds = []
            self.angles = []
            self.dihedrals = []
            self.rigidwaters = []
            self.exclusions = []
            self.pairs = []
            self.cmaps = []
            self.virtuals = defaultdict(list)

    def __init__(self, top_file, gro_file, system=None, include_dir=None, defines=None):
        """Initializes the parser with all required metadata.

        Args:
            defines: Sets of default defines to use while parsing.
        """
        self.top_filename = top_file
        self.gro_file = gro_file
        if not system:
            system = System()
        self.system = system

        if include_dir is None:
            include_dir = default_gromacs_include_dir()
        self.include_dirs = (os.path.dirname(top_file), include_dir)
        # Most of the gromacs water itp files for different forcefields,
        # unless the preprocessor #define FLEXIBLE is given, don't define
        # bonds between the water hydrogen and oxygens, but only give the
        # constraint distances and exclusions.
        self.defines = dict()
        if defines is not None:
            self.defines.update(defines)

    def read(self):
        """Load the files into InterMol's abstract representation.

        Returns:
            system
        """
        self.current_directive = None
        self.if_stack = list()
        self.else_stack = list()
        self.molecule_types = OrderedDict()
        self.molecules = list()
        self.current_molecule_type = None
        self.current_molecule = None
        self.bondtypes = dict()
        self.angletypes = dict()
        self.dihedraltypes = dict()
        self.implicittypes = dict()
        self.pairtypes = dict()
        self.cmaptypes = dict()
        self.nonbondedtypes = dict()

        # Parse the top_filename into a set of plain text, intermediate
        # TopMoleculeType objects.
        self.process_file(self.top_filename)

        # Open the corresponding gro file and push all the information to the
        # InterMol system.
        self.gro = GromacsGroParser(self.gro_file)
        self.gro.read()
        self.system.box_vector = self.gro.box_vector
        self.system.n_atoms = self.gro.positions.shape[0]

        self.n_atoms_added = 0
        for mol_name, mol_count in self.molecules:
            if mol_name not in self.molecule_types:
                raise GromacsError("Unknown molecule type: {0}".format(mol_name))
            # Grab the relevent plain text molecule type.
            top_moltype = self.molecule_types[mol_name]
            self.create_moleculetype(top_moltype, mol_name, mol_count)

        return self.system

    # =========== System writing =========== #
    def write(self):
        """Write this topology in GROMACS file format.

        Args:
            filename: the name of the file to write out to
        """
        gro = GromacsGroParser(self.gro_file)
        gro.write(self.system)

        with open(self.top_filename, 'w') as top:
            self.write_defaults(top)
            self.write_atomtypes(top)
            if self.system.nonbonded_types:
                self.write_nonbonded_types(top)

            self.write_moleculetypes(top)

            self.write_system(top)
            self.write_molecules(top)

    def write_defaults(self, top):
        top.write('[ defaults ]\n')
        top.write('; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n')
        top.write('{0:6d} {1:6s} {2:6s} {3:8.6f} {4:8.6f}\n\n'.format(
                   self.system.nonbonded_function,
                   self.lookup_gromacs_combination_rules[self.system.combination_rule],
                   self.system.genpairs,
                   self.system.lj_correction,
                   self.system.coulomb_correction))

    def write_atomtypes(self, top):
        top.write('[ atomtypes ]\n')
        top.write(';type, bondingtype, atomic_number, mass, charge, ptype, sigma, epsilon\n')
        for atomtype in sorted(self.system.atomtypes.values(), key=lambda x: x.atomtype):
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

    def write_nonbonded_types(self, top):
        top.write('[ nonbond_params ]\n')
        top.write(';i    j    func    sigma     epsilon\n')
        for nbtype in sorted(self.system.nonbonded_types.values(), key=lambda x: (x.atom1, x.atom2)):
            # TODO: support for buckingham NB types
            top.write('{0:6s} {1:6s} {2:3d}'.format(
                    nbtype.atom1, nbtype.atom2, nbtype.type))
            if self.system.combination_rule == 'Multiply-C6C12':
                top.write('{0:18.8e} {1:18.8e}\n'.format(
                    nbtype.C6.value_in_unit(units.kilojoules_per_mole * units.nanometers**(6)),
                    nbtype.C12.value_in_unit(units.kilojoules_per_mole * units.nanometers**(12))))
            elif self.system.combination_rule in ['Lorentz-Berthelot', 'Multiply-Sigeps']:
                top.write('{0:18.8e} {1:18.8e}\n'.format(
                    nbtype.sigma.value_in_unit(units.nanometers),
                    nbtype.epsilon.value_in_unit(units.kilojoules_per_mole)))
        top.write('\n')

    def write_moleculetypes(self, top):
        for mol_name, mol_type in self.system.molecule_types.items():
            self.current_molecule_type = mol_type
            top.write('[ moleculetype ]\n')
            # Gromacs can't handle spaces in the molecule name.
            printname = mol_name
            printname = printname.replace(' ', '_')
            printname = printname.replace('"', '')
            top.write('{0:s} {1:10d}\n\n'.format(printname, mol_type.nrexcl))

            self.write_atoms(top)

            if self.current_molecule_type.pair_forces:
                self.write_pairs(top)
            if self.current_molecule_type.bond_forces and not self.current_molecule_type.rigidwaters:
                self.write_bonds(top)
            if self.current_molecule_type.angle_forces and not self.current_molecule_type.rigidwaters:
                self.write_angles(top)
            if self.current_molecule_type.dihedral_forces:
                self.write_dihedrals(top)
            if self.current_molecule_type.virtual_forces:
                self.write_virtual_sites(top)
            if self.current_molecule_type.rigidwaters:
                self.write_rigidwaters(top)
            if self.current_molecule_type.exclusions:
                self.write_exclusions(top)

    def write_system(self, top):
        top.write('[ system ]\n')
        top.write('{0}\n\n'.format(self.system.name))

    def write_molecules(self, top):
        top.write('[ molecules ]\n')
        top.write('; Compound        nmols\n')
        for mol_name, mol_type in self.system.molecule_types.items():
            n_molecules = len(mol_type.molecules)
            # The following lines are more 'chemical'.
            printname = mol_name
            printname = printname.replace(' ', '_')
            printname = printname.replace('"', '')
            top.write('{0:<15s} {1:8d}\n'.format(printname, n_molecules))

    def write_atoms(self, top):
        top.write('[ atoms ]\n')
        top.write(';num, type, resnum, resname, atomname, cgnr, q, m\n')

        # Start iterating the set to get the first entry (somewhat kludgy...)
        for i, atom in enumerate(next(iter(self.current_molecule_type.molecules)).atoms):
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

    def write_pairs(self, top):
        top.write('[ pairs ]\n')
        top.write(';  ai    aj   funct\n')
        pairlist = sorted(self.current_molecule_type.pair_forces,
                          key=lambda x: (x.atom1, x.atom2))
        for pair in pairlist:
            p_type = self.lookup_gromacs_pairs[pair.__class__]
            if p_type:
                # Gromacs type is the first character
                top.write('{0:6d} {1:7d} {2:4d}'.format(
                    pair.atom1, pair.atom2, int(p_type[0])))

                pair_params = self.get_parameter_list_from_force(pair)
                # Don't want to write over actual array.
                param_units = list(self.unitvars[pair.__class__.__name__])
                if p_type[0] == '2' and pair.scaleQQ:
                    # We have a scaleQQ as well, which has no units.
                    pair_params.insert(0, pair.scaleQQ)
                    param_units.insert(0, units.dimensionless)
                for i, param in enumerate(pair_params):
                        top.write("{0:18.8e}".format(
                                param.value_in_unit(param_units[i])))
                top.write('\n')
            else:
                logger.warning("Found unsupported pair type {0}".format(
                        pair.__class__.__name__))

        top.write('\n')

    def write_virtual_sites(self, top):
        virtuals = defaultdict(list)
        for force in self.current_molecule_type.virtual_forces:
            if hasattr(force, 'atom5'):
                virtuals[4].append(force)
            elif hasattr(force, 'atom4'):
                virtuals[3].append(force)
            else:
                virtuals[2].append(force)

        virtuals[2] = sorted(virtuals[2], key=lambda x: (x.atom1, x.atom2, x.atom3))
        virtuals[3] = sorted(virtuals[3], key=lambda x: (x.atom1, x.atom2, x.atom3, x.atom4))
        virtuals[4] = sorted(virtuals[4], key=lambda x: (x.atom1, x.atom2, x.atom3, x.atom4, x.atom5))

        for n_body_type, vsites in virtuals.items():
            top.write('[ virtual_sites{0} ]\n'.format(n_body_type))
            top.write(';from   atoms({0})   func      params\n'.format(n_body_type))
            for vsite in vsites:
                for n in range(1, n_body_type + 2):
                    atom = getattr(vsite, 'atom{}'.format(n))
                    top.write('{0:7d} '.format(atom))

                top.write('{:4s}'.format(self.lookup_gromacs_virtuals[vsite.__class__][-1]))

                vsite_params = self.get_parameter_list_from_force(vsite)
                param_units = self.unitvars[vsite.__class__.__name__]
                for param, unit in zip(vsite_params, param_units):
                    top.write('{0:18.8e}'.format(param.value_in_unit(unit)))
                top.write('\n')
            top.write('\n')

    def write_bonds(self, top):
        top.write('[ bonds ]\n')
        top.write(';   ai     aj funct  r               k\n')
        bondlist = sorted(self.current_molecule_type.bond_forces,
                          key=lambda x: (x.atom1, x.atom2))
        for bond in bondlist:
            bond_params = self.get_parameter_list_from_force(bond)
            b_type, bond_params = self.canonical_bond(bond_params, bond, direction='from')
            top.write('{0:7d} {1:7d} {2:4s}'.format(
                bond.atom1, bond.atom2, b_type))

            param_units = self.unitvars[bond.__class__.__name__]
            for param, param_unit in zip(bond_params, param_units):
                top.write('{0:18.8e}'.format(param.value_in_unit(param_unit)))
            top.write('\n')
        top.write('\n')

    def write_angles(self, top):
        top.write('[ angles ]\n')
        top.write(';   ai     aj     ak     funct  theta    cth\n')
        anglelist = sorted(self.current_molecule_type.angle_forces,
                           key=lambda x: (x.atom1, x.atom2, x.atom3))
        for angle in anglelist:
            angle_params = self.get_parameter_list_from_force(angle)
            a_type, angle_params = self.canonical_angle(angle_params, angle, direction='from')
            top.write('{0:7d} {1:7d} {2:7d} {3:4s}'.format(
                angle.atom1, angle.atom2, angle.atom3, a_type))

            param_units = self.unitvars[angle.__class__.__name__]
            for param, param_unit in zip(angle_params, param_units):
                top.write('{0:18.8e}'.format(param.value_in_unit(param_unit)))
            top.write('\n')
        top.write('\n')

    def write_dihedrals(self, top):
        top.write('[ dihedrals ]\n')
        top.write(';    i      j      k      l   func\n')
        dihedrallist = sorted(self.current_molecule_type.dihedral_forces,
                              key=lambda x: (x.atom1, x.atom2, x.atom3, x.atom4))
        for dihedral in dihedrallist:
            atoms = dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4
            top.write("{0:7d} {1:7d} {2:7d} {3:7d}".format(
                   atoms[0], atoms[1], atoms[2], atoms[3]))

            kwds = self.get_parameter_kwds_from_force(dihedral)
            d_type, paramlist = self.canonical_dihedral(kwds, dihedral, direction='from')
            converted_dihedral = self.gromacs_dihedrals[d_type](*atoms, **paramlist[0])

            top.write("{0:6d}".format(int(d_type)))

            paramlist = self.get_parameter_list_from_force(converted_dihedral)
            param_units = self.unitvars[converted_dihedral.__class__.__name__]
            for param, param_unit in zip(paramlist, param_units):
                top.write('{0:18.8e}'.format(param.value_in_unit(param_unit)))
            top.write('\n')
        top.write('\n')

    def write_rigidwaters(self, top):
        for rigidwater in self.current_molecule_type.rigidwaters:
            top.write('[ settles ]\n')
            top.write('; i  funct   dOH  dHH\n')
            s_type = 1
            # gromacs only uses the first atom of the rigid waters as the O, expects the others follow in sequence
            top.write('{0:6d} {1:6d} {2:18.8f} {3:18.8f}\n'.format(
                    rigidwater.atom1,
                    s_type,
                    rigidwater.dOH.value_in_unit(units.nanometers),
                    rigidwater.dHH.value_in_unit(units.nanometers)))
            top.write('\n')

    def write_exclusions(self, top):
        top.write('[ exclusions ]\n')
        exclusionlist = sorted(self.current_molecule_type.exclusions,
                               key=lambda x: (x[0], x[1]))
        for exclusion in exclusionlist:
            top.write('{0:7d} {1:7d}\n'.format(exclusion[0], exclusion[1]))
        top.write('\n')

    # =========== System creation =========== #
    def create_moleculetype(self, top_moltype, mol_name, mol_count):
        # Check if the moleculetype already exists
        if self.system.molecule_types.get(mol_name):
            self.current_molecule_type = self.system.molecule_types[mol_name]
        else:
            # Create an intermol moleculetype.
            moltype = MoleculeType(mol_name)
            moltype.nrexcl = top_moltype.nrexcl
            self.system.add_molecule_type(moltype)
            self.current_molecule_type = moltype

        # Create all the intermol molecules of the current type.
        for n_mol in range(mol_count):
            self.create_molecule(top_moltype, mol_name)
        for pair in top_moltype.pairs:
            self.create_pair(pair)
        for bond in top_moltype.bonds:
            self.create_bond(bond)
        for angle in top_moltype.angles:
            self.create_angle(angle)
        for dihedral in top_moltype.dihedrals:
            self.create_dihedral(dihedral)
        for rigidwater in top_moltype.rigidwaters:
            self.create_rigidwater(rigidwater)
        for exclusion in top_moltype.exclusions:
            self.create_exclusion(exclusion)
        for vsite_type, vsites in top_moltype.virtuals.items():
            for vsite in vsites:
                self.create_virtual_site(vsite, vsite_type)

    def create_molecule(self, top_moltype, mol_name):
        molecule = Molecule(mol_name)
        self.system.add_molecule(molecule)
        self.current_molecule = molecule
        for atom in top_moltype.atoms:
            self.create_atom(atom)

    def create_atom(self, temp_atom):
        index = self.n_atoms_added + 1
        atomtype = temp_atom[1]
        res_id = self.gro.residue_ids[self.n_atoms_added]
        res_name = self.gro.residue_names[self.n_atoms_added]
        atom_name = temp_atom[4]
        cgnr = int(temp_atom[5])
        charge = float(temp_atom[6]) * units.elementary_charge
        if len(temp_atom) in [8, 11]:
            mass = float(temp_atom[7]) * units.amu
        else:
            mass = -1 * units.amu

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
        atom.velocity = self.gro.velocities[self.n_atoms_added]

        for state, atomtype in atom.atomtype.items():
            intermol_atomtype = self.system.atomtypes.get(atomtype)
            if not intermol_atomtype:
                logger.warning('A corresponding AtomType for {0} was not'
                            ' found.'.format(atom))
                continue
            atom.atomic_number = intermol_atomtype.atomic_number
            if not atom.bondingtype:
                if intermol_atomtype.bondtype:
                    atom.bondingtype = intermol_atomtype.bondtype
                else:
                    atom.bondingtype = atomtype
            if atom.mass.get(state)._value < 0:
                if intermol_atomtype.mass._value >= 0:
                    atom.mass = (state, intermol_atomtype.mass)
                else:
                    logger.warning("Suspicious mass parameter found for atom "
                                "{0}. Visually inspect before using.".format(atom))
            atom.sigma = (state, intermol_atomtype.sigma)
            atom.epsilon = (state, intermol_atomtype.epsilon)

        self.current_molecule.add_atom(atom)
        self.n_atoms_added += 1

    def create_bond(self, bond):
        n_atoms = 2
        numeric_bondtype = bond[n_atoms]
        atoms = [int(n) for n in bond[:n_atoms]]
        btypes = tuple([self.lookup_atom_bondingtype(int(x))
                        for x in bond[:n_atoms]])

        # Get forcefield parameters.
        if len(bond) == n_atoms + 1:
            bond_type = self.find_forcetype(btypes, self.bondtypes)
        else:
            bond[0] = btypes[0]
            bond[1] = btypes[1]
            bond = " ".join(bond)
            bond_type = self.process_forcetype(btypes, 'bond', bond, n_atoms,
                self.gromacs_bond_types, self.canonical_bond)
            bond = bond.split()

        # Create the actual force.
        if numeric_bondtype in self.gromacs_bonds:
            gromacs_bond = self.gromacs_bonds[numeric_bondtype]
            # Connection bonds don't have bondtypes.
            if gromacs_bond == ConnectionBond:
                kwds = dict()
            else:
                kwds = self.choose_parameter_kwds_from_forces(
                    bond, n_atoms, bond_type, gromacs_bond)
            # Give it canonical form parameters.
            canonical_bond, kwds = self.canonical_bond(kwds, gromacs_bond,
                                                        direction='into')
            new_bond = canonical_bond(*atoms, **kwds)
        else:
            logger.warning("Unsupported Gromacs bondtype: {0}".format(numeric_bondtype))

        if not new_bond:
            logger.warning("Undefined bond formatting.")
        else:
            self.current_molecule_type.bond_forces.add(new_bond)

    def create_pair(self, pair):
        """Create a pair force object based on a [ pairs ] entry"""

        n_entries = len(pair)
        numeric_pairtype = pair[2]
        atoms = [int(pair[0]), int(pair[1])]
        atomtypes = tuple([self.lookup_atom_atomtype(int(pair[0])),
                          self.lookup_atom_atomtype(int(pair[1]))])
        if n_entries == 3:
            pairtype = self.find_forcetype(atomtypes, self.pairtypes)
        else:
            atomtypes = [None, None]

        pairvars = [atoms[0], atoms[1], atomtypes[0], atomtypes[1]]
        optpairvars = dict()
        if numeric_pairtype == '1':
            if self.system.combination_rule == "Multiply-C6C12":
                thispair = LjCPair
            elif self.system.combination_rule in ['Multiply-Sigeps', 'Lorentz-Berthelot']:
                thispair = LjSigepsPair
            thispairtype = thispair.__base__  # what's the base class

            u = self.unitvars[thispairtype.__name__]
            if n_entries > 3:
                pairvars.extend([float(pair[3]) * u[0], float(pair[4]) * u[1]])
            elif n_entries == 3:
                if not pairtype:
                    # assume the values will be created by system defaults
                    thispair = LjDefaultPair
                else:
                    assert isinstance(pairtype, thispairtype)
                    pairvars.extend(self.get_parameter_list_from_force(pairtype))
            new_pair = thispair(*pairvars)
        elif numeric_pairtype == '2':
            if self.system.combination_rule == "Multiply-C6C12":
                thispair = LjqCPair
            elif self.system.combination_rule in ['Multiply-Sigeps', 'Lorentz-Berthelot']:
                thispair = LjqSigepsPair
            thispairtype = thispair.__base__  # what's the parent?

            u = self.unitvars[thispairtype.__name__]
            if n_entries > 3:
                pairvars.extend([float(pair[4]) * u[0], float(pair[5]) * u[1],
                                 float(pair[6]) * u[2], float(pair[7]) * u[3]])
                # Generate a default filled dictionary, then fill in the pair.
                optpairvars = ff.optforceparams('pair')
                optpairvars['scaleQQ'] = float(pair[3]) * units.dimensionless
            elif n_entries == 3:
                if not pairtype:
                    # Assume the values will be created by system defaults.
                    thispair = LjqDefaultPair
                else:
                    assert isinstance(pairtype, thispairtype)
                    # Bring the data from this pairtype.
                    optpairvars['scaleQQ'] = pairtype.scaleQQ
                    pairvars.extend(self.get_parameter_list_from_force(pairtype))
            new_pair = thispair(*pairvars, **optpairvars)
        else:
            logger.warning("Unsupported Gromacs pairtype: {0}".format(
                numeric_pairtype))

        if not new_pair:
            logger.warning("Undefined pair formatting.")
        else:
            self.current_molecule_type.pair_forces.add(new_pair)

    def create_rigidwater(self, rigidwater):
        # gromacs just assumes the first atom is the oxygen.
        atom1 = int(rigidwater[0])
        atom2 = atom1 + 1
        atom3 = atom1 + 2
        new_rigidwater = RigidWater(atom1, atom2, atom3,
                            float(rigidwater[2]) * units.nanometers,
                            float(rigidwater[3]) * units.nanometers)
        self.current_molecule_type.rigidwaters.add(new_rigidwater)

        waterbondrefk = 900*units.kilojoules_per_mole * units.nanometers**(-2)
        wateranglerefk = 400*units.kilojoules_per_mole * units.degrees**(-2)
        angle = 2.0 * math.asin(0.5 * float(rigidwater[3]) / float(rigidwater[2])) * units.radians
        dOH = float(rigidwater[2]) * units.nanometers

        new_bond = HarmonicBond(atom1, atom2, None, None, dOH, waterbondrefk, c=True)
        self.current_molecule_type.bond_forces.add(new_bond)

        new_bond = HarmonicBond(atom1, atom3, None, None, dOH, waterbondrefk, c=True)
        self.current_molecule_type.bond_forces.add(new_bond)

        new_angle = HarmonicAngle(atom2, atom1, atom3, None, None, None, angle, wateranglerefk, c=True)
        self.current_molecule_type.angle_forces.add(new_angle)

    def create_exclusion(self, exclusion):
        first = exclusion[0]
        for index in exclusion:
            if first < index:
                self.current_molecule_type.exclusions.add((int(first), int(index)))

    def create_virtual_site(self, vsite, n_body_type):
        n_entries = len(vsite)
        n_atoms = int(n_body_type) + 1
        numeric_vsite_type = vsite[n_atoms]
        force_lookup_key = '{}-{}'.format(n_body_type, numeric_vsite_type)
        VSite = self.gromacs_virtuals[force_lookup_key]
        VSiteType = VSite.__base__

        atoms = [int(n) for n in vsite[:n_atoms]]
        btypes = [None] * n_atoms  # TODO: Unused? Can we remove from classes?
        params = self.create_kwds_from_entries(vsite, VSiteType, offset=n_atoms + 1)

        # Can't unpack multiple args in Python <3.5 so combine into one.
        atoms.extend(btypes)

        new_vsite = VSite(*atoms, **params)
        self.current_molecule_type.virtual_forces.add(new_vsite)

    def create_angle(self, angle):
        n_atoms = 3
        atoms = [int(n) for n in angle[:n_atoms]]
        btypes = tuple([self.lookup_atom_bondingtype(int(x))
                        for x in angle[:n_atoms]])
        numeric_angletype = angle[n_atoms]

        # Get forcefield parameters.
        if len(angle) == n_atoms + 1:
            angle_type = self.find_forcetype(btypes, self.angletypes)
        else:
            angle[0] = btypes[0]
            angle[1] = btypes[1]
            angle[2] = btypes[2]
            angle = " ".join(angle)
            angle_type = self.process_forcetype(btypes, 'angle', angle, n_atoms,
                self.gromacs_angle_types, self.canonical_angle)
            angle = angle.split()

        # Create the actual force.
        if numeric_angletype in self.gromacs_angles:
            gromacs_angle = self.gromacs_angles[numeric_angletype]
            kwds = self.choose_parameter_kwds_from_forces(
                angle, n_atoms, angle_type, gromacs_angle)
            # Give it canonical form parameters.
            canonical_angle, kwds = self.canonical_angle(kwds, gromacs_angle,
                                                         direction='into')
            new_angle = canonical_angle(*atoms, **kwds)
        else:
            logger.warning("Unsupported Gromacs angletype: {0}".format(numeric_angletype))

        if not new_angle:
            logger.warning("Undefined angle formatting.")
        else:
            self.current_molecule_type.angle_forces.add(new_angle)

    def create_dihedral(self, dihedral):
        """Create a dihedral object based on a [ dihedrals ] entry. """
        n_entries = len(dihedral)
        n_atoms = 4
        atoms = [int(i) for i in dihedral[0:n_atoms]]
        numeric_dihedraltype = dihedral[n_atoms]

        improper = numeric_dihedraltype in ['2', '4']

        dihedral_types = [None]
        if n_entries == n_atoms + 1:
            btypes = [self.lookup_atom_bondingtype(int(x))
                      for x in dihedral[:n_atoms]]

            # Use the returned btypes that we get a match with!
            dihedral_types = self.find_dihedraltype(btypes, improper=improper)
            # This dihedraltype has been found before and already converted.
            if numeric_dihedraltype in ['1', '3', '4', '5', '9']:
                gromacs_dihedral = TrigDihedral
            else:
                gromacs_dihedral = self.gromacs_dihedrals[numeric_dihedraltype]
        elif n_entries == n_atoms + 2:
            # This case handles special dihedral given via a #define.
            if self.defines.get(dihedral[-1]):
                params = self.defines[dihedral[-1]].split()
                dihedral = dihedral[:-1] + params
            gromacs_dihedral = self.gromacs_dihedrals[numeric_dihedraltype]
        else:
            # Some gromacs parameters don't include sufficient entries for all
            # types, so add some zeros. A bit of a kludge...
            dihedral += ['0.0'] * 3
            gromacs_dihedral = self.gromacs_dihedrals[numeric_dihedraltype]

        for d_type in dihedral_types:
            kwds = self.choose_parameter_kwds_from_forces(
                    dihedral, n_atoms, d_type, gromacs_dihedral)
            canonical_dihedral, kwds = self.canonical_dihedral(
                    kwds, gromacs_dihedral, direction="into")

            kwds['improper'] = improper
            new_dihedral = canonical_dihedral(*atoms, **kwds)
            self.current_molecule_type.dihedral_forces.add(new_dihedral)

    def find_dihedraltype(self, bondingtypes, improper):
        """Determine the type of dihedral interaction between four atoms. """
        a1, a2, a3, a4 = bondingtypes
        # All possible ways to match a dihedraltype
        atom_orders = [[a1, a2, a3, a4],    # original order
                       [a4, a3, a2, a1],    # flip it
                       [a1, a2, a3, 'X'],   # single wildcard 1
                       ['X', a3, a2, a1],   # flipped single wildcard 1
                       [a4, a3, a2, 'X'],   # flipped single wildcard 2
                       ['X', a2, a3, a4],   # single wildcard 2
                       ['X', a2, a3, 'X'],  # double wildcard
                       ['X', 'X', a3, a4],  # front end double wildcard
                       [a1, a2, 'X', 'X'],  # rear end double wildcard
                       ['X', 'X', a2, a1],  # rear end double wildcard
                       [a1, 'X', 'X', a4],  # middle double wildcard
                       ['X', a3, a2, 'X'],  # flipped double wildcard
                       [a4, a3, 'X', 'X'],  # flipped front end double wildcard
                       [a4, 'X', 'X', a1],  # flipped middle double wildcard
                       ]

        dihedral_types = set()
        for i, atoms in enumerate(atom_orders):
            a1, a2, a3, a4 = atoms
            key = tuple([a1, a2, a3, a4, improper])
            dihedral_type = self.dihedraltypes.get(key)
            if dihedral_type:
                for to_be_added in dihedral_type:
                    for already_added in dihedral_types:
                        if not self.type_parameters_are_unique(to_be_added,
                                                               already_added):
                            break
                    else:  # The loop completed without breaking.
                        dihedral_types.add(to_be_added)
                break
        if not dihedral_types:
            logger.warning("Lookup failed for dihedral: {0}".format(bondingtypes))
            return []
        else:
            return list(dihedral_types)

    @staticmethod
    def type_parameters_are_unique(a, b):
        """Check if two force types are unique.

        Currently only tests TrigDihedralType and ImproperHarmonicDihedralType
        because these are the only two forcetypes that we currently allow to
        to have multiple values for the same set of 4 atom bondingtypes.
        """
        if (isinstance(a, TrigDihedralType) and
                isinstance(b, TrigDihedralType)):
            return not (a.fc0 == b.fc0 and
                        a.fc1 == b.fc1 and
                        a.fc2 == b.fc2 and
                        a.fc3 == b.fc3 and
                        a.fc4 == b.fc4 and
                        a.fc5 == b.fc5 and
                        a.fc6 == b.fc6 and
                        a.improper == b.improper and
                        a.phi == b.phi)
        elif (isinstance(a, ImproperHarmonicDihedralType) and
                isinstance(b, ImproperHarmonicDihedralType)):
            return not (a.xi == b.xi and
                        a.k == b.k and
                        a.improper == b.improper)
        else:
            return True

    def lookup_atom_bondingtype(self, index):
        return self.current_molecule.atoms[index - 1].bondingtype

    def lookup_atom_atomtype(self, index, state=0):
        return self.current_molecule.atoms[index - 1].atomtype[state]

    def find_forcetype(self, bondingtypes, types_of_kind):
        forcetype = types_of_kind.get(bondingtypes)
        if not forcetype:
            forcetype = types_of_kind.get(bondingtypes[::-1])

        if not forcetype:
            logger.debug("Lookup failed for atom bonding types'{0}' in {1}".format(
                    bondingtypes, types_of_kind.keys()))
        return forcetype

    # =========== Pre-processing and forcetype creation =========== #
    def process_file(self, top_filename):
        append = ''
        with open(top_filename) as top_file:
            for line in top_file:
                if line.strip().endswith('\\'):
                    append = '{0} {1}'.format(append, line[:line.rfind('\\')])
                else:
                    self.process_line(top_filename, '{0} {1}'.format(append, line))
                    append = ''

    def process_line(self, top_filename, line):
        """Process one line from a file."""
        if ';' in line:
            line = line[:line.index(';')]
        stripped = line.strip()
        ignore = not all(self.if_stack)
        if stripped.startswith('*') or len(stripped) == 0:
            # A comment or empty line.
            return

        elif stripped.startswith('[') and not ignore:
            # The start of a category.
            if not stripped.endswith(']'):
                raise GromacsError('Illegal line in .top file: '+line)
            self.current_directive = stripped[1:-1].strip()
            logger.debug("Parsing {0}...".format(self.current_directive))

        elif stripped.startswith('#'):
            # A preprocessor command.
            fields = stripped.split()
            command = fields[0]
            if len(self.if_stack) != len(self.else_stack):
                raise GromacsError('#if/#else stack out of sync')

            if command == '#include' and not ignore:
                # Locate the file to include
                name = stripped[len(command):].strip(' \t"<>')
                search_dirs = self.include_dirs+(os.path.dirname(top_filename),)
                for sub_dir in search_dirs:
                    top_filename = os.path.join(sub_dir, name)
                    if os.path.isfile(top_filename):
                        # We found the file, so process it.
                        self.process_file(top_filename)
                        break
                else:
                    raise GromacsError('Could not locate #include file: {}\n\n'
                                       'Did you add the GROMACS share directory'
                                       ' to "GMXDATA"?'.format(name))

            elif command == '#define' and not ignore:
                # Add a value to our list of defines.
                if len(fields) < 2:
                    raise GromacsError('Illegal line in .top file: '+line)
                name = fields[1]
                value_start = stripped.find(name, len(command))+len(name)+1
                value = line[value_start:].strip()
                self.defines[name] = value
            elif command == '#ifdef':
                # See whether this block should be ignored.
                if len(fields) < 2:
                    raise GromacsError('Illegal line in .top file: '+line)
                name = fields[1]
                self.if_stack.append(name in self.defines)
                self.else_stack.append(False)
            elif command == '#ifndef':
                # See whether this block should be ignored.
                if len(fields) < 2:
                    raise GromacsError('Illegal line in .top file: '+line)
                name = fields[1]
                self.if_stack.append(name not in self.defines)
                self.else_stack.append(False)
            elif command == '#endif':
                # Pop an entry off the if stack
                if len(self.if_stack) == 0:
                    raise GromacsError('Unexpected line in .top file: '+line)
                del(self.if_stack[-1])
                del(self.else_stack[-1])
            elif command == '#else':
                # Reverse the last entry on the if stack
                if len(self.if_stack) == 0:
                    raise GromacsError('Unexpected line in .top file: '+line)
                if self.else_stack[-1]:
                    raise GromacsError('Unexpected line in .top file: #else has'
                                     ' already been used ' + line)
                self.if_stack[-1] = (not self.if_stack[-1])
                self.else_stack[-1] = True

        elif not ignore:
            # A line of data for the current category
            if self.current_directive is None:
                raise GromacsError('Unexpected line in .top file: "{0}"'.format(line))
            if self.current_directive == 'defaults':
                self.process_defaults(line)
            elif self.current_directive == 'moleculetype':
                self.process_moleculetype(line)
            elif self.current_directive == 'molecules':
                self.process_molecule(line)
            elif self.current_directive == 'atoms':
                self.process_atom(line)
            elif self.current_directive == 'bonds':
                self.process_bond(line)
            elif self.current_directive == 'angles':
                self.process_angle(line)
            elif self.current_directive == 'dihedrals':
                self.process_dihedral(line)
            elif self.current_directive == 'settles':
                self.process_settle(line)
            elif self.current_directive == 'exclusions':
                self.process_exclusion(line)
            elif self.current_directive == 'pairs':
                self.process_pair(line)
            elif self.current_directive == 'cmap':
                self.process_cmap(line)
            elif self.current_directive == 'atomtypes':
                self.process_atomtype(line)
            elif self.current_directive == 'bondtypes':
                self.process_bondtype(line)
            elif self.current_directive == 'angletypes':
                self.process_angletype(line)
            elif self.current_directive == 'dihedraltypes':
                self.process_dihedraltype(line)
            elif self.current_directive == 'implicit_genborn_params':
                self.process_implicittype(line)
            elif self.current_directive == 'pairtypes':# and not self.system.genpairs:
                self.process_pairtype(line)
            elif self.current_directive == 'cmaptypes':
                self.process_cmaptype(line)
            elif self.current_directive == 'nonbond_params':
                self.process_nonbond_params(line)
            elif self.current_directive.startswith('virtual_sites'):
                vsite_type = self.current_directive[-1]
                self.process_virtual_sites(line, vsite_type)

    def process_defaults(self, line):
        """Process the [ defaults ] line."""
        fields = line.split()
        if len(fields) < 4:
            self.too_few_fields(line)
        self.system.nonbonded_function = int(fields[0])
        self.system.combination_rule = self.gromacs_combination_rules[fields[1]]
        self.system.genpairs = fields[2]
        self.system.lj_correction = float(fields[3])
        self.system.coulomb_correction = float(fields[4])

    def process_moleculetype(self, line):
        """Process a line in the [ moleculetypes ] category."""
        fields = line.split()
        if len(fields) < 1:
            self.too_few_fields(line)
        mol_type = self.TopMoleculeType()
        mol_type.nrexcl = int(fields[1])
        self.molecule_types[fields[0]] = mol_type
        self.current_molecule_type = mol_type

    def process_molecule(self, line):
        """Process a line in the [ molecules ] category."""
        fields = line.split()
        if len(fields) < 2:
            self.too_few_fields(line)
        self.molecules.append((fields[0], int(fields[1])))

    def process_atom(self, line):
        """Process a line in the [ atoms ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 5:
            self.too_few_fields(line)
        if len(fields) not in [7, 8, 11]:
            self.invalid_line(line)
        self.current_molecule_type.atoms.append(fields)

    def process_bond(self, line):
        """Process a line in the [ bonds ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 3:
            self.too_few_fields(line)
        self.current_molecule_type.bonds.append(fields)

    def process_angle(self, line):
        """Process a line in the [ angles ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 4:
            self.too_few_fields(line)
        self.current_molecule_type.angles.append(fields)

    def process_dihedral(self, line):
        """Process a line in the [ dihedrals ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 5:
            self.too_few_fields(line)
        self.current_molecule_type.dihedrals.append(fields)

    def process_settle(self, line):
        """Process a line in the [ settles ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 4:
            self.too_few_fields(line)
        self.current_molecule_type.rigidwaters.append(fields)

    def process_exclusion(self, line):
        """Process a line in the [ exclusions ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 2:
            self.too_few_fields(line)
        self.current_molecule_type.exclusions.append(fields)

    def process_pair(self, line):
        """Process a line in the [ pairs ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 3:
            self.too_few_fields(line)
        self.current_molecule_type.pairs.append(fields)

    def process_cmap(self, line):
        """Process a line in the [ cmaps ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype('cmap')
        fields = line.split()
        if len(fields) < 6:
            self.too_few_fields(line)
        self.current_molecule_type.cmaps.append(fields)

    def process_atomtype(self, line):
        """Process a line in the [ atomtypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            self.too_few_fields(line)
        if len(fields[3]) == 1:
            # Bonded type and atomic number are both missing.
            fields.insert(1, None)
            fields.insert(1, None)
        elif len(fields[4]) == 1 and len(fields[5]) >= 1:
            if fields[1][0].isalpha():
                # Atomic number is missing.
                fields.insert(2, None)
            else:
                # Bonded type is missing.
                fields.insert(1, None)

        atomtype = fields[0]
        if fields[1] == None:
            bondingtype = atomtype
        else:
            bondingtype = fields[1]
        if fields[2]:
            atomic_number = int(fields[2])
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
            raise InterMolError("Unknown combination rule: {0}".format(self.system.combination_rule))
        new_atom_type = AtomtypeClass(atomtype, bondingtype, atomic_number,
                                      mass, charge, ptype, lj_param1, lj_param2)
        self.system.add_atomtype(new_atom_type)

    def process_bondtype(self, line):
        """Process a line in the [ bondtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            self.too_few_fields(line)

        btypes = fields[:2]
        bond_type = self.process_forcetype(btypes, 'bond', line, 2,
                self.gromacs_bond_types, self.canonical_bond)
        self.bondtypes[tuple(fields[:2])] = bond_type

    def process_angletype(self, line):
        """Process a line in the [ angletypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            self.too_few_fields(line)
        btypes = fields[:3]
        angle_type = self.process_forcetype(btypes, 'angle', line, 3,
                self.gromacs_angle_types, self.canonical_angle)
        self.angletypes[tuple(fields[:3])] = angle_type

    def process_dihedraltype(self, line):
        """Process a line in the [ dihedraltypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            self.too_few_fields(line)

        # Some gromacs parameters don't include sufficient numbers of types.
        # Add some zeros (bit of a kludge).
        line += ' 0.0 0.0 0.0'
        fields = line.split()

        # Check whether they are using 2 or 4 atom types
        if fields[2].isdigit():
            btypes = ['X', fields[0], fields[1], 'X']
            n_atoms_specified = 2
        elif fields[4].isdigit() and not fields[3].isdigit(): # assumes gromacs types are not all digits.
            btypes = fields[:4]
            n_atoms_specified = 4
        else:
            # TODO: Come up with remaining cases (are there any?) and a proper
            #       failure case.
            logger.warning('Should never have gotten here.')
        dihedral_type = self.process_forcetype(
            btypes, 'dihedral', line, n_atoms_specified,
            self.gromacs_dihedral_types, self.canonical_dihedral)

        # Still need a bit more information
        numeric_dihedraltype = fields[n_atoms_specified]
        dihedral_type.improper = numeric_dihedraltype in ['2', '4']

        key = tuple([btypes[0], btypes[1], btypes[2], btypes[3],
                     dihedral_type.improper])

        if key in self.dihedraltypes:
            # There are multiple dihedrals defined for these atom types.
            self.dihedraltypes[key].add(dihedral_type)
        else:
            self.dihedraltypes[key] = {dihedral_type}

    def process_forcetype(self, bondingtypes, forcename, line, n_atoms,
                          gromacs_force_types, canonical_force):
        """ """
        fields = line.split()

        numeric_forcetype = fields[n_atoms]
        gromacs_force_type = gromacs_force_types[numeric_forcetype]
        kwds = self.create_kwds_from_entries(fields, gromacs_force_type, offset=n_atoms+1)
        CanonicalForceType, kwds = canonical_force(
            kwds, gromacs_force_type, direction='into')

        force_type = CanonicalForceType(*bondingtypes, **kwds)

        if not force_type:
            logger.warning("{0} is not a supported {1} type".format(fields[2], forcename))
            return
        else:
            return force_type

    def process_implicittype(self, line):
        """Process a line in the [ implicit_genborn_params ] category."""
        fields = line.split()
        if len(fields) < 6:
            self.too_few_fields(line)
        self.implicittypes[fields[0]] = fields

    def process_pairtype(self, line):
        """Process a line in the [ pairtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            self.too_few_fields(line)

        pair_type = None
        PairFunc = None
        combination_rule = self.system.combination_rule
        kwds = dict()
        numeric_pairtype = fields[2]
        if numeric_pairtype == '1':
            # LJ/Coul. 1-4 (Type 1)
            if len(fields) == 5:
                if combination_rule == "Multiply-C6C12":
                    PairFunc = LjCPairType
                elif combination_rule in ['Multiply-Sigeps', 'Lorentz-Berthelot']:
                    PairFunc = LjSigepsPairType
            offset = 3
        elif numeric_pairtype == '2':
            if combination_rule == "Multiply-C6C12":
                PairFunc = LjqCPairType
            elif combination_rule in ['Multiply-Sigeps', 'Lorentz-Berthelot']:
                PairFunc = LjqSigepsPairType
            offset = 4
        else:
            logger.warning("Could not find pair type for line: {0}".format(line))

        if PairFunc:
            pairvars = [fields[0], fields[1]]
            kwds = self.create_kwds_from_entries(fields, PairFunc, offset=offset)
            # kludge because of placement of scaleQQ...
            if numeric_pairtype == '2':
                # try to get this out ...
                kwds['scaleQQ'] = float(fields[3]) * units.dimensionless
            pair_type = PairFunc(*pairvars, **kwds)

        self.pairtypes[tuple(fields[:2])] = pair_type

    def process_cmaptype(self, line):
        """Process a line in the [ cmaptypes ] category."""
        fields = line.split()
        if len(fields) < 8 or len(fields) < 8+int(fields[6])*int(fields[7]):
            self.too_few_fields(line)
        self.cmaptypes[tuple(fields[:5])] = fields

    def process_nonbond_params(self, line):
        """Process a line in the [ nonbond_param ] category."""
        fields = line.split()
        NonbondedFunc = None
        combination_rule = self.system.combination_rule

        if fields[2] == '1':
            if combination_rule == 'Multiply-C6C12':
                NonbondedFunc = LjCNonbondedType
            elif combination_rule in ['Lorentz-Berthelot', 'Multiply-Sigeps']:
                NonbondedFunc = LjSigepsNonbondedType
        elif fields[2] == '2':
            if combination_rule == 'Buckingham':
                NonbondedFunc = BuckinghamNonbondedType
        else:
            logger.warning("Could not find nonbonded type for line: {0}".format(line))

        nonbonded_vars = [fields[0], fields[1]]
        kwds = self.create_kwds_from_entries(fields, NonbondedFunc, offset=3)
        nonbonded_type = NonbondedFunc(*nonbonded_vars, **kwds)
        # TODO: figure out what to do with the gromacs numeric type
        nonbonded_type.type = int(fields[2])
        self.system.nonbonded_types[tuple(nonbonded_vars)] = nonbonded_type

    def process_virtual_sites(self, line, v_site_type):
        """Process a line in a [ virtual_sites? ] category."""
        if v_site_type == 'n':
            raise UnimplementedSetting('Parsing of [ virtual_sitesn ] directives'
                                       ' is not yet implemented')
        fields = line.split()
        self.current_molecule_type.virtuals[v_site_type].append(fields)

    # =========== Pre-processing errors =========== #
    def too_few_fields(self, line):
        raise GromacsError('Too few fields in [ {0} ] line: {1}'.format(
            self.current_directive, line))

    def invalid_line(self, line):
        raise GromacsError('Invalid format in [ {0} ] line: {1}'.format(
            self.current_directive, line))

    def directive_before_moleculetype(self):
        raise GromacsError('Found [ {0} ] directive before [ moleculetype ]'.format(
            self.current_directive))
