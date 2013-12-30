from os.path import splitext
import warnings
import re
import pdb

import intermol.unit as units
from intermol.Atom import Atom
from intermol.Molecule import Molecule
from intermol.System import System
from intermol.Types import *
from intermol.Force import *
from intermol.HashMap import *

def writeData(filename):
    """Write system to LAMMPS data file

    Args:
        filename (str): name of output file


    """
    # Real units
    # TODO: allow user to specify different unit sets
    LENGTH = units.angstroms
    ENERGY = units.kilocalorie
    MASS = units.grams / units.mole
    MOLE = units.mole

    mass_list = list()
    mass_list.append('\n')
    mass_list.append('Masses\n')
    mass_list.append('\n')

    pair_coeff_list = list()
    pair_coeff_list.append('\n')
    pair_coeff_list.append('Pair Coeffs\n')
    pair_coeff_list.append('\n')

    bond_coeff_list = list()
    bond_coeff_list.append('\n')
    bond_coeff_list.append('Bond Coeffs\n')
    bond_coeff_list.append('\n')

    angle_coeff_list = list()
    angle_coeff_list.append('\n')
    angle_coeff_list.append('Angle Coeffs\n')
    angle_coeff_list.append('\n')

    dihedral_coeff_list = list()
    dihedral_coeff_list.append('\n')
    dihedral_coeff_list.append('Dihedral Coeffs\n')
    dihedral_coeff_list.append('\n')

    improper_coeff_list = list()
    improper_coeff_list.append('\n')
    improper_coeff_list.append('Improper Coeffs\n')
    improper_coeff_list.append('\n')

    atom_list = list()
    atom_list.append('\n')
    atom_list.append('Atoms\n')
    atom_list.append('\n')

    vel_list = list()
    vel_list.append('\n')
    vel_list.append('Velocities\n')
    vel_list.append('\n')

    atom_type_dict = dict()
    a_i = 1  # counter for atom types
    for moleculeType in System._sys._molecules.itervalues():
        for molecule in moleculeType.moleculeSet:
            for atom in molecule._atoms:
                # type, mass and pair coeffs
                if atom._atomtype[0] in atom_type_dict:
                    pass
                else:
                    atom_type_dict[atom._atomtype[0]] = a_i
                    mass_list.append('%d %8.4f\n'
                                % (a_i,
                                   atom._mass[0].in_units_of(MASS)._value))
                    pair_coeff_list.append('%d %8.4f %8.4f\n'
                                % (a_i,
                                   atom._sigma[0].in_units_of(LENGTH)._value,
                                   atom._epsilon[0].in_units_of(ENERGY/MOLE)._value))
                    a_i += 1
                # atom
                atom_list.append('%-6d %-6d %-6d %5.8f %8.3f %8.3f %8.3f\n'
                    % (atom.atomIndex,
                       atom.residueIndex,
                       atom_type_dict[atom._atomtype[0]],
                       atom._charge[0]._value,
                       atom._position[0]._value,
                       atom._position[1]._value,
                       atom._position[2]._value))
                # velocity
                vel_list.append('%-6d %8.4f %8.4f %8.4f\n'
                    % (atom.atomIndex,
                       atom._velocity[0]._value,
                       atom._velocity[1]._value,
                       atom._velocity[2]._value))

    bond_style = 0
    bond_type_set = set()
    bond_type_dict = dict()
    for moleculeType in System._sys._molecules.itervalues():
        # bond types
        if moleculeType.bondForceSet:
            for bond in moleculeType.bondForceSet.itervalues():
                if isinstance(bond, Bond):
                    if bond_style == 0:
                        bond_style = 'harmonic'
                    elif bond_style != 'harmonic':
                        warnings.warn("More than one bond style found!")
                else:
                    warnings.warn("Found unsupported bond type for LAMMPS!")

                identifier = 0
                bond_type_set.add(identifier)
                temp = (type(bond),
                        bond.length._value,
                        bond.k._value)
                bond_type_dict[identifier] = temp

    bond_type_int = dict(zip(bond_type_set, range(1, len(bond_type_set) + 1)))
    pdb.set_trace()

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

    for moleculeType in System._sys._molecules.itervalues():
        if moleculeType.bondForceSet:
            for bond in moleculeType.bondForceSet.itervalues():
                if isinstance(bond, Bond):
                    b_type = 1  # LOOKUP NUMBER
                    line = ('%4d%7d%7d\n'
                            % (b_type,
                               bond.atom1,
                               bond.atom2))

                bond_list.append(line)
        if moleculeType.angleForceSet:
            for angle in moleculeType.angleForceSet.itervalues():
                if isinstance(angle, Angle):
                    a_type = 1  # LOOKUP NUMBER
                    line = ('%4d%7d%7d%7d\n'
                            % (a_type,
                               angle.atom1,
                               angle.atom2,
                               angle.atom3))

                angle_list.append(line)

    # actual file writing
    with open(filename, 'w') as f:
        f.write(System._sys._name + '\n')
        f.write('\n')

        n_atoms = len(atom_list) - 3
        n_bonds = len(bond_list) - 3
        n_angles = len(angle_list) - 3
        n_dihedrals = len(dihedral_list) - 3
        n_impropers = len(improper_list) - 3

        n_atom_types = len(pair_coeff_list) - 3
        n_bond_types = len(bond_coeff_list) - 3
        n_angle_types = len(angle_coeff_list) - 3
        n_dihedral_types = len(dihedral_coeff_list) - 3
        n_improper_types = len(improper_coeff_list) - 3

        f.write(str(n_atoms) + ' atoms\n')
        f.write(str(n_bonds) + ' bonds\n')
        f.write(str(n_angles) + ' angles\n')
        f.write(str(n_dihedrals) + ' dihedrals\n')
        f.write(str(n_impropers) + ' impropers\n')
        f.write('\n')

        f.write(str(n_atom_types) + ' atom types\n')
        if n_bonds > 0:
            f.write(str(n_bond_types) + ' bond types\n')
        if n_angles > 0:
            f.write(str(n_angle_types) + ' angle types\n')
        if n_dihedrals > 0:
            f.write(str(n_dihedral_types) + ' dihedral types\n')
        if n_impropers > 0:
            f.write(str(n_improper_types) + ' improper types\n')
        f.write('\n')

        f.write('0.0 %10.6f xlo xhi\n' % (System._sys._v1x._value))
        f.write('0.0 %10.6f ylo yhi\n' % (System._sys._v2y._value))
        f.write('0.0 %10.6f zlo zhi\n' % (System._sys._v3z._value))

        # write masses
        for mass in mass_list:
            f.write(mass)

        # write coefficients
        for pair in pair_coeff_list:
            f.write(pair)
        if len(bond_coeff_list) > 3:
            for bond in bond_coeff_list:
                f.write(bond)
        if len(angle_coeff_list) > 3:
            for angle in angle_coeff_list:
                f.write(pair)
        if len(dihedral_coeff_list) > 3:
            for dihedral in dihedral_coeff_list:
                f.write(dihedral)
        if len(improper_coeff_list) > 3:
            for improper in improper_coeff_list:
                f.write(improper)

        # write atoms and velocities
        for atom in atom_list:
            f.write(atom)
        for vel in vel_list:
            f.write(vel)

        # write connectivity
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

    basename = splitext(filename)[0]
    input_filename = basename + '.input'
    with open(input_filename, 'w') as f:
        # units
        # atomstyle
        # pair style
        # pair modify
        # bond style
        # angle style
        # dihedral style
        # improper style
        # special bonds?
        pass
