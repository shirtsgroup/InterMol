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
    """
    if len(gbb.pair_types) > 0:
        f.write('\n')
        f.write('Pair Coeffs\n')
        f.write('\n')
        for i, value in sorted(gbb.pair_types.items()):
            f.write("%d %8.4f %8.4f\n" % (i, value[0], value[1]))

    if len(gbb.bond_types) > 0:
        f.write('\n')
        f.write('Bond Coeffs\n')
        f.write('\n')
        for i, value in sorted(gbb.bond_types.items()):
            f.write("%d %8.4f %8.4f\n" % (i, value[0], value[1]))

    if len(gbb.angle_types) > 0:
        f.write('\n')
        f.write('Angle Coeffs\n')
        f.write('\n')
        for i, value in sorted(gbb.angle_types.items()):
            f.write("%d %8.4f %8.4f\n" % (i, value[0], value[1]))

    if len(gbb.dihedral_types) > 0:
        f.write('\n')
        f.write('Dihedral Coeffs\n')
        f.write('\n')
        for i, value in sorted(gbb.dihedral_types.items()):
            f.write("%d %8.4f %8.4f %8.4f %8.4f\n" % (i, value[0], value[1], value[2], value[3]))
    """

    n_atom_types = 0
    n_bond_types = 0
    n_angle_types = 0
    n_dihedral_types = 0
    n_improper_types = 0

    # add mapping for numeric types
    atom_list = list()
    atom_list.append('\n')
    atom_list.append('Atoms\n')
    atom_list.append('\n')

    vel_list = list()
    vel_list.append('\n')
    vel_list.append('Velocities\n')
    vel_list.append('\n')
    for moleculetype in System._sys._molecules.values():
        for molecule in moleculetype.moleculeSet:
            for atom in molecule._atoms:
                atom_list.append('%-6d %-6d %-6d %5.8f %8.3f %8.3f %8.3f\n'
                    % (atom.atomIndex,
                       atom.residueIndex,
                       #atom.numericType,
                       1,
                       atom._charge[0]._value,
                       atom._position[0]._value,
                       atom._position[1]._value,
                       atom._position[2]._value))
                vel_list.append('%-6d %8.4f %8.4f %8.4f\n'
                    % (atom.atomIndex,
                       atom._velocity[0]._value,
                       atom._velocity[1]._value,
                       atom._velocity[2]._value))

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

    with open(filename, 'w') as f:
        f.write(System._sys._name + '\n')
        f.write('\n')

        n_atoms = len(atom_list)
        n_bonds = len(bond_list)
        n_angles = len(angle_list)
        n_dihedrals = len(dihedral_list)
        n_impropers = len(improper_list)

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

        f.write('\n')
        f.write('Masses\n')
        f.write('\n')

        # find unique masses and corresponding atomtypes
        masses = set()
        for moleculetype in System._sys._molecules.values():
            for molecule in moleculetype.moleculeSet:
                for atom in molecule._atoms:
                    #masses.add((atom.numericType, atom._mass))
                    masses.add((1, atom._mass[0]._value))
        for mass in sorted(masses):
            f.write(" ".join(map(str, mass)) + '\n')

        for atom in atom_list:
            f.write(atom)

        for vel in vel_list:
            f.write(vel)

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
