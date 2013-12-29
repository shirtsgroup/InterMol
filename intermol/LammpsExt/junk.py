import warnings
import re
import pdb

import intermol.unit as units
from intermol.System import System
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
    with open(filename, 'w') as f:
        f.write(System._sys._name + '\n')
        f.write('\n')

        n_atoms = 0
        n_bonds = 0
        n_angles = 0
        n_dihedrals = 0
        n_impropers = 0

        f.write(str(gbb.xyz.shape[0]) + ' atoms\n')
        f.write(str(n_bonds) + ' bonds\n')
        f.write(str(n_angles) + ' angles\n')
        f.write(str(n_dihedrals) + ' dihedrals\n')
        f.write('\n')

        f.write(str(int(gbb.types.max())) + ' atom types\n')
        if n_bonds > 0:
            f.write(str(int(gbb.bonds[:, 0].max())) + ' bond types\n')
        if n_angles > 0:
            f.write(str(int(gbb.angles[:, 0].max())) + ' angle types\n')
        if n_dihedrals > 0:
            f.write(str(int(gbb.dihedrals[:, 0].max())) + ' dihedral types\n')
        f.write('\n')

        f.write('0.0 %8.4f xlo xhi\n' % (System._sys._v1x._value))
        f.write('0.0 %8.4f ylo yhi\n' % (System._sys._v2y._value))
        f.write('0.0 %8.4f zlo zhi\n' % (System._sys._v3z._value))

        f.write('\n')
        f.write('Masses\n')
        f.write('\n')

        # find unique masses and corresponding atomtypes
        masses = set()
        for i, mass in enumerate(gbb.masses):
            masses.add((int(gbb.types[i]), mass))
        for mass in sorted(masses):
            f.write(" ".join(map(str, mass)) + '\n')

        for moleculetype in System.__sys.molecules.values():
            for molecule in moleculetype.moleculeSet:
                for atom in molecule._atoms:
                    masses.add((atom.numericType, atom._mass))

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

        # add mapping for numeric types
        f.write('\n')
        f.write('Atoms\n')
        f.write('\n')
        for moleculetype in System.__sys.molecules.values():
            for molecule in moleculetype.moleculeSet:
                for atom in molecule._atoms:
                    f.write('%-6d %-6d %-6d %5.8f %8.3f %8.3f %8.3f\n'
                        % (atom.atomIndex,
                           atom.residueIndex,
                           atom.numericType,
                           atom._charge[0]._value,
                           atom._position[0]._value,
                           atom._position[1]._value,
                           atom._position[2]._value))

        # velocities

        if n_bonds > 0:
            f.write('\n')
            f.write('Bonds\n')
            f.write('\n')
            for i, bond in enumerate(gbb.bonds):
                f.write(str(i+1) + " " + " ".join(map(str, bond)) + '\n')

        if n_angles > 0:
            f.write('\n')
            f.write('Angles\n')
            f.write('\n')
            for i, angle in enumerate(gbb.angles):
                f.write(str(i+1) + " " + " ".join(map(str, angle)) + '\n')

        if n_dihedrals > 0:
            f.write('\n')
            f.write('Dihedrals\n')
            f.write('\n')
            for i, dihedral in enumerate(gbb.dihedrals):
                f.write(str(i+1) + " " + " ".join(map(str, dihedral)) + '\n')

        if n_impropers > 0:
            if prototype:
                f_dihedral = open(sys_name + '_dihedral.txt', 'w')
            f.write('\n')
            f.write('Impropers\n')
            f.write('\n')
            for i, improper in enumerate(gbb.impropers):
                f.write(str(i+1) + " " + " ".join(map(str, improper)) + '\n')

    print "Wrote file '" + filename + "'"
