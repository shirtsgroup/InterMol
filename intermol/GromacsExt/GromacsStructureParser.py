import pdb
import numpy as np
from intermol.Atom import *
from intermol.Molecule import *
from intermol.System import System


def readStructure(filename):
    """Read in a Gromacs structure file

    Args:
        filename (str): the file to be read in
    """

    fd = open(filename, 'r')
    lines = list(fd)
    fd.close()

    i = 2
    for moleculetype in System._sys._molecules.values():
        for molecule in moleculetype.moleculeSet:
            for atom in molecule._atoms:
                if lines[i]:
                    atom.residueIndex = int(lines[i][0:5])
                    atom.residueName = lines[i][5:10].strip()
                    atom.atomName = lines[i][10:15].strip()
                    atom.atomIndex = int(lines[i][15:20])
                    variables = (lines[i][20:]).split()
                    position = np.zeros([3], float) * units.nanometers
                    velocity = np.zeros([3], float) * units.nanometers / units.picoseconds
                    if len(variables) >= 3:
                        positionElements = variables[0:3]
                        for k in range(3):
                            position[k] = float(positionElements[k]) * units.nanometers
                    atom.setPosition(position[0], position[1], position[2])
                    if len(variables) >= 6:
                        velocityElements = variables[3:6]
                        for k in range(3):
                            velocity[k] = float(velocityElements[k]) * units.nanometers / units.picoseconds
                    atom.setVelocity(velocity[0], velocity[1], velocity[2])
                    i += 1
                else:
                    sys.exit()

    rawBoxVector = lines[i].split()
    v = np.zeros([3, 3], float) * units.nanometers
    if len(rawBoxVector) == 3:
        for i in range(3):
            v[i, i] = float(rawBoxVector[i]) * units.nanometers

    elif len(rawBoxVector) == 9:
        for i in range(3):
            for j in range(3):
                v[i, j] = float(rawBoxVector[3*i+j]) * units.nanometers
    # need to make this numpy sensitive so we can just pass the vector
    System._sys.setBoxVector(v)


def writeStructure(filename):
    """Write the system out  in a Gromacs 4.6 format

    Args:
        filename (str): the file to write out to
    """
    lines = list()
    n = 0
    for moleculetype in System._sys._molecules.values():
        for molecule in moleculetype.moleculeSet:
            for atom in molecule._atoms:
                if atom.atomName.isdigit():
                    atom.atomName = "LMP_" + atom.atomName
                n += 1
                lines.append('%5d%-4s%6s%5d%13.8f%13.8f%13.8f%13.8f%13.8f%13.8f\n'
                             % (atom.residueIndex,
                                atom.residueName,
                                atom.atomName,
                                atom.atomIndex,
                                atom._position[0]._value,
                                atom._position[1]._value,
                                atom._position[2]._value,
                                atom._velocity[0]._value,
                                atom._velocity[1]._value,
                                atom._velocity[2]._value))
    # print the box vector
    # check for rectangular; should be symmetric, so we don't have to check 6 values
    if (System._sys._boxVector[0, 1]._value == 0 and
        System._sys._boxVector[0, 2]._value == 0 and 
        System._sys._boxVector[1, 2]._value == 0):
            for i in range(3):
                lines.append('%11.7f ' % System._sys._boxVector[i, i].in_units_of(units.nanometers)._value)
    else:
        for i in range(3):
            for j in range(3):
                lines.append('%11.7f ' % System._sys._boxVector[i, j].in_units_of(units.nanometers)._value)
    lines.append('\n')

    lines.insert(0, (System._sys._name + '\n'))
    lines.insert(1, str(n) + '\n')

    fout = open(filename, 'w')
    for line in lines:
        fout.write(line)
    fout.close()
