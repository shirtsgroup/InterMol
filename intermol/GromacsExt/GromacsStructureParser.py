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
                    atom.residueIndex=int(lines[i][0:5])
                    atom.residueName = lines[i][6:10].strip()
                    atom.atomName = lines[i][11:15].strip()
                    atom.atomIndex = int(lines[i][16:20])
                    variables = (lines[i][20:]).split()
                    position = np.zeros([3],float) * units.nanometers
                    velocity = np.zeros([3],float) * units.nanometers / units.picoseconds
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
    v = np.zeros([3,3],float) * units.nanometers
    if len(rawBoxVector) == 3:
        for i in range(3):
            v[i,i] = float(rawBoxVector[i]) * units.nanometers

    elif len(rawBoxVector) == 9:
        n = -1
        BoxVectorElements = split.rawBoxVector()
        for i in range(3):
            for j in range(3):
                v[i,j] = float(BoxVectorElements[3*i+j]) * units.nanometers
    # need to make this numpy sensitive so we can just pass the vector
    System._sys.setBoxVector(v[0,0],v[0,1],v[0,2],v[1,0],v[1,1],v[1,2],v[2,0],v[2,1],v[2,2])

def writeStructure(filename):
    """Write the system out  in a Gromacs 4.5.4 format

    Args:
        filename (str): the file to write out to
    """
    lines = list()
    i=0
    for moleculetype in System._sys._molecules.values():
        for molecule in moleculetype.moleculeSet:
            for atom in molecule._atoms:
                i += 1
                lines.append('%5d%-4s%6s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n'
                        %(atom.residueIndex,
                        atom.residueName,
                        atom.atomName,
                        atom.atomIndex,
                        atom._position[0]._value,
                        atom._position[1]._value,
                        atom._position[2]._value,
                        atom._velocity[0]._value,
                        atom._velocity[1]._value,
                        atom._velocity[2]._value))
    lines.append('%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n'
          %(System._sys._v1x._value,
            System._sys._v2y._value,
            System._sys._v3z._value,
            System._sys._v1y._value,
            System._sys._v1z._value,
            System._sys._v2x._value,
            System._sys._v2z._value,
            System._sys._v3x._value,
            System._sys._v3y._value))

    lines.insert(0, (System._sys._name+'\n'))
    lines.insert(1, str(i)+'\n')

    fout = open(filename, 'w')
    for line in lines:
        fout.write(line)
    fout.close()


