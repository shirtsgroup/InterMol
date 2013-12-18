import sys
import os
import string
from collections import deque

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
                    atom.residueName=lines[i][6:10].strip()
                    atom.atomName=lines[i][11:15].strip()
                    atom.atomIndex=int(lines[i][16:20])
                    position = [None, None, None]
                    position[0]=float(lines[i][20:28]) * units.nanometers
                    position[1]=float(lines[i][28:36]) * units.nanometers
                    position[2]=float(lines[i][36:44]) * units.nanometers
                    atom.setPosition(position[0], position[1], position[2])
                    velocity  = [None, None, None]
                    try:
                        velocity[0]=float(lines[i][44:52]) * units.nanometers / units.picoseconds
                        velocity[1]=float(lines[i][52:60]) * units.nanometers / units.picoseconds
                        velocity[2]=float(lines[i][60:68]) * units.nanometers / units.picoseconds
                    except:
                        velocity[0] = 0.0 * units.nanometers / units.picoseconds
                        velocity[1] = 0.0 * units.nanometers / units.picoseconds
                        velocity[2] = 0.0 * units.nanometers / units.picoseconds
                    atom.setVelocity(velocity[0], velocity[1], velocity[2])
                    i += 1

                else:
                    sys.exit()

    rawBoxVector = lines[i]
    v1x = None
    v2x = None
    v3x = None
    v1y = None
    v2y = None
    v3y = None
    v1z = None
    v2z = None
    v3z = None
    if len(rawBoxVector.split()) == 3:
        v1x = float(rawBoxVector[0:10]) * units.nanometers
        v2x = 0.0 * units.nanometers
        v3x = 0.0 * units.nanometers
        v1y = 0.0 * units.nanometers
        v2y = float(rawBoxVector[11:22]) * units.nanometers
        v3y = 0.0 * units.nanometers
        v1z = 0.0 * units.nanometers
        v2z = 0.0 * units.nanometers
        v3z = float(rawBoxVector[23:34]) * units.nanometers

    elif len(rawBoxVector.split()) == 9:
        v1x = float(rawBoxVector[0:10]) * units.nanometers
        v2x = float(rawBoxVector[11:22]) * units.nanometers
        v3x = float(rawBoxVector[23:34]) * units.nanometers
        v1y = float(rawBoxVector[35:46]) * units.nanometers
        v2y = float(rawBoxVector[47:58]) * units.nanometers
        v3y = float(rawBoxVector[59:70]) * units.nanometers
        v1z = float(rawBoxVector[71:82]) * units.nanometers
        v2z = float(rawBoxVector[83:94]) * units.nanometers
        v3z = float(rawBoxVector[95:106]) * units.nanometers
    System._sys.setBoxVector(v1x, v2x, v3x, v1y, v2y, v3y, v1z, v2z, v3z)

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


