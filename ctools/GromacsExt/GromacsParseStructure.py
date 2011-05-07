import pdb
import sys
import os
import string
import re
from collections import deque
from ctools.Atom import *
from ctools.Molecule import *



def readStructure(filename, system):
    """
    Read in a Gromacs structure file
    """
    
    fd = open(filename, 'r')
    lines = list(fd)
    fd.close()
    
    i = 2
    for moleculetype in system.molecules.values():
        for molecule in moleculetype.moleculeSet:
            for atom in molecule.atoms:
                if lines[i]:
                    atom.resNum=int(lines[i][0:5])
                    atom.resName=lines[i][6:10].strip()
                    atom.atomName=lines[i][11:15].strip()
                    atom.atomNum=int(lines[i][16:20])
                    atom.position[0]=float(lines[i][20:28])*units.nanometers
                    atom.position[1]=float(lines[i][28:36])*units.nanometers
                    atom.position[2]=float(lines[i][36:44])*units.nanometers
                    try:
                        atom.velocity[0]=float(lines[i][44:52])*units.nanometers * units.picoseconds**(-1)
                        atom.velocity[1]=float(lines[i][52:60])*units.nanometers * units.picoseconds**(-1)
                        atom.velocity[2]=float(lines[i][60:68])*units.nanometers * units.picoseconds**(-1)
                    except:
                        atom.velocity[0] = 0.0*units.nanometers * units.picoseconds**(-1)
                        atom.velocity[1] = 0.0*units.nanometers * units.picoseconds**(-1)
                        atom.velocity[2] = 0.0*units.nanometers * units.picoseconds**(-1)

                    i += 1

                else:
                    sys.exit()
 
    rawBoxVector = lines[i]      
    if len(rawBoxVector.split()) == 3:
        system.v1x = float(rawBoxVector[0:10]) * units.nanometers
        system.v2x = 0.0 * units.nanometers
        system.v3x = 0.0 * units.nanometers
        system.v1y = 0.0 * units.nanometers
        system.v2y = float(rawBoxVector[11:22]) * units.nanometers
        system.v3y = 0.0 * units.nanometers
        system.v1z = 0.0 * units.nanometers
        system.v2z = 0.0 * units.nanometers
        system.v3z = float(rawBoxVector[23:34]) * units.nanometers

    if len(rawBoxVector.split()) == 9:
        system.v1x = float(rawBoxVector[0:10]) * units.nanometers
        system.v2x = float(rawBoxVector[11:22]) * units.nanometers
        system.v3x = float(rawBoxVector[23:34]) * units.nanometers
        system.v1y = float(rawBoxVector[35:46]) * units.nanometers
        system.v2y = float(rawBoxVector[47:58]) * units.nanometers
        system.v3y = float(rawBoxVector[59:70]) * units.nanometers
        system.v1z = float(rawBoxVector[71:82]) * units.nanometers
        system.v2z = float(rawBoxVector[83:94]) * units.nanometers
        system.v3z = float(rawBoxVector[95:106]) * units.nanometers

def writeStructure(system, filename):

    lines = list()

    i=0
    for moleculetype in system.molecules.values():
        for molecule in moleculetype.moleculeSet:
            for atom in molecule.atoms:
                i += 1
                lines.append('%5d%-4s%6s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n'
                        %(atom.resNum,
                        atom.resName,
                        atom.atomName,
                        atom.atomNum,
                        atom.position[0]._value,
                        atom.position[1]._value,
                        atom.position[2]._value,
                        atom.velocity[0]._value,
                        atom.velocity[1]._value,
                        atom.velocity[2]._value))
    lines.append('%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n'
          %(system.v1x._value,
            system.v2y._value,
            system.v3z._value,
            system.v1y._value,
            system.v1z._value,
            system.v2x._value,
            system.v2z._value,
            system.v3x._value,
            system.v3y._value))

    lines.insert(0, (system.name+'\n'))
    lines.insert(1, str(i)+'\n')



    fout = open(filename, 'w')
    for line in lines:
        fout.write(line)
    fout.close()


