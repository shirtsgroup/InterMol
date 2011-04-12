import pdb
import sys
import os
import string
import re
from collections import deque
from Topology.Atom import *
from Topology.Molecule import *



def readStructure(filename):
    """
    Read in a Gromacs structure file
    """
    
    def readLines(filename):
        fin = open(filename, 'r')
        lines = fin.readlines()
        fin.close()
        return lines

    def preprocess(lines):

        processedLines = []

        for line in lines:
            if len(line.strip()) > 0:
                processedLines.append(line)
        return processedLines

    def createMolecules(processedStructureLines):
        """
        Converts the 'processedStructureLines'  lines into 'Atoms' and creates 'Molecules'

        """

        systemName = processedStructureLines.pop(0).strip()
        numAtoms = processedStructureLines.pop(0)
        rawBoxVector = processedStructureLines.pop(-1).split()

        tempList = list()
        atomList = list()
        moleculeList = list()

        # Create Atoms
        for line in processedStructureLines:
            resNum=int(line[0:5])
            resName=line[6:10].strip()
            atomName=line[11:15].strip()
            atomNum=int(line[16:20])
            x=float(line[20:28])
            y=float(line[28:36])
            z=float(line[36:44])
            try:
                vx=float(line[44:52])
                vy=float(line[52:60])
                vz=float(line[60:68])
            except: # if there aren't any velocities yet
                vx=0.0
                vy=0.0
                vz=0.0
            atom=Atom(atomNum, 
                    resNum, 
                    resName, 
                    atomName, 
                    x * units.nanometers, 
                    y * units.nanometers,
                    z * units.nanometers, 
                    vx * units.nanometers * units.picoseconds**(-1), 
                    vy * units.nanometers * units.picoseconds**(-1), 
                    vz * units.nanometers * units.picoseconds**(-1))
            atomList.append(atom)
            tempList.append(atom)

        # Extract protein
        residuetypes = os.path.join(os.environ['GMXLIB'], 'residuetypes.dat')
        if os.path.exists(residuetypes):
            proteinList = set()
            fin = open(residuetypes, 'r')
            for line in fin.readlines():
                line = line.split()
                res = line[0].strip()
                if line[1].strip() == 'Protein':
                    proteinList.add(res)
        else:
            # sys.stderr.write()
            print "ERROR: 'residuetypes.dat' not found in "+os.getenv('GMXLIB')
        currentMolecule = Molecule('Protein')
       
        i = 0
        while i < len(tempList):           
            if tempList[i].resName.strip() in proteinList:
                if tempList[i].atomNum == 129:
                    pdb.set_trace()
                x = tempList.pop(i)
                currentMolecule.atoms.add(x)
                print currentMolecule.atoms
            else:
                i=i+1
        
        if currentMolecule.atoms:
            moleculeList.append(currentMolecule)
        """
        for atom in tempList:
            resName = atom.resName
            if resName in proteinList:
                currentMolecule.atoms.add(atom)
        if currentMolecule.atoms:
            moleculeList.append(currentMolecule)
        
        for atom in currentMolecule.atoms:
            if atom in tempList:
                tempList.remove(atom)
        """
        # Extract non-protein molecules
        if len(tempList) > 0:
            oldResNum = tempList[0].resNum
            currentMolecule = Molecule(tempList[0].resName)
        for atom in tempList:
            resNum = atom.resNum
            resName = atom.resName
            if resNum != oldResNum:
                moleculeList.append(currentMolecule)
                currentMolecule = Molecule()
                currentMolecule.atoms.add(atom)
                currentMolecule.name = resName
                oldResNum = resNum
            else:
                    currentMolecule.atoms.add(atom)

        if currentMolecule not in moleculeList:
            moleculeList.append(currentMolecule)

        # Process box vector
        if len(rawBoxVector) == 3:
            v1x = float(rawBoxVector[0]) * units.nanometers
            v2x = 0.0 * units.nanometers
            v3x = 0.0 * units.nanometers
            v1y = 0.0 * units.nanometers
            v2y = float(rawBoxVector[1]) * units.nanometers
            v3y = 0.0 * units.nanometers
            v1z = 0.0 * units.nanometers
            v2z = 0.0 * units.nanometers
            v3z = float(rawBoxVector[2]) * units.nanometers
 
        if len(rawBoxVector) == 9:
            v1x = float(rawBoxVector[0]) * units.nanometers
            v2x = float(rawBoxVector[5]) * units.nanometers
            v3x = float(rawBoxVector[7]) * units.nanometers
            v1y = float(rawBoxVector[3]) * units.nanometers
            v2y = float(rawBoxVector[1]) * units.nanometers
            v3y = float(rawBoxVector[8]) * units.nanometers
            v1z = float(rawBoxVector[4]) * units.nanometers
            v2z = float(rawBoxVector[6]) * units.nanometers
            v3z = float(rawBoxVector[2]) * units.nanometers

        return systemName, v1x, v2x, v3x, v1y, v2y, v3y, v1z, v2z, v3z, moleculeList, atomList

    # read in the lines of a .gro file
    lines = readLines(filename)
    # preprocess the lines
    processedStructureLines = preprocess(lines)
    
    # create 'Atom' and 'Molecule' objects and populate the 'molecules' dictionary in System
    return createMolecules(processedStructureLines)

def writeStructure(sys, filename):

    lines = list()

    lines.append(sys.name+'\n')
    lines.append('%d\n' %(len(sys.atoms)))

    for atom in sys.atoms:
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
   
    lines.append('%11.7f%11.7f%11.7f%11.7f%11.7f%11.7f%11.7f%11.7f%11.7f\n'
          %(sys.v1x._value,
            sys.v2y._value,
            sys.v3z._value,
            sys.v1y._value,
            sys.v1z._value,
            sys.v2x._value,
            sys.v2z._value,
            sys.v3x._value,
            sys.v3y._value))


    fout = open(filename, 'w')
    for line in lines:
        fout.write(line)
    fout.close()
    
    
