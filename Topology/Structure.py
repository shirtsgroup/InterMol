tlc2olc = {
"ALA" : "A" ,
"CYS" : "C" ,
"CYN" : "C" ,
"CYX" : "C" ,
"ASP" : "D" ,
"GLU" : "E" ,
"PHE" : "F" ,
"CPHE" : "F" ,
"GLY" : "G" ,
"HIS" : "H" ,
"HIE" : "H" ,
"HIP" : "H" ,
"ILE" : "I" ,
"NLE" : "J" ,
"LYS" : "K" ,
"LYP" : "K" ,
"LEU" : "L" ,
"NLEU" : "L" ,
"MET" : "M" ,
"ASN" : "N" ,
"PRO" : "P" ,
"GLN" : "Q" ,
"ARG" : "R" ,
"SER" : "S" ,
"THR" : "T" ,
"VAL" : "V" ,
"TRP" : "W" ,
"TYR" : "Y"
}

import pdb

class Structure(object):
    """Parent class for structure objects"""
    def __init__(self):
        self.molecules = []

    def addMolecule(self):
        pass

    def convertFromGromacsStructure(self, gromacsStructure):
        for molecule in gromacsStructure.molecules:
            newMolecule = []
            for atom in molecule.atomList:
                positions = [atom.x, atom.y, atom.z]
                velocities = [atom.vx, atom.vy, atom.vz]
                newAtom = Atom(positions, velocities)
                newMolecule.append(newAtom)
            self.molecules.append(newMolecule)
        return

    def convertToGromacsStructure(self):

        return



class Atom(object):

    def __init__(self, positions=[0.0, 0.0, 0.0], velocities=[0.0, 0.0, 0.0], forces=[0.0, 0.0, 0.0]):
        self.positions = positions
        self.velocities = velocities
        self.forces = forces



class GromacsStructure(Structure):
#===================WARNING==========================
# ONLY WORKS FOR .GRO FILES WITH ZERO OR ONE PROTEIN(S)
#====================================================

    def __init__(self, structure=None, grofile=None, name=None):

        self.title = ''
        self.nSystem = ''
        self.molecules = list()
        self.boxvector = list()

        if (grofile):
            self.GromacsStructureFileObject = GromacsStructureFile(grofile=grofile)
        else:
            self.GromacsStructureFileObject = None

    def setName(self, name):
        self.name = name
        return

    def readStructureFile(self, grofile):
        self.GromacsStructureFileObject = GromacsStructureFile(grofile=grofile)
        self.title = self.GromacsStructureFileObject.title
        self.nSystem = self.GromacsStructureFileObject.nSystem
        self.molecules = self.processMolecules()
        self.boxvector = self.processBoxVector()
        return

    def processMolecules(self):

        tempList = list()
        moleculeList = list()

        #create atoms from "preprocessedMolecules"
        for line in self.GromacsStructureFileObject.preprocessedMolecules:
            resnum=int(line[0:5])
            resname=line[5:9]
            atomname=line[9:15]
            atomnum=int(line[15:20])
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
            atom=GromacsAtom(resnum,resname,atomname,atomnum,x,y,z,vx,vy,vz)
            tempList.append(atom)

        #extract protein
        currentMolecule = GromacsMolecule()
        for atom in tempList:
            resname = atom.getResname()
            if resname.strip() in tlc2olc:
                currentMolecule.name = 'Protein'
                currentMolecule.atomList.append(atom)
        if currentMolecule.atomList: moleculeList.append(currentMolecule)

        for atom in currentMolecule.atomList:
            if atom in tempList:
                tempList.remove(atom)

        #extract other molecules
        if len(tempList) > 0:
            oldResnum = tempList[0].getResnum()
            currentMolecule = GromacsMolecule()
            currentMolecule.name = tempList[0].getResname()
        for atom in tempList:
            resnum = atom.getResnum()
            resname = atom.getResname()
            if resnum != oldResnum:
                moleculeList.append(currentMolecule)
                currentMolecule = GromacsMolecule()
                currentMolecule.atomList.append(atom)
                currentMolecule.name = resname
                oldResnum = resnum
            else:
                currentMolecule.atomList.append(atom)
        if currentMolecule not in moleculeList:
            moleculeList.append(currentMolecule)

        return moleculeList

    def processBoxVector(self):
        fields = self.GromacsStructureFileObject.boxvector.split()
        nBoxVector = len(fields)

        if nBoxVector == 3:
            # diagonal
            v1x = float(fields[0])
            v2y = float(fields[1])
            v3z = float(fields[2])
            boxvector = [[v1x,0,0],[0,v2y,0],[0,0,v3z]]
        if nBoxVector == 9:
            # not diagonal
            v1x = float(fields[0])
            v2y = float(fields[1])
            v3z = float(fields[2])
            v1y = float(fields[3])
            v1z = float(fields[4])
            v2x = float(fields[5])
            v2z = float(fields[6])
            v3x = float(fields[7])
            v3y = float(fields[8])
            boxvector = [[v1x,v2x,v3x],[v1y,v2y,v3y],[v1z,v2z,v3z]]
        return boxvector

    def writeStructureFile(self, filename):
        self.GromacsStructureFileObject.write(filename)
        return

class GromacsAtom(object):
    """A type which contains all the known information about one atom in a gromacs structure,
    topology, or trajectory. (For now just works with the gro structure file.)

    attributes:
    -----------
    - resnum
    - resname
    - atomname
    - atomnum
    - x, y, and z
    - vx, vy, and vz

    methods:
    --------
    - __copy__()        - makes a copy of the atom (TEST!)
    - position()        - returns the position as a 3-d tuple (make this work on the velocity also)
    - printgroatom()- prints the atom as it's formatted in the gro file (used by the GromacsStructure type)
    """

    def __init__(self, resnum=0, resname="", atomname="", atomnum=0, x=0.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0):
        self.resnum=resnum
        self.resname=resname
        self.atomname=atomname
        self.atomnum=atomnum
        self.x=x
        self.y=y
        self.z=z
        self.vx=vx
        self.vy=vy
        self.vz=vz

    def getResnum(self):
        return self.resnum

    def getResname(self):
        return self.resname

    def __copy__(self):
        return GromacsAtom(self.resnum, self.resname, self.atomname, self.atomnum, self.x, self.y, self.z)

    def position(self):
        # is there a way to make this return the velocity vector as well?
        return (self.x, self.y, self.z)

    def printGroAtom(self):
        """prints the atom as it would be formatted in the .gro file"""
        resnum="%5d" % self.resnum
        resname="%-4s" % self.resname # should be left justified
        atomname="%6s" % self.atomname
        atomnum="%5d" % self.atomnum
        pos="% 8.3f% 8.3f% 8.3f" % (self.x, self.y, self.z)
        vel="% 8.4f% 8.4f% 8.4f" % (self.vx, self.vy, self.vz)
        return resnum+resname+atomname+atomnum+pos+vel

class GromacsMolecule(object):

    def __init__(self, name = None):
        self.name = name
        if self.name == None: 
            self.name = "Untitled"

        self.atomList = list()

    def getLength(self):
        return len(self.atomList)

    def __copy__(self): return GromacsMolecule(self)

    def printGroMolecule(self):
        outtxt = ''
        for atom in self.atomList:
            outtxt += atom.printGroAtom() + '\n'
        return outtxt

class GromacsStructureFile(object):

    def __init__(self, grofile = None):

        self.lines = []
        self.processedLines = []

        if not grofile == None:
            self.read(grofile)

        # Here is where title, nSystem, molecules and boxvector will be defined
        """
        self.title (string)
        self.nSystem (int)
        self.molecules (list of lists of string/int/float)
            molecules will be composed of various parts
            residue_number residue_name atom_name atom_number x y z vx vy vz
        self.boxvector (list/matrix of float)
            similar matrix of either 3 or 6 values
        """

    def read(self, filename):

        # read in the *.gro file lines
        self.lines = self.readLines(filename)

        # Perform pre-processing
        self.processedLines = self.preprocess()

        self.title = self.processedLines.pop(0)
        self.nSystem = self.processedLines.pop(0)
        self.boxvector = self.processedLines.pop(-1)
        self.preprocessedMolecules = self.processedLines
        return


    def readLines(self, filename):

        fin = open(filename, 'r')
        lines = fin.readlines()
        fin.close()
        return lines

    def preprocess(self):

        processedLines = []

        for line in self.lines:
            if len(line.strip()) > 0:
                processedLines.append(line)
        return processedLines

    def write(self, filename):
        fout = open(filename, 'w')
        fout.write(repr(self))
        fout.close()
        return

    def __repr__(self):
        outtxt = self.title
        outtxt += self.nSystem
        for mol in self.preprocessedMolecules:
            outtxt += mol
        outtxt += self.boxvector
        return outtxt


class PDBStructure(Structure):
    pass

class AmberStructure(Structure):
    pass

class DesmondStructure(Structure):
    pass

class SchroedingerStructure(Structure):
    pass

