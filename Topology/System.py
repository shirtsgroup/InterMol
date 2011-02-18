from Decorators import *
from MoleculeType import MoleculeType
from OrderedSet import OrderedSet
from OrderedDict import OrderedDict

class System(object):
    def __init__(self, name = None):
        """
        Create a new System object.
        """
        if name:
            self.name = name
        else:
            self.name = "Untitled"
        
        self.v1x = 0.0
        self.v2x = 0.0
        self.v3x = 0.0
        self.v1y = 0.0
        self.v2y = 0.0
        self.v3y = 0.0
        self.v1z = 0.0
        self.v2z = 0.0
        self.v3z = 0.0

        self.nbFunc = 0
        self.combinationRule = 0
        self.genpairs = 'yes'
        self.ljCorrection = 0
        self.coulombCorrection = 0  
        self.molecules = OrderedDict()
        self.atoms = list()
        self.forces = OrderedSet()


    def addMolecule(self, molecule, name):
        # if key is in the dictionary, return its value. If not, insert key with a value of default and return default.        
        self.molecules.setdefault(molecule.name,MoleculeType(name)).addMolecule(molecule)
    
    def delMoleculeType(self, name, molecule):
        self.molecules[name].remove(molecule)

    def getBoxVector(self):
        return self.boxVector

    @accepts_compatible_units(units.nanometers, units.nanometers, units.nanometers, units.nanometers, units.nanometers, units.nanometers,units.nanometers, units.nanometers, units.nanometers)

    def setBoxVector(self, v1x, v2x, v3x, v1y, v2y, v3y, v1z, v2z, v3z):
        """
        Sets the boxvector for the system. Assumes the box vector is in the correct form.
        [x, y, z]
        or
        [[v1x,v2x,v3x],[v1y,v2y,v3y],[v1z,v2z,v3z]]
        """
        self.v1x = v1x
        self.v2x = v2x
        self.v3x = v3x
        self.v1y = v1y
        self.v2y = v2y
        self.v3y = v3y
        self.v1z = v1z
        self.v2z = v2z
        self.v3z = v3z

    def getNBFunc(self):
        return self.nbFunc

    def setNBFunc(self, nbFunc):
        self.nbFunc = nbFunc
    
    def getCombinationRule(self):
        return self.combinationRule
    
    def setCombinationRule(self, combinationRule):
        self.combinationRule = combinationRule

    def getLJCorrection(self):
        return self.ljCorrection
    
    def setLJCorrection(self, ljCorrection):
        self.ljCorrection = ljCorrection 

    def __repr__(self):
        return self.name
