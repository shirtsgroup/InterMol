"""
.. module:: System
    :platform: UNIX
"""
from ctools.Converter import *
from ctools.MoleculeType import MoleculeType
from ctools.OrderedSet import OrderedSet
from ctools.OrderedDict import OrderedDict
from ctools.HashMap import HashMap

class System(object):
    _sys = None
    def __init__(self, name = None):
        """Initialize a new System object. This must be run before the system can be used.

        Args:
            name (str): The name of the system

        >>> print __init__(name='sysname')
        """
        if name:
            self.name = name
        else:
            self.name = "Untitled"
        
        self._v1x = 0.0
        self._v2x = 0.0
        self._v3x = 0.0
        self._v1y = 0.0
        self._v2y = 0.0
        self._v3y = 0.0
        self._v1z = 0.0
        self._v2z = 0.0
        self._v3z = 0.0

        self._nbFunc = 0
        self._combinationRule = 0
        self._genpairs = True
        self._ljCorrection = 0
        self._coulombCorrection = 0  
        self._molecules = OrderedDict()
        self._atomtypes = HashMap()
        self._forces = OrderedSet()


    def addMolecule(self, molecule):
        """Append a molecule into the System.

        Args:
            molecule (Molecule): The molecule object to be appended
        """
        # if key is in the dictionary, return its value. If not, insert key with a value of default and return default.        
        self._molecules.setdefault(molecule.name,MoleculeType(molecule.name)).addMolecule(molecule)
    
    def removeMoleculeType(self, molecule):
        """Remove a molecule from the System.
    
        Args:
           molecule (Molecule): The molecule object to be removed 
        """
        self._molecules[molecule.name].remove(molecule)

    def getBoxVector(self):
        """Get the box vector coordinates
        """
        return self._boxVector


    def setBoxVector(self, v1x, v2x, v3x, v1y, v2y, v3y, v1z, v2z, v3z):
        """Sets the boxvector for the system. Assumes the box vector is in the correct form. [[v1x,v2x,v3x],[v1y,v2y,v3y],[v1z,v2z,v3z]]
        """
        unit = units.nanometers
        self._v1x = v1x
        self._v2x = v2x
        self._v3x = v3x
        self._v1y = v1y
        self._v2y = v2y
        self._v3y = v3y
        self._v1z = v1z
        self._v2z = v2z
        self._v3z = v3z

    def __str__(self):
        """String representation of a System object
        """
        return "System: " + self.name
        
    def __repr__(self):
        """String representation of a System object
        """
        return "System: " + self.name
