"""
.. module:: System
    :platform: UNIX
"""

import numpy as np
from intermol.OrderedDict import OrderedDict
from intermol.Converter import *
from intermol.MoleculeType import MoleculeType
from intermol.OrderedSet import OrderedSet
from intermol.HashMap import HashMap

class System(object):
    _sys = None
    def __init__(self, name = None):
        """Initialize a new System object. This must be run before the system can be used.

        Args:
            name (str): The name of the system

        >>> print __init__(name='sysname')
        """
        if name:
            self._name = name
        else:
            self._name = "Untitled"

        self._nbFunc = 0
        self._combinationRule = 0
        self._genpairs = True
        self._ljCorrection = 0
        self._coulombCorrection = 0
        self._molecules = OrderedDict()
        self._atomtypes = HashMap()
        self._forces = OrderedSet()
        self._boxVector = np.zeros([3,3],float)


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


    def setBoxVector(self, v):
        """Sets the boxvector for the system. Assumes the box vector is in the correct form. [[v1x,v2x,v3x],[v1y,v2y,v3y],[v1z,v2z,v3z]]
        """
        self._boxVector = v

    def __str__(self):
        """String representation of a System object
        """
        return "System: " + self.name

    def __repr__(self):
        """String representation of a System object
        """
        return "System: " + self.name
