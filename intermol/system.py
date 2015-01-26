"""
.. module:: System
    :platform: UNIX
"""
from collections import OrderedDict
import numpy as np

from intermol.moleculetype import MoleculeType
from intermol.hashmap import HashMap


class System(object):
    _sys = None

    def __init__(self, name=None):
        """Initialize a new System object.

        This must be run before the system can be used.

        Args:
            name (str): The name of the system

        >>> print __init__(name='sysname')
        """
        if name:
            self._name = name
        else:
            self._name = "Untitled"

        self.nonbonded_function = 0
        self.combination_rule = 0
        self.genpairs = 'yes'
        self.lj_correction = 0
        self.coulomb_correction = 0
        self._molecules = OrderedDict()
        self._atomtypes = HashMap()
        self._nonbonded = HashMap()
        self._box_vector = np.zeros([3, 3])

        self._components = list()


    def add_molecule(self, molecule):
        """Append a molecule into the System.

        Args:
            molecule (Molecule): The molecule object to be appended
        """
        # If key is in the dictionary, return its value.
        # If not, insert key with a value of default and return default.
        self._molecules.setdefault(molecule.name,
                MoleculeType(molecule.name)).add_molecule(molecule)

    def remove_molecule_type(self, molecule):
        """Remove a molecule from the System.

        Args:
           molecule (Molecule): The molecule object to be removed
        """
        self._molecules[molecule.name].remove(molecule)

    @property
    def box_vector(self):
        """Get the box vector coordinates
        """
        return self._box_vector

    @box_vector.setter
    def box_vector(self, v):
        """Sets the boxvector for the system.

        Assumes the box vector is in the correct form:
            [[v1x,v2x,v3x],[v1y,v2y,v3y],[v1z,v2z,v3z]]
        """
        self._box_vector = v

    def __str__(self):
        """String representation of a System object
        """
        return "System: " + self._name

    def __repr__(self):
        """String representation of a System object
        """
        return "System: " + self._name
