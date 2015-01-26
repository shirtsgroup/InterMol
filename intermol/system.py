from collections import OrderedDict
import logging

import numpy as np
import simtk.unit as units

logger = logging.getLogger('InterMolLog')  # TODO: do we need the logger here?


class System(object):
    """  """

    def __init__(self, name=None):
        """Initialize a new System object.

        Args:
            name (str): The name of the system
        """
        if name:
            self.name = name
        else:
            self.name = "Untitled"

        self.nonbonded_function = 0
        self.combination_rule = 0
        self.genpairs = 'yes'
        self.lj_correction = 0
        self.coulomb_correction = 0

        self._box_vector = np.zeros([3, 3]) * units.nanometers

        self._n_atoms = None
        self._molecule_types = OrderedDict()
        self._atomtypes = dict()
        self._nonbonded_types = dict()

    def add_molecule(self, molecule):
        """Add a molecule into the System. """
        self._molecule_types[molecule.name].add_molecule(molecule)

    def add_molecule_type(self, molecule_type):
        """Add a molecule_type into the System. """
        self._molecule_types[molecule_type.name] = molecule_type

    def add_atomtype(self, atomtype):
        """ """
        self._atomtypes[atomtype.atomtype] = atomtype

    @property
    def atomtypes(self):
        return self._atomtypes

    @property
    def nonbonded_types(self):
        return self._nonbonded_types

    @property
    def molecule_types(self):
        return self._molecule_types

    @property
    def atoms(self):
        for mol_type in self.molecule_types.itervalues():
            for mol in mol_type.molecules:
                for atom in mol.atoms:
                    yield atom

    @property
    def n_atoms(self):
        if not self._n_atoms:
            self._n_atoms = len(list(self.atoms))
        return self._n_atoms

    @n_atoms.setter
    def n_atoms(self, n):
        self._n_atoms = n

    @property
    def box_vector(self):
        """Return the box vector. """
        return self._box_vector

    @box_vector.setter
    def box_vector(self, v):
        """Sets the box vector for the system.

        Assumes the box vector is in the correct form:
            [[v1x,v2x,v3x],[v1y,v2y,v3y],[v1z,v2z,v3z]]
        """
        if v.shape != (3, 3):
            e = ValueError("Box vector with incorrect format: {0}".format(v))
            logger.exception(e)
        self._box_vector = np.array(v)

    def __repr__(self):
        return "System '{}' ".format(self.name)

    def __str__(self):
        return "System{} '{}'".format(id(self), self.name)


