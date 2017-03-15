from collections import OrderedDict
import logging

import numpy as np
import parmed.unit as units

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

        self._bondgraph = None

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
        for mol_type in self.molecule_types.values():
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
        # handle two different 'unit' module behaviors   
        if hasattr(v,'unit'):    
            self._box_vector = np.array(v)*v[0,0].unit
        else:
            self._box_vector = np.array(v)

    @property
    def connected_pairs(self):
        for mol_type in self.molecule_types.values():
            for molecule in mol_type.molecules:
                for bond in mol_type.bond_forces:
                    atom1 = molecule.atoms[bond.atom1 - 1]
                    atom2 = molecule.atoms[bond.atom2 - 1]
                    yield atom1, atom2

    @property
    def bondgraph(self):
        """Create a NetworkX graph from the atoms and bonds in this system.

        Returns
        -------
        g : nx.Graph
            A graph whose nodes are the Atoms in this topology, and
            whose edges are the bonds

        See Also
        --------
        connected_pairs

        Notes
        -----
        This method requires the NetworkX python package.
        """
        if self._bondgraph is not None:
            return self._bondgraph
        else:
            import networkx as nx
            self._bondgraph = nx.Graph()
            self._bondgraph.add_edges_from(self.connected_pairs)
            return self._bondgraph

    # def gen_pairs(self, n_excl=4):
    #
    #     # loop over moleculetypes
    #     #   loop over dihedral forces
    #     #   loop over pairs and add 1-nexcl pairs
    #     if n_excl == 4:
    #         for dihedral in self.dihedrals:
    #             self.pair.add((dihedral.atom1, dihedral.atom4))
    #     else:
    #         raise ValueError('Unsupported number of pair exclusions.')

    def __repr__(self):
        return "System '{}' ".format(self.name)

    def __str__(self):
        return "System{} '{}'".format(id(self), self.name)


