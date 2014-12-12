class Molecule(object):
    """An abstract molecule object. """

    def __init__(self, name=None):
        """Initialize the molecule

        Args:
            name (str): name of the molecule
        """
        if not name:
            name = "MOL"
        self.name = name
        self._atoms = list()

    def add_atom(self, atom):
        """Add an atom

        Args:
            atom (Atom): the atom to add into the molecule
        """
        self._atoms.append(atom)

    @property
    def atoms(self):
        """Return an orderedset of atoms. """
        return self._atoms

    def __repr__(self):
        return "Molecule '{}' with {} atoms".format(self.name, len(self.atoms))

    def __str__(self):
        return "Molecule{} '{}' with {} atoms".format(
            id(self), self.name, len(self.atoms))
