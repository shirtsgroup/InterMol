from orderedset import OrderedSet


class MoleculeType(object):
    """An abstract container for molecules of one type. """
    def __init__(self, name=None):
        """Initialize the MoleculeType container.

        Args:
            name (str): the name of the moleculetype to add
        """
        if not name:
            name = 'MOL'
        self.name = name
        self.molecules = OrderedSet()

        self.bond_forces = set()
        self.pair_forces = set()
        self.angle_forces = set()
        self.dihedral_forces = set()
        self.virtual_forces = set()
        self.torsiontorsion_forces = set()
        self.constraints = set()
        self.exclusions = set()
        self.settles = None
        self.nrexcl = None

    def add_molecule(self, molecule):
        """Add a molecule into the moleculetype. """
        self.molecules.add(molecule)

    def __repr__(self):
        return "MoleculeType '{}' with {} molecules".format(
            self.name, len(self.molecules))

    def __str__(self):
        return "MoleculeType{} '{}' with {} molecules".format(
            id(self), self.name, len(self.molecules))
