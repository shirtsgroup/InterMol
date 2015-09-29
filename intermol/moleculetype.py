from intermol.orderedset import OrderedSet


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
        self.rigidwaters = set()
        self.nrexcl = None

    def add_molecule(self, molecule):
        """Add a molecule into the moleculetype. """
        self.molecules.add(molecule)

    # These three functions could probably be made faster through some sort of list comprehension
    # rather than iterating explicitly over the lists

    def _match_two_atoms(self, newforce, oldforces):
        """ find and return any force with the same three atoms. For now, respect ordering in matches. """
        newatoms = [newforce.atom1, newforce.atom2]
        for force in oldforces:
            oldatoms = [force.atom1, force.atom2]
            if newatoms == oldatoms:
                return force
        return False

    def _match_three_atoms(self, newforce, oldforces):
        """ find and return any force with the same three atoms. For now, respect ordering in matches. """
        newatoms = [newforce.atom1, newforce.atom2, newforce.atom3]
        for force in oldforces:
            oldatoms = [force.atom1, force.atom2, force.atom3]
            if newatoms == oldatoms:
                return force
        return False

    def _match_four_atoms(self, newforce, oldforces):
        """ find and return any force with the same three atoms. For now, respect ordering in matches. """
        newatoms = [newforce.atom1, newforce.atom2, newforce.atom3, newforce.atom4]
        for force in oldforces:
            oldatoms = [force.atom1, force.atom2, force.atom3, force.atom4]
            if newatoms == oldatoms:
                return force
        return False

    def match_bonds(self, bond):
        return self._match_two_atoms(bond, self.bond_forces)

    def match_pairs(self, pair):
        return self._match_two_atoms(pair, self.pair_forces)

    def match_angles(self, angle):
        return self._match_three_atoms(angle, self.angle_forces)

    def match_dihedrals(self, dihedral):
        return self._match_four_atoms(dihedral, self.dihedral_forces)

    def __repr__(self):
        return "MoleculeType '{}' with {} molecules".format(
            self.name, len(self.molecules))

    def __str__(self):
        return "MoleculeType{} '{}' with {} molecules".format(
            id(self), self.name, len(self.molecules))
