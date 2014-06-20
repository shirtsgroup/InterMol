class AbstractBondType(object):
    __slots__ = ['atom1', 'atom2']
    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2

    def __eq__(self, bond_type):
        return ((self.atom1 == bond_type.atom1 and
                 self.atom2 == bond_type.atom2)
                or
                (self.atom2 == bond_type.atom2 and
                 self.atom3 == bond_type.atom1))

    def __hash__(self):
        return hash(tuple([self.atom1, self.atom2]))

