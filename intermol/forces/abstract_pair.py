class AbstractPair(object):
    __slots__ = ['atom1', 'atom2', 'pairtype']

    def __init__(self, atom1, atom2, pairtype = None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.pairtype = pairtype

    def __eq__(self, object):
        if ((self.atom1 == object.atom1 and self.atom2 == object.atom2) or (self.atom1 == object.atom2 and self.atom2 == object.atom1)):
            return True
        else:
            return False

    def __hash__(self):
        return hash(tuple(sorted([self.atom1, self.atom2])))

