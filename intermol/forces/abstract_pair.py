class AbstractPair(object):
    __slots__ = ['atom1', 'atom2', 'type']

    def __init__(self, atom1, atom2, type = None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.type = type

    def __eq__(self, object):
        if (self.atom1 == (object.atom1 or object.atom2)) and (self.atom2 == (object.atom2 or object.atom1)):
            return True
        else:
            return False

    def __hash__(self):
        return hash(tuple([self.atom1, self.atom2]))

