class AbstractNonbondedType(object):
    __slots__ = ['atom1', 'atom2', 'type']

    def __init__(self, atom1, atom2, type):
        self.atom1 = atom1
        self.atom2 = atom2
        self.type = type

    def __repr__(self):
        return  "{0}({1})".format(self.__class__.__name__, ", ".join(["{0}={1}".format(x, getattr(self, x, "Undefined")) for x in dir(self) if not(x.endswith("__") or hasattr(getattr(self, x, "Undefined"), '__call__'))]))
