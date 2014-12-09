class AbstractAngleType(object):
    __slots__ = ['bondingtype1', 'bondingtype2', 'bondingtype3', 'c']

    def __init__(self, bondingtype1, bondingtype2, bondingtype3, c = False):
        """An abstract representation of a generic angle type."""
        self.bondingtype1 = bondingtype1
        self.bondingtype2 = bondingtype2
        self.bondingtype3 = bondingtype3
        self.c     = c  # constrained angle, for desmond
        

    def __repr__(self):
        return  "{0}({1})".format(self.__class__.__name__, ", ".join(["{0}={1}".format(x, getattr(self, x, "Undefined")) for x in dir(self) if not(x.endswith("__") or hasattr(getattr(self, x, "Undefined"), '__call__'))]))
