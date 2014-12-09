class AbstractPairType(object):
    __slots__ = ['bondingtype1', 'bondingtype2', 'scaleLJ', 'scaleQQ', 'long']

    def __init__(self, bondingtype1, bondingtype2, scaleLJ = None, scaleQQ = None, long = False):
        self.bondingtype1 = bondingtype1
        self.bondingtype2 = bondingtype2
        self.scaleLJ = scaleLJ
        self.scaleQQ = scaleQQ
        self.long = long

    def __repr__(self):
        return  "{0}({1})".format(self.__class__.__name__, ", ".join(["{0}={1}".format(x, getattr(self, x, "Undefined")) for x in dir(self) if not(x.endswith("__") or hasattr(getattr(self, x, "Undefined"), '__call__'))]))
