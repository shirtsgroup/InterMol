class Abstract4VirtualType(object):
    __slots__ = ['bondingtype1', 'bondingtype2', 'bondingtype3', 'bondingtype4', 'bondingtype4', 'placeholder']

    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, bondingtype5, placeholder=False):
        self.bondingtype1 = bondingtype1
        self.bondingtype2 = bondingtype2
        self.bondingtype3 = bondingtype3
        self.bondingtype4 = bondingtype4
        self.bondingtype5 = bondingtype5

    def __repr__(self):
        return  "{0}({1})".format(self.__class__.__name__, ", ".join(["{0}={1}".format(x, getattr(self, x, "Undefined")) for x in dir(self) if not(x.endswith("__") or hasattr(getattr(self, x, "Undefined"), '__call__'))]))
