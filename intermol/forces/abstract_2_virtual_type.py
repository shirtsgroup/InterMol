class Abstract2VirtualType(object):
    __slots__ = ['bondingtype1', 'bondingtype2', 'bondingtype2', 'placeholder']

    def __init__(self, bondingtype1, bondingtype2, bondingtype3, placeholder=False):
        self.bondingtype1 = bondingtype1
        self.bondingtype2 = bondingtype2
        self.bondingtype3 = bondingtype3

    def __repr__(self):
        attributes = ["{0}={1}".format(x, getattr(self, x, "Undefined")) for x
                in dir(self)
                if not( x.endswith("__") or hasattr(getattr(self, x, "Undefined"), '__call__'))]
        return  "{0}({1})".format(self.__class__.__name__, ", ".join(attributes))
