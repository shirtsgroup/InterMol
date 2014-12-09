class AbstractBondType(object):
    __slots__ = ['bondingtype1', 'bondingtype2', 'order', 'c']

    def __init__(self, bondingtype1, bondingtype2, order=1, c=False):
        self.bondingtype1 = bondingtype1
        self.bondingtype2 = bondingtype2
        self.order = order  # what is the bond order, Desmond only
        self.c = c  # is the bond constrained or not, Desmond only


    def __repr__(self):
        attributes = ["{0}={1}".format(x, getattr(self, x, "Undefined")) for x
                      in dir(self)
                      if not (x.endswith("__") or hasattr(getattr(self, x, "Undefined"), '__call__'))]
        return "{0}({1})".format(self.__class__.__name__, ", ".join(attributes))
