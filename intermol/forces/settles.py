from intermol.decorators import *


class Settles(object):

    @accepts_compatible_units(None, units.nanometers, units.nanometers)
    def __init__(self, atom1, dOH, dHH):
        """
        """
        if atom1:
            self.atom1 = atom1
        self.dOH = dOH
        self.dHH = dHH

    def __repr__(self):
        return "{0}({1})".format(self.__class__.__name__,
                        ", ".join(["{0}={1}".format(x, getattr(self, x, "Undefined"))
                        for x in dir(self) if not(x.endswith("__") or
                                                  hasattr(getattr(self, x, "Undefined"), '__call__'))]))
