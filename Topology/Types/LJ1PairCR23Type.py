from cctools.Decorators import *
from cctools.Types.AbstractPairType import *

class LJ1PairCR23Type(AbstractPairType):
    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole)
    def __init__(self, atom1, atom2, type, V, W):
        AbstractPairType.__init__(self, atom1, atom2, type)
        self.V = V
        self.W = W

