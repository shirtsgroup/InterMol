from cctools.Decorators import *
from cctools.Types.AbstractPairType import *

class LJ1PairCR1Type(AbstractPairType):
    @accepts_compatible_units(None, None, None, units.kilojoules_per_mole * units.nanometers**(6), units.kilojoules_per_mole * units.nanometers**(12))
    def __init__(self, atom1, atom2, type, V, W):
        AbstractPairType.__init__(self, atom1, atom2, type)
        self.V = V
        self.W = W

