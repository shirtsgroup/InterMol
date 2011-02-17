from cctools.Decorators import *
from cctools.Types.AbstractPairType import *

class LJNBPairCR23Type(AbstractPairType):
    @accepts_compatible_units(None, None, None, units.elementary_charge, units.elementary_charge, units.nanometers, units.kilojoules_per_mole)
    def __init__(self, atom1, atom2, type, qi, qj, V, W):
        AbstractPairType.__init__(self, atom1, atom2, type)
        self.qi = qi
        self.qj = qj
        self.V = V
        self.W = W
