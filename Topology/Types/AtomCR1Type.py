from cctools.Decorators import *
from AbstractAtomType import AbstractAtomType

class AtomCR1Type(AbstractAtomType):
    @accepts_compatible_units(None, None, units.amu, units.elementary_charge, None, units.kilojoules_per_mole * units.nanometers**(6), units.kilojoules_per_mole * units.nanometers**(12))
    def __init__(self, atomtype, Z, m, q, ptype, V, W):
        AbstractAtomType.__init__(self, atomtype, Z, m, q, ptype)
        self.V = V
        self.W = W

