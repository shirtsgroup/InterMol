from cctools.Decorators import *
from AbstractAtomType import AbstractAtomType


class AtomCR23Type(AbstractAtomType):
    @accepts_compatible_units(None, None, units.amu, units.elementary_charge, None, units.nanometers, units.kilojoules_per_mole)

    def __init__(self, atomtype, Z, m, q, ptype, V, W):
        AbstractAtomType.__init__(self, atomtype, Z, m, q, ptype)
        self.V = V
        self.W = W


