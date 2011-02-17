from cctools.Decorators import *
from cctools.Types.AbstractBondType import *

class MorseBondType(AbstractBondType):
    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole,  units.nanometers**(-1))
    def __init__(self, atom1, atom2, type, length, D, beta):
        AbstractBondType.__init__(self, atom1, atom2, type)
        self.length = length
        self.D = D
        self.beta = beta
