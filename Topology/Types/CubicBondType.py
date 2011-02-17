from cctools.Decorators import *
from cctools.Types.AbstractBondType import *

class CubicBondType(AbstractBondType):
    @accepts_compatible_units(None, 
            None, 
            None, 
            units.nanometers, 
            units.kilojoules_per_mole * units.nanometers**(-2), 
            units.kilojoules_per_mole * units.nanometers**(-3))
    def __init__(self, atom1, atom2, type, length, C2, C3):
        AbstractBondType.__init__(self, atom1, atom2, type)
        self.length = length
        self.C2 = C2
        self.C3 = C3

