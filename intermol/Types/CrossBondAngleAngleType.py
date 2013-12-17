import sys
sys.path.append('..')

from intermol.Decorators import *
from AbstractAngleType import *

class CrossBondAngleAngleType(AbstractAngleType):
    @accepts_compatible_units(None,
            None,
            None,
            None,
            units.nanometers,
            units.nanometers,
            units.nanometers,
            units.kilojoules_per_mole * units.nanometers**(-2))
    def __init__(self, atom1, atom2, atom3, type, r1, r2, r3, k):
        AbstractAngleType.__init__(self, atom1, atom2, atom3, type)
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.k = k

