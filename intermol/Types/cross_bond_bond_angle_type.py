from intermol.decorators import *
from abstract_angle_type import *

class CrossBondBondAngleType(AbstractAngleType):
    @accepts_compatible_units(None,
            None,
            None,
            units.nanometers,
            units.nanometers,
            units.kilojoules_per_mole * units.nanometers**(-2))
    def __init__(self, atom1, atom2, atom3, r1, r2, k):
        AbstractAngleType.__init__(self, atom1, atom2, atom3)
        self.r1 = r1
        self.r2 = r2
        self.k = k

