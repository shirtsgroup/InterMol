import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_angle_type import AbstractAngleType


class CosineAngleType(AbstractAngleType):
    __slots__ = ['k', 'c']

    @accepts_compatible_units(None, None, None, 
                              k=units.kilojoules_per_mole,
                              c=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, 
                 k=0.0 * units.kilojoules_per_mole,
                 c=False):
        AbstractAngleType.__init__(self, bondingtype1, bondingtype2, bondingtype3, c)
        self.k = k
