import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_angle_type import AbstractAngleType


class HarmonicAngleType(AbstractAngleType):
    __slots__ = ['theta', 'k', 'c']

    @accepts_compatible_units(None, None, None, 
                              theta=units.degrees,
                              k=units.kilojoules_per_mole * units.radians **(-2),
                              c=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, 
                 theta=0.0 * units.degrees,
                 k=0.0 * units.kilojoules_per_mole * units.radians **(-2),
                 c=False):
        AbstractAngleType.__init__(self, bondingtype1, bondingtype2, bondingtype3, c)
        self.theta = theta
        self.k = k
