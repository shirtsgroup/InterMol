import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_angle_type import AbstractAngleType


class UreyBradleyAngleType(AbstractAngleType):
    __slots__ = ['theta', 'k', 'r', 'kUB', 'c']

    @accepts_compatible_units(None, None, None, 
                              theta=units.degrees,
                              k=units.kilojoules_per_mole * units.radians ** (-2),
                              r=units.nanometers,
                              kUB=units.kilojoules_per_mole * units.nanometers ** (-2),
                              c=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, 
                 theta=0.0 * units.degrees,
                 k=0.0 * units.kilojoules_per_mole * units.radians ** (-2),
                 r=0.0 * units.nanometers,
                 kUB=0.0 * units.kilojoules_per_mole * units.nanometers ** (-2),
                 c=False):
        AbstractAngleType.__init__(self, bondingtype1, bondingtype2, bondingtype3, c)
        self.theta = theta
        self.k = k
        self.r = r
        self.kUB = kUB


