import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_angle_type import AbstractAngleType


class UreyBradleyNoharmAngleType(AbstractAngleType):
    __slots__ = ['r', 'kUB', 'c']

    @accepts_compatible_units(None, None, None, 
                              r=units.nanometers,
                              kUB=units.kilojoules_per_mole * units.nanometers ** (-2),
                              c=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, 
                 r=0.0 * units.nanometers,
                 kUB=0.0 * units.kilojoules_per_mole * units.nanometers ** (-2),
                 c=False):
        AbstractAngleType.__init__(self, bondingtype1, bondingtype2, bondingtype3, c)
        self.r = r
        self.kUB = kUB


class UreyBradleyNoharmAngle(UreyBradleyNoharmAngleType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, bondingtype1=None, bondingtype2=None, bondingtype3=None, 
                 r=0.0 * units.nanometers,
                 kUB=0.0 * units.kilojoules_per_mole * units.nanometers ** (-2),
                 c=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        UreyBradleyNoharmAngleType.__init__(self, bondingtype1, bondingtype2, bondingtype3, 
                r=r,
                kUB=kUB,
                c=c)