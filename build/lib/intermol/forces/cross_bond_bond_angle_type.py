import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_angle_type import AbstractAngleType


class CrossBondBondAngleType(AbstractAngleType):
    __slots__ = ['r1', 'r2', 'k', 'c']

    @accepts_compatible_units(None, None, None, 
                              r1=units.nanometers,
                              r2=units.nanometers,
                              k=units.kilojoules_per_mole * units.nanometers ** (-2),
                              c=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, 
                 r1=0.0 * units.nanometers,
                 r2=0.0 * units.nanometers,
                 k=0.0 * units.kilojoules_per_mole * units.nanometers ** (-2),
                 c=False):
        AbstractAngleType.__init__(self, bondingtype1, bondingtype2, bondingtype3, c)
        self.r1 = r1
        self.r2 = r2
        self.k = k


class CrossBondBondAngle(CrossBondBondAngleType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, bondingtype1=None, bondingtype2=None, bondingtype3=None, 
                 r1=0.0 * units.nanometers,
                 r2=0.0 * units.nanometers,
                 k=0.0 * units.kilojoules_per_mole * units.nanometers ** (-2),
                 c=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        CrossBondBondAngleType.__init__(self, bondingtype1, bondingtype2, bondingtype3, 
                r1=r1,
                r2=r2,
                k=k,
                c=c)