import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_angle_type import AbstractAngleType


class RestrictedBendingAngleType(AbstractAngleType):
    __slots__ = ['theta', 'k', 'c']

    @accepts_compatible_units(None, None, None, 
                              theta=units.degrees,
                              k=units.kilojoules_per_mole,
                              c=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, 
                 theta=0.0 * units.degrees,
                 k=0.0 * units.kilojoules_per_mole,
                 c=False):
        AbstractAngleType.__init__(self, bondingtype1, bondingtype2, bondingtype3, c)
        self.theta = theta
        self.k = k


class RestrictedBendingAngle(RestrictedBendingAngleType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, bondingtype1=None, bondingtype2=None, bondingtype3=None, 
                 theta=0.0 * units.degrees,
                 k=0.0 * units.kilojoules_per_mole,
                 c=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        RestrictedBendingAngleType.__init__(self, bondingtype1, bondingtype2, bondingtype3, 
                theta=theta,
                k=k,
                c=c)