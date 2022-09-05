import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_angle_type import AbstractAngleType


class QuarticAngleType(AbstractAngleType):
    __slots__ = ['theta', 'C0', 'C1', 'C2', 'C3', 'C4', 'c']

    @accepts_compatible_units(None, None, None, 
                              theta=units.degrees,
                              C0=units.kilojoules_per_mole,
                              C1=units.kilojoules_per_mole * units.radians ** (-1),
                              C2=units.kilojoules_per_mole * units.radians ** (-2),
                              C3=units.kilojoules_per_mole * units.radians ** (-3),
                              C4=units.kilojoules_per_mole * units.radians ** (-4),
                              c=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, 
                 theta=0.0 * units.degrees,
                 C0=0.0 * units.kilojoules_per_mole,
                 C1=0.0 * units.kilojoules_per_mole * units.radians ** (-1),
                 C2=0.0 * units.kilojoules_per_mole * units.radians ** (-2),
                 C3=0.0 * units.kilojoules_per_mole * units.radians ** (-3),
                 C4=0.0 * units.kilojoules_per_mole * units.radians ** (-4),
                 c=False):
        AbstractAngleType.__init__(self, bondingtype1, bondingtype2, bondingtype3, c)
        self.theta = theta
        self.C0 = C0
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4


class QuarticAngle(QuarticAngleType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, bondingtype1=None, bondingtype2=None, bondingtype3=None, 
                 theta=0.0 * units.degrees,
                 C0=0.0 * units.kilojoules_per_mole,
                 C1=0.0 * units.kilojoules_per_mole * units.radians ** (-1),
                 C2=0.0 * units.kilojoules_per_mole * units.radians ** (-2),
                 C3=0.0 * units.kilojoules_per_mole * units.radians ** (-3),
                 C4=0.0 * units.kilojoules_per_mole * units.radians ** (-4),
                 c=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        QuarticAngleType.__init__(self, bondingtype1, bondingtype2, bondingtype3, 
                theta=theta,
                C0=C0,
                C1=C1,
                C2=C2,
                C3=C3,
                C4=C4,
                c=c)