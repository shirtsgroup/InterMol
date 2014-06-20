from intermol.decorators import *
from abstract_angle_type import *

class QuarticAngleType(AbstractAngleType):
    @accepts_compatible_units(None,
            None,
            None,
            units.degrees,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole * units.radians**(-1),
            units.kilojoules_per_mole * units.radians**(-2),
            units.kilojoules_per_mole * units.radians**(-3),
            units.kilojoules_per_mole * units.radians**(-4))
    def __init__(self, atom1, atom2, atom3, theta, C0, C1, C2, C3, C4):
        AbstractAngleType.__init__(self, atom1, atom2, atom3)
        self.theta = theta
        self.C0 = C0
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4
