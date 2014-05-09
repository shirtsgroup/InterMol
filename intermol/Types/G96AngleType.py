from intermol.Decorators import *
from AbstractAngleType import *

class G96AngleType(AbstractAngleType):
    @accepts_compatible_units(None,
            None,
            None,
            units.degrees,
            units.kilojoules_per_mole)
    def __init__(self, atom1, atom2, atom3, theta, k):
        AbstractAngleType.__init__(self, atom1, atom2, atom3)
        self.theta = theta
        self.k = k

