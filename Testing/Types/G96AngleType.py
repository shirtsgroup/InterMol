import sys
sys.path.append('..')

from Decorators import *
from Types.AbstractAngleType import *

class G96AngleType(AbstractAngleType):
    @accepts_compatible_units(None,
            None,
            None,
            None,
            units.degrees,
            units.kilojoules_per_mole)
    def __init__(self, atom1, atom2, atom3, type, theta, k):
        AbstractAngleType.__init__(self, atom1, atom2, atom3, type)
        self.theta = theta
        self.k = k

