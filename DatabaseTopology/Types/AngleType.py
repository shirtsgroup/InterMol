from Topology.Decorators import *
from Topology.Types.AbstractAngleType import *

class AngleType(AbstractAngleType):
    @accepts_compatible_units(None,
            None,
            None,
            None,
            units.degrees,
            units.kilojoules_per_mole * units.radians**(-2))
    def __init__(self, atom1, atom2, atom3, type, theta, k):
        """Angle Type blah blah blah
        """
        AbstractAngleType.__init__(self, atom1, atom2, atom3, type)
        self.theta = theta
        self.k = k

