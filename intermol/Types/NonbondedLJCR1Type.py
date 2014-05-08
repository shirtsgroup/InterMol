from intermol.Decorators import *
from AbstractNonbondedType import *

class NonbondedLJCR1Type(AbstractNonbondedType):
    @accepts_compatible_units(None, None, None,
            units.kilojoules_per_mole * units.nanometers**(6),
            units.kilojoules_per_mole * units.nanometers**(12))
    def __init__(self, atom1, atom2, type, sigma, epsilon):
        AbstractNonbondedType.__init__(self, atom1, atom2, type)
        self.sigma = sigma
        self.epsilon = epsilon

