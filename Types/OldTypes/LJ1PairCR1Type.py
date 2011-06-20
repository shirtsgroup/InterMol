from ctools.Decorators import *
from ctools.Types.AbstractPairType import *

class LJ1PairCR1Type(AbstractPairType):
    @accepts_compatible_units(None, None, None, units.kilojoules_per_mole * units.nanometers**(6), units.kilojoules_per_mole * units.nanometers**(12))
    def __init__(self, atom1, atom2, type, sigma, epsilon):
        AbstractPairType.__init__(self, atom1, atom2, type)
        self.sigma = sigma
        self.epsilon = epsilon

