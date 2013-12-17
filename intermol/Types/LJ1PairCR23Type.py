import sys
sys.path.append('..')

from intermol.Decorators import *
from AbstractPairType import *

class LJ1PairCR23Type(AbstractPairType):
    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole)
    def __init__(self, atom1, atom2, type, sigma, epsilon):
        AbstractPairType.__init__(self, atom1, atom2, type)
        self.sigma = sigma
        self.epsilon = epsilon

