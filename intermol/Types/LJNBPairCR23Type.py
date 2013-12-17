import sys
sys.path.append('..')

from intermol.Decorators import *
from AbstractPairType import *

class LJNBPairCR23Type(AbstractPairType):
    @accepts_compatible_units(None, None, None, units.elementary_charge, units.elementary_charge, units.nanometers, units.kilojoules_per_mole)
    def __init__(self, atom1, atom2, type, qi, qj, sigma, epsilon):
        AbstractPairType.__init__(self, atom1, atom2, type)
        self.qi = qi
        self.qj = qj
        self.sigma = sigma
        self.epsilon = epsilon
