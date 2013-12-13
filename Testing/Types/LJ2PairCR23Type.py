import sys
sys.path.append('..')

from Decorators import *
from Types.AbstractPairType import *

class LJ2PairCR23Type(AbstractPairType):
    @accepts_compatible_units(None, None, None, None, units.elementary_charge, units.elementary_charge, units.nanometers, units.kilojoules_per_mole)
    def __init__(self, atom1, atom2, type, fudgeQQ, qi, qj, sigma, epsilon):
        AbstractPairType.__init__(self, atom1, atom2, type)
        self.fudgeQQ = fudgeQQ
        self.qi = qi
        self.qj = qj
        self.sigma = sigma
        self.epsilon = epsilon

