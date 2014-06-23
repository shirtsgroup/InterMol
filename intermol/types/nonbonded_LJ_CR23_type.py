from intermol.decorators import *
from abstract_nonbonded_type import *

class NonbondedLJCR23Type(AbstractNonbondedType):
    @accepts_compatible_units(None, None, None,
            units.nanometers,
            units.kilojoules_per_mole)
    def __init__(self, atom1, atom2, type, sigma, epsilon):
        AbstractNonbondedType.__init__(self, atom1, atom2, type)
        self.sigma = sigma
        self.epsilon = epsilon

