import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_pair_type import AbstractPairType


class LjSigepsPairType(AbstractPairType):
    __slots__ = ['sigma', 'epsilon', 'scaleLJ', 'scaleQQ', 'long']

    @accepts_compatible_units(None, None, 
                              sigma=units.nanometers,
                              epsilon=units.kilojoules_per_mole,
                              scaleLJ=None,
                              scaleQQ=None,
                              long=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 sigma=0.0 * units.nanometers,
                 epsilon=0.0 * units.kilojoules_per_mole,
                 scaleLJ=None, scaleQQ=None, long=False):
        AbstractPairType.__init__(self, bondingtype1, bondingtype2, scaleLJ, scaleQQ, long)
        self.sigma = sigma
        self.epsilon = epsilon
