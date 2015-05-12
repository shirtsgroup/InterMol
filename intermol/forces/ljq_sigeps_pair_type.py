import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_pair_type import AbstractPairType


class LjqSigepsPairType(AbstractPairType):
    __slots__ = ['qi', 'qj', 'sigma', 'epsilon', 'scaleLJ', 'scaleQQ', 'long']

    @accepts_compatible_units(None, None, 
                              qi=units.elementary_charge,
                              qj=units.elementary_charge,
                              sigma=units.nanometers,
                              epsilon=units.kilojoules_per_mole,
                              scaleLJ=None,
                              scaleQQ=None,
                              long=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 qi=0.0 * units.elementary_charge,
                 qj=0.0 * units.elementary_charge,
                 sigma=0.0 * units.nanometers,
                 epsilon=0.0 * units.kilojoules_per_mole,
                 scaleLJ=None, scaleQQ=None, long=False):
        AbstractPairType.__init__(self, bondingtype1, bondingtype2, scaleLJ, scaleQQ, long)
        self.qi = qi
        self.qj = qj
        self.sigma = sigma
        self.epsilon = epsilon
