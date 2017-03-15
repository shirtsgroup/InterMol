import parmed.unit as units

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


class LjqSigepsPair(LjqSigepsPairType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 qi=0.0 * units.elementary_charge,
                 qj=0.0 * units.elementary_charge,
                 sigma=0.0 * units.nanometers,
                 epsilon=0.0 * units.kilojoules_per_mole,
                 scaleLJ=None, scaleQQ=None, long=False):
        self.atom1 = atom1
        self.atom2 = atom2
        LjqSigepsPairType.__init__(self, bondingtype1, bondingtype2, 
                qi=qi,
                qj=qj,
                sigma=sigma,
                epsilon=epsilon,
                scaleLJ=scaleLJ, scaleQQ=scaleQQ, long=long)