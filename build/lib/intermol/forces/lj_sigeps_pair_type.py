import parmed.unit as units

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


class LjSigepsPair(LjSigepsPairType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 sigma=0.0 * units.nanometers,
                 epsilon=0.0 * units.kilojoules_per_mole,
                 scaleLJ=None, scaleQQ=None, long=False):
        self.atom1 = atom1
        self.atom2 = atom2
        LjSigepsPairType.__init__(self, bondingtype1, bondingtype2, 
                sigma=sigma,
                epsilon=epsilon,
                scaleLJ=scaleLJ, scaleQQ=scaleQQ, long=long)