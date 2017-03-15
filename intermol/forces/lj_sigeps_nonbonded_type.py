import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_nonbonded_type import AbstractNonbondedType


class LjSigepsNonbondedType(AbstractNonbondedType):
    __slots__ = ['sigma', 'epsilon', 'type']

    @accepts_compatible_units(None, None, 
                              sigma=units.nanometers,
                              epsilon=units.kilojoules_per_mole,
                              type=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 sigma=0.0 * units.nanometers,
                 epsilon=0.0 * units.kilojoules_per_mole,
                 type=False):
        AbstractNonbondedType.__init__(self, bondingtype1, bondingtype2, type)
        self.sigma = sigma
        self.epsilon = epsilon


class LjSigepsNonbonded(LjSigepsNonbondedType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 sigma=0.0 * units.nanometers,
                 epsilon=0.0 * units.kilojoules_per_mole,
                 type=False):
        self.atom1 = atom1
        self.atom2 = atom2
        LjSigepsNonbondedType.__init__(self, bondingtype1, bondingtype2, 
                sigma=sigma,
                epsilon=epsilon,
                type=type)