import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_bond_type import AbstractBondType


class FeneExpandableBondType(AbstractBondType):
    __slots__ = ['k', 'length', 'epsilon', 'sigma', 'delta', 'order', 'c']

    @accepts_compatible_units(None, None, 
                              k=units.kilojoules_per_mole * units.nanometers ** (-2),
                              length=units.nanometers,
                              epsilon=units.kilojoules_per_mole,
                              sigma=units.nanometers,
                              delta=units.nanometers,
                              order=None,
                              c=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 k=0.0 * units.kilojoules_per_mole * units.nanometers ** (-2),
                 length=0.0 * units.nanometers,
                 epsilon=0.0 * units.kilojoules_per_mole,
                 sigma=0.0 * units.nanometers,
                 delta=0.0 * units.nanometers,
                 order=1, c=False):
        AbstractBondType.__init__(self, bondingtype1, bondingtype2, order, c)
        self.k = k
        self.length = length
        self.epsilon = epsilon
        self.sigma = sigma
        self.delta = delta


class FeneExpandableBond(FeneExpandableBondType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 k=0.0 * units.kilojoules_per_mole * units.nanometers ** (-2),
                 length=0.0 * units.nanometers,
                 epsilon=0.0 * units.kilojoules_per_mole,
                 sigma=0.0 * units.nanometers,
                 delta=0.0 * units.nanometers,
                 order=1, c=False):
        self.atom1 = atom1
        self.atom2 = atom2
        FeneExpandableBondType.__init__(self, bondingtype1, bondingtype2, 
                k=k,
                length=length,
                epsilon=epsilon,
                sigma=sigma,
                delta=delta,
                order=order, c=c)