import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_bond_type import AbstractBondType


class MorseBondType(AbstractBondType):
    __slots__ = ['length', 'D', 'beta', 'order', 'c']

    @accepts_compatible_units(None, None, 
                              length=units.nanometers,
                              D=units.kilojoules_per_mole,
                              beta=units.nanometers ** (-1),
                              order=None,
                              c=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 length=0.0 * units.nanometers,
                 D=0.0 * units.kilojoules_per_mole,
                 beta=0.0 * units.nanometers ** (-1),
                 order=1, c=False):
        AbstractBondType.__init__(self, bondingtype1, bondingtype2, order, c)
        self.length = length
        self.D = D
        self.beta = beta


class MorseBond(MorseBondType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 length=0.0 * units.nanometers,
                 D=0.0 * units.kilojoules_per_mole,
                 beta=0.0 * units.nanometers ** (-1),
                 order=1, c=False):
        self.atom1 = atom1
        self.atom2 = atom2
        MorseBondType.__init__(self, bondingtype1, bondingtype2, 
                length=length,
                D=D,
                beta=beta,
                order=order, c=c)