import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_bond_type import AbstractBondType


class FeneBondType(AbstractBondType):
    __slots__ = ['length', 'kb', 'order', 'c']

    @accepts_compatible_units(None, None, 
                              length=units.nanometers,
                              kb=units.kilojoules_per_mole * units.nanometers ** (-2),
                              order=None,
                              c=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 length=0.0 * units.nanometers,
                 kb=0.0 * units.kilojoules_per_mole * units.nanometers ** (-2),
                 order=1, c=False):
        AbstractBondType.__init__(self, bondingtype1, bondingtype2, order, c)
        self.length = length
        self.kb = kb


class FeneBond(FeneBondType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 length=0.0 * units.nanometers,
                 kb=0.0 * units.kilojoules_per_mole * units.nanometers ** (-2),
                 order=1, c=False):
        self.atom1 = atom1
        self.atom2 = atom2
        FeneBondType.__init__(self, bondingtype1, bondingtype2, 
                length=length,
                kb=kb,
                order=order, c=c)