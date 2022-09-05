import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_bond_type import AbstractBondType


class G96BondType(AbstractBondType):
    __slots__ = ['length', 'k', 'order', 'c']

    @accepts_compatible_units(None, None, 
                              length=units.nanometers,
                              k=units.kilojoules_per_mole * units.nanometers ** (-4),
                              order=None,
                              c=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 length=0.0 * units.nanometers,
                 k=0.0 * units.kilojoules_per_mole * units.nanometers ** (-4),
                 order=1, c=False):
        AbstractBondType.__init__(self, bondingtype1, bondingtype2, order, c)
        self.length = length
        self.k = k


class G96Bond(G96BondType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 length=0.0 * units.nanometers,
                 k=0.0 * units.kilojoules_per_mole * units.nanometers ** (-4),
                 order=1, c=False):
        self.atom1 = atom1
        self.atom2 = atom2
        G96BondType.__init__(self, bondingtype1, bondingtype2, 
                length=length,
                k=k,
                order=order, c=c)