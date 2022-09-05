import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_bond_type import AbstractBondType


class ConnectionBondType(AbstractBondType):
    __slots__ = ['order', 'c']

    @accepts_compatible_units(None, None, 
                              order=None,
                              c=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 order=1, c=False):
        AbstractBondType.__init__(self, bondingtype1, bondingtype2, order, c)


class ConnectionBond(ConnectionBondType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 order=1, c=False):
        self.atom1 = atom1
        self.atom2 = atom2
        ConnectionBondType.__init__(self, bondingtype1, bondingtype2, 
                order=order, c=c)