import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_2_virtual_type import Abstract2VirtualType


class TwoVirtualType(Abstract2VirtualType):
    __slots__ = ['a', 'placeholder']

    @accepts_compatible_units(None, None, None, 
                              a=units.dimensionless,
                              placeholder=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, 
                 a=0.0 * units.dimensionless,
                 placeholder=False):
        Abstract2VirtualType.__init__(self, bondingtype1, bondingtype2, bondingtype3, placeholder)
        self.a = a


class TwoVirtual(TwoVirtualType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, bondingtype1=None, bondingtype2=None, bondingtype3=None, 
                 a=0.0 * units.dimensionless,
                 placeholder=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        TwoVirtualType.__init__(self, bondingtype1, bondingtype2, bondingtype3,
                a=a,
                placeholder=placeholder)
