import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_3_virtual_type import Abstract3VirtualType


class ThreeFdVirtualType(Abstract3VirtualType):
    __slots__ = ['a', 'd', 'placeholder']

    @accepts_compatible_units(None, None, None, None, 
                              a=units.dimensionless,
                              d=units.nanometers,
                              placeholder=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                 a=0.0 * units.dimensionless,
                 d=0.0 * units.nanometers,
                 placeholder=False):
        Abstract3VirtualType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, placeholder)
        self.a = a
        self.d = d


class ThreeFdVirtual(ThreeFdVirtualType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, atom4, bondingtype1=None, bondingtype2=None, bondingtype3=None, bondingtype4=None, 
                 a=0.0 * units.dimensionless,
                 d=0.0 * units.nanometers,
                 placeholder=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        ThreeFdVirtualType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                a=a,
                d=d,
                placeholder=placeholder)