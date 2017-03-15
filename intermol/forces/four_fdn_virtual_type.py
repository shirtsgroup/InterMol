import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_4_virtual_type import Abstract4VirtualType


class FourFdnVirtualType(Abstract4VirtualType):
    __slots__ = ['a', 'b', 'c', 'placeholder']

    @accepts_compatible_units(None, None, None, None, None, 
                              a=units.dimensionless,
                              b=units.dimensionless,
                              c=units.nanometers,
                              placeholder=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, bondingtype5, 
                 a=0.0 * units.dimensionless,
                 b=0.0 * units.dimensionless,
                 c=0.0 * units.nanometers,
                 placeholder=False):
        Abstract4VirtualType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, bondingtype5, placeholder)
        self.a = a
        self.b = b
        self.c = c


class FourFdnVirtual(FourFdnVirtualType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, atom4, atom5, bondingtype1=None, bondingtype2=None, bondingtype3=None, bondingtype4=None, bondingtype5=None, 
                 a=0.0 * units.dimensionless,
                 b=0.0 * units.dimensionless,
                 c=0.0 * units.nanometers,
                 placeholder=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.atom5 = atom5
        FourFdnVirtualType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, bondingtype5, 
                a=a,
                b=b,
                c=c,
                placeholder=placeholder)