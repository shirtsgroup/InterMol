import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_3_virtual_type import Abstract3VirtualType


class ThreeLinearVirtualType(Abstract3VirtualType):
    __slots__ = ['a', 'b', 'placeholder']

    @accepts_compatible_units(None, None, None, None, 
                              a=units.dimensionless,
                              b=units.dimensionless,
                              placeholder=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                 a=0.0 * units.dimensionless,
                 b=0.0 * units.dimensionless,
                 placeholder=False):
        Abstract3VirtualType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, placeholder)
        self.a = a
        self.b = b


class ThreeLinearVirtual(ThreeLinearVirtualType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, atom4, bondingtype1=None, bondingtype2=None, bondingtype3=None, bondingtype4=None, 
                 a=0.0 * units.dimensionless,
                 b=0.0 * units.dimensionless,
                 placeholder=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        ThreeLinearVirtualType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                a=a,
                b=b,
                placeholder=placeholder)