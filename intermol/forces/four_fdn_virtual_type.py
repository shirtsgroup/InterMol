import simtk.unit as units

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
