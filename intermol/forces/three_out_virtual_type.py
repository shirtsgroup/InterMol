import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_3_virtual_type import Abstract3VirtualType


class ThreeOutVirtualType(Abstract3VirtualType):
    __slots__ = ['a', 'b', 'c', 'placeholder']

    @accepts_compatible_units(None, None, None, None, 
                              a=units.dimensionless,
                              b=units.dimensionless,
                              c=units.nanometers ** (-1),
                              placeholder=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                 a=0.0 * units.dimensionless,
                 b=0.0 * units.dimensionless,
                 c=0.0 * units.nanometers ** (-1),
                 placeholder=False):
        Abstract3VirtualType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, placeholder)
        self.a = a
        self.b = b
        self.c = c

