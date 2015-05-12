import simtk.unit as units

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

