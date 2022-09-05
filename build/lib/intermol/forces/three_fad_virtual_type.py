import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_3_virtual_type import Abstract3VirtualType


class ThreeFadVirtualType(Abstract3VirtualType):
    __slots__ = ['theta', 'd', 'placeholder']

    @accepts_compatible_units(None, None, None, None, 
                              theta=units.degrees,
                              d=units.nanometers,
                              placeholder=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                 theta=0.0 * units.degrees,
                 d=0.0 * units.nanometers,
                 placeholder=False):
        Abstract3VirtualType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, placeholder)
        self.theta = theta
        self.d = d


class ThreeFadVirtual(ThreeFadVirtualType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, atom4, bondingtype1=None, bondingtype2=None, bondingtype3=None, bondingtype4=None, 
                 theta=0.0 * units.degrees,
                 d=0.0 * units.nanometers,
                 placeholder=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        ThreeFadVirtualType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                theta=theta,
                d=d,
                placeholder=placeholder)