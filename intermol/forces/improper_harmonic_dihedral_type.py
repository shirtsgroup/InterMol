import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_dihedral_type import AbstractDihedralType


class ImproperHarmonicDihedralType(AbstractDihedralType):
    __slots__ = ['xi', 'k', 'improper']

    @accepts_compatible_units(None, None, None, None, 
                              xi=units.degrees,
                              k=units.kilojoules_per_mole * units.radians **(-2),
                              improper=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                 xi=0.0 * units.degrees,
                 k=0.0 * units.kilojoules_per_mole * units.radians **(-2),
                 improper=False):
        AbstractDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, improper)
        self.xi = xi
        self.k = k

