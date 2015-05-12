import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_dihedral_type import AbstractDihedralType


class FourierDihedralType(AbstractDihedralType):
    __slots__ = ['c1', 'c2', 'c3', 'c4', 'c5', 'improper']

    @accepts_compatible_units(None, None, None, None, 
                              c1=units.kilojoules_per_mole,
                              c2=units.kilojoules_per_mole,
                              c3=units.kilojoules_per_mole,
                              c4=units.kilojoules_per_mole,
                              c5=units.kilojoules_per_mole,
                              improper=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                 c1=0.0 * units.kilojoules_per_mole,
                 c2=0.0 * units.kilojoules_per_mole,
                 c3=0.0 * units.kilojoules_per_mole,
                 c4=0.0 * units.kilojoules_per_mole,
                 c5=0.0 * units.kilojoules_per_mole,
                 improper=False):
        AbstractDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, improper)
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4
        self.c5 = c5
