import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_dihedral_type import AbstractDihedralType


class BendingTorsionDihedralType(AbstractDihedralType):
    __slots__ = ['a0', 'a1', 'a2', 'a3', 'a4', 'improper']

    @accepts_compatible_units(None, None, None, None, 
                              a0=units.kilojoules_per_mole,
                              a1=units.kilojoules_per_mole,
                              a2=units.kilojoules_per_mole,
                              a3=units.kilojoules_per_mole,
                              a4=units.kilojoules_per_mole,
                              improper=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                 a0=0.0 * units.kilojoules_per_mole,
                 a1=0.0 * units.kilojoules_per_mole,
                 a2=0.0 * units.kilojoules_per_mole,
                 a3=0.0 * units.kilojoules_per_mole,
                 a4=0.0 * units.kilojoules_per_mole,
                 improper=False):
        AbstractDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, improper)
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4


class BendingTorsionDihedral(BendingTorsionDihedralType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, atom4, bondingtype1=None, bondingtype2=None, bondingtype3=None, bondingtype4=None, 
                 a0=0.0 * units.kilojoules_per_mole,
                 a1=0.0 * units.kilojoules_per_mole,
                 a2=0.0 * units.kilojoules_per_mole,
                 a3=0.0 * units.kilojoules_per_mole,
                 a4=0.0 * units.kilojoules_per_mole,
                 improper=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        BendingTorsionDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                a0=a0,
                a1=a1,
                a2=a2,
                a3=a3,
                a4=a4,
                improper=improper)