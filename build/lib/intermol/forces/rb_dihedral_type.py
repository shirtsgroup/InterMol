import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_dihedral_type import AbstractDihedralType


class RbDihedralType(AbstractDihedralType):
    __slots__ = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'improper']

    @accepts_compatible_units(None, None, None, None, 
                              C0=units.kilojoules_per_mole,
                              C1=units.kilojoules_per_mole,
                              C2=units.kilojoules_per_mole,
                              C3=units.kilojoules_per_mole,
                              C4=units.kilojoules_per_mole,
                              C5=units.kilojoules_per_mole,
                              C6=units.kilojoules_per_mole,
                              improper=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                 C0=0.0 * units.kilojoules_per_mole,
                 C1=0.0 * units.kilojoules_per_mole,
                 C2=0.0 * units.kilojoules_per_mole,
                 C3=0.0 * units.kilojoules_per_mole,
                 C4=0.0 * units.kilojoules_per_mole,
                 C5=0.0 * units.kilojoules_per_mole,
                 C6=0.0 * units.kilojoules_per_mole,
                 improper=False):
        AbstractDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, improper)
        self.C0 = C0
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4
        self.C5 = C5
        self.C6 = C6


class RbDihedral(RbDihedralType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, atom4, bondingtype1=None, bondingtype2=None, bondingtype3=None, bondingtype4=None, 
                 C0=0.0 * units.kilojoules_per_mole,
                 C1=0.0 * units.kilojoules_per_mole,
                 C2=0.0 * units.kilojoules_per_mole,
                 C3=0.0 * units.kilojoules_per_mole,
                 C4=0.0 * units.kilojoules_per_mole,
                 C5=0.0 * units.kilojoules_per_mole,
                 C6=0.0 * units.kilojoules_per_mole,
                 improper=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        RbDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                C0=C0,
                C1=C1,
                C2=C2,
                C3=C3,
                C4=C4,
                C5=C5,
                C6=C6,
                improper=improper)