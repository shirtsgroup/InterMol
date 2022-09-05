import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_dihedral_type import AbstractDihedralType


class RestrictedBendingDihedralType(AbstractDihedralType):
    __slots__ = ['theta', 'k', 'improper']

    @accepts_compatible_units(None, None, None, None, 
                              theta=units.degrees,
                              k=units.kilojoules_per_mole,
                              improper=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                 theta=0.0 * units.degrees,
                 k=0.0 * units.kilojoules_per_mole,
                 improper=False):
        AbstractDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, improper)
        self.theta = theta
        self.k = k


class RestrictedBendingDihedral(RestrictedBendingDihedralType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, atom4, bondingtype1=None, bondingtype2=None, bondingtype3=None, bondingtype4=None, 
                 theta=0.0 * units.degrees,
                 k=0.0 * units.kilojoules_per_mole,
                 improper=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        RestrictedBendingDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                theta=theta,
                k=k,
                improper=improper)