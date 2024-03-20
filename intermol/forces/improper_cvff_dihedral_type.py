import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_dihedral_type import AbstractDihedralType


class ImproperCvffDihedralType(AbstractDihedralType):
    __slots__ = ['k', 'sign', 'multiplicity']

    @accepts_compatible_units(None, None, None, None, 
                              k=units.kilojoules_per_mole,
                              sign=units.dimensionless,
                              multiplicity=units.dimensionless,)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                 k=0.0 * units.kilojoules_per_mole,
                 sign=1.0 * units.dimensionless,
                 multiplicity=0.0 * units.dimensionless,
                 ):
        AbstractDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4)
        self.k = k
        self.sign = sign
        self.multiplicity = multiplicity


class ImproperCvffDihedral(ImproperCvffDihedralType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, atom4, bondingtype1=None, bondingtype2=None, bondingtype3=None, bondingtype4=None, 
                 k=0.0 * units.kilojoules_per_mole,
                 sign=1.0 * units.dimensionless,
                 multiplicity=0.0 * units.dimensionless):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        ImproperCvffDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                k=k,
                sign=sign,
                multiplicity=multiplicity)