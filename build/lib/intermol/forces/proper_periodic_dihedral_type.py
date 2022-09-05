import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_dihedral_type import AbstractDihedralType


class ProperPeriodicDihedralType(AbstractDihedralType):
    __slots__ = ['phi', 'k', 'multiplicity', 'weight', 'improper']

    @accepts_compatible_units(None, None, None, None, 
                              phi=units.degrees,
                              k=units.kilojoules_per_mole,
                              multiplicity=units.dimensionless,
                              weight=units.dimensionless,
                              improper=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                 phi=0.0 * units.degrees,
                 k=0.0 * units.kilojoules_per_mole,
                 multiplicity=0.0 * units.dimensionless,
                 weight=0.0 * units.dimensionless,
                 improper=False):
        AbstractDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, improper)
        self.phi = phi
        self.k = k
        self.multiplicity = multiplicity
        self.weight = weight


class ProperPeriodicDihedral(ProperPeriodicDihedralType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, atom4, bondingtype1=None, bondingtype2=None, bondingtype3=None, bondingtype4=None, 
                 phi=0.0 * units.degrees,
                 k=0.0 * units.kilojoules_per_mole,
                 multiplicity=0.0 * units.dimensionless,
                 weight=0.0 * units.dimensionless,
                 improper=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        ProperPeriodicDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                phi=phi,
                k=k,
                multiplicity=multiplicity,
                weight=weight,
                improper=improper)