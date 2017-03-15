import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_angle_type import AbstractAngleType


class CosineAngleType(AbstractAngleType):
    __slots__ = ['k', 'c']

    @accepts_compatible_units(None, None, None, 
                              k=units.kilojoules_per_mole,
                              c=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, 
                 k=0.0 * units.kilojoules_per_mole,
                 c=False):
        AbstractAngleType.__init__(self, bondingtype1, bondingtype2, bondingtype3, c)
        self.k = k


class CosineAngle(CosineAngleType):
    """
    http://lammps.sandia.gov/doc/angle_cosine.html
    """
    def __init__(self, atom1, atom2, atom3, bondingtype1=None, bondingtype2=None, bondingtype3=None, 
                 k=0.0 * units.kilojoules_per_mole,
                 c=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        CosineAngleType.__init__(self, bondingtype1, bondingtype2, bondingtype3, 
                k=k,
                c=c)