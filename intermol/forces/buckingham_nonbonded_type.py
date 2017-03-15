import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_nonbonded_type import AbstractNonbondedType


class BuckinghamNonbondedType(AbstractNonbondedType):
    __slots__ = ['a', 'b', 'C6', 'type']

    @accepts_compatible_units(None, None, 
                              a=units.kilojoules_per_mole,
                              b=units.nanometers ** (-1),
                              C6=units.kilojoules_per_mole * units.nanometers ** (6),
                              type=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 a=0.0 * units.kilojoules_per_mole,
                 b=0.0 * units.nanometers ** (-1),
                 C6=0.0 * units.kilojoules_per_mole * units.nanometers ** (6),
                 type=False):
        AbstractNonbondedType.__init__(self, bondingtype1, bondingtype2, type)
        self.a = a
        self.b = b
        self.C6 = C6


class BuckinghamNonbonded(BuckinghamNonbondedType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 a=0.0 * units.kilojoules_per_mole,
                 b=0.0 * units.nanometers ** (-1),
                 C6=0.0 * units.kilojoules_per_mole * units.nanometers ** (6),
                 type=False):
        self.atom1 = atom1
        self.atom2 = atom2
        BuckinghamNonbondedType.__init__(self, bondingtype1, bondingtype2, 
                a=a,
                b=b,
                C6=C6,
                type=type)