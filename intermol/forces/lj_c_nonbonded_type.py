import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_nonbonded_type import AbstractNonbondedType


class LjCNonbondedType(AbstractNonbondedType):
    __slots__ = ['C6', 'C12', 'type']

    @accepts_compatible_units(None, None, 
                              C6=units.kilojoules_per_mole * units.nanometers ** (6),
                              C12=units.kilojoules_per_mole * units.nanometers ** (12),
                              type=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 C6=0.0 * units.kilojoules_per_mole * units.nanometers ** (6),
                 C12=0.0 * units.kilojoules_per_mole * units.nanometers ** (12)):
        AbstractNonbondedType.__init__(self, bondingtype1, bondingtype2)
        self.C6 = C6
        self.C12 = C12

