import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_pair_type import AbstractPairType


class LjCPairType(AbstractPairType):
    __slots__ = ['C6', 'C12', 'scaleLJ', 'scaleQQ', 'long']

    @accepts_compatible_units(None, None, 
                              C6=units.kilojoules_per_mole * units.nanometers ** (6),
                              C12=units.kilojoules_per_mole * units.nanometers ** (12),
                              scaleLJ=None,
                              scaleQQ=None,
                              long=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 C6=0.0 * units.kilojoules_per_mole * units.nanometers ** (6),
                 C12=0.0 * units.kilojoules_per_mole * units.nanometers ** (12),
                 scaleLJ=None, scaleQQ=None, long=False):
        AbstractPairType.__init__(self, bondingtype1, bondingtype2, scaleLJ, scaleQQ, long)
        self.C6 = C6
        self.C12 = C12


class LjCPair(LjCPairType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 C6=0.0 * units.kilojoules_per_mole * units.nanometers ** (6),
                 C12=0.0 * units.kilojoules_per_mole * units.nanometers ** (12),
                 scaleLJ=None, scaleQQ=None, long=False):
        self.atom1 = atom1
        self.atom2 = atom2
        LjCPairType.__init__(self, bondingtype1, bondingtype2, 
                C6=C6,
                C12=C12,
                scaleLJ=scaleLJ, scaleQQ=scaleQQ, long=long)