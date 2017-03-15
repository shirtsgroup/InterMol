import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_pair_type import AbstractPairType


class LjDefaultPairType(AbstractPairType):
    __slots__ = ['scaleLJ', 'scaleQQ', 'long']

    @accepts_compatible_units(None, None, 
                              scaleLJ=None,
                              scaleQQ=None,
                              long=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 scaleLJ=None, scaleQQ=None, long=False):
        AbstractPairType.__init__(self, bondingtype1, bondingtype2, scaleLJ, scaleQQ, long)


class LjDefaultPair(LjDefaultPairType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 scaleLJ=None, scaleQQ=None, long=False):
        self.atom1 = atom1
        self.atom2 = atom2
        LjDefaultPairType.__init__(self, bondingtype1, bondingtype2, 
                scaleLJ=scaleLJ, scaleQQ=scaleQQ, long=long)