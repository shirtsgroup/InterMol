import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_pair_type import AbstractPairType


class LjqDefaultPairType(AbstractPairType):
    __slots__ = ['scaleLJ', 'scaleQQ', 'long']

    @accepts_compatible_units(None, None, 
                              scaleLJ=None,
                              scaleQQ=None,
                              long=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 scaleLJ=None, scaleQQ=None, long=False):
        AbstractPairType.__init__(self, bondingtype1, bondingtype2, scaleLJ, scaleQQ, long)


class LjqDefaultPair(LjqDefaultPairType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 scaleLJ=None, scaleQQ=None, long=False):
        self.atom1 = atom1
        self.atom2 = atom2
        LjqDefaultPairType.__init__(self, bondingtype1, bondingtype2, 
                scaleLJ=scaleLJ, scaleQQ=scaleQQ, long=long)