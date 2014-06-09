from intermol.Decorators import *
from AbstractDihedralType import *

class ImproperHarmonicDihedralType(AbstractDihedralType):

    @accepts_compatible_units(None,
            None,
            None,
            None,
            units.degrees,
            units.kilojoules_per_mole * units.radians**(-2))
    def __init__(self, atom1, atom2, atom3, atom4, xi, k):
        """
        """
        AbstractDihedralType.__init__(self, atom1, atom2, atom3, atom4, improper=True)
        self.xi = xi
        self.k = k

