import sys
sys.path.append('..')

from intermol.Decorators import *
from AbstractDihedralType import *

class ImproperDihedral2Type(AbstractDihedralType):

    @accepts_compatible_units(None,
            None,
            None,
            None,
            None,
            units.degrees,
            units.kilojoules_per_mole * units.radians**(-2))
    def __init__(self, atom1, atom2, atom3, atom4, type, xi, k):
        """
        """
        AbstractDihedralType.__init__(self, atom1, atom2, atom3, atom4, type)
        self.xi = xi
        self.k = k

