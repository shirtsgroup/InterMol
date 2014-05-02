import sys
sys.path.append('..')

from intermol.Decorators import *
from AbstractDihedralType import *

class ProperPeriodicDihedralType(AbstractDihedralType):

    @accepts_compatible_units(None,
            None,
            None,
            None,
            None,
            units.degrees,
            units.kilojoules_per_mole,
            None)
    def __init__(self, atom1, atom2, atom3, atom4, type, phi, k, multiplicity):
        """
        """
        AbstractDihedralType.__init__(self, atom1, atom2, atom3, atom4, type)
        self.phi = phi
        self.k = k
        self.multiplicity = multiplicity

