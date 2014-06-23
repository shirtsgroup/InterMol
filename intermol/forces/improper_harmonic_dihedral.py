from intermol.decorators import *
from abstract_dihedral import *

class ImproperHarmonicDihedral(AbstractDihedral):
    """
    """
    @accepts_compatible_units(None, None, None, None,
            units.degrees,
            units.kilojoules_per_mole * units.radians**(-2))
    def __init__(self, atom1, atom2, atom3, atom4, xi, k):
        AbstractDihedral.__init__(self, atom1, atom2, atom3, atom4, improper=True)
        self.xi = xi
        self.k = k

    def get_parameters(self):
        return (self.atom1, self.atom2, self.atom3, self.atom4, self.xi, self.k)

    def __repr__(self):
        return "{0}, {1}, {2}, {3}: {4}, {5}".format(
                self.atom1, self.atom2, self.atom3, self.atom4,
                self.xi, self.k)

    def __str__(self):
         return "{0}, {1}, {2}, {3}: {4}, {5}".format(
                self.atom1, self.atom2, self.atom3, self.atom4,
                self.xi, self.k)
