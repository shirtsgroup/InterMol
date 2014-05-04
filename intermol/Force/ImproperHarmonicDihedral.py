from intermol.Decorators import *
from AbstractDihedral import *


class ImproperHarmonicDihedral(AbstractDihedral):

    @accepts_compatible_units(None, None, None, None, units.degrees, units.kilojoules_per_mole * units.radians**(-2))
    def __init__(self, atom1, atom2, atom3, atom4, xi, k):
        """
        """
        AbstractDihedral.__init__(self, atom1, atom2, atom3, atom4)
        self.xi = xi
        self.k = k

    def get_parameters(self):
        return (self.atom1, self.atom2, self.atom3, self.atom4, self.xi, self.k)

    def __repr__(self):
        print self.atom1+'  '+self.atom2+'  '+ self.atom3+'  '+self.atom4+'  '+self.xi+'  '+self.k

    def __str__(self):
        print self.atom1+'  '+self.atom2+'  '+ self.atom3+'  '+self.atom4+'  '+self.xi+'  '+self.k


