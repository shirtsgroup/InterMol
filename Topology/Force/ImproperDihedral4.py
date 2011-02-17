from cctools.Decorators import *
from AbstractDihedral import AbstractDihedral

class ImproperDihedral4(AbstractDihedral):

    @accepts_compatible_units(None, None, None, None, units.degrees, units.kilojoules_per_mole, None)
    def __init__(self, atom1, atom2, atom3, atom4, phi, k, multiplicity):
        """
        """
        AbstractDihedral.__init__(self, atom1, atom2, atom3, atom4)
        self.phi = phi
        self.k = k
        self.multiplicity = multiplicity
    def getForceParameters(self):
        return (self.atom1, self.atom2, self.atom3, self.atom4, self.phi, self.k, self.multiplicity)

    def __repr__(self):
        print self.atom1+'  '+ self.atom2+'  '+ self.atom3+'  '+ self.atom4+'  '+self.phi+'  '+ self.k+'  '+ self.multiplicity

    def __str__(self):
        print self.atom1+'  '+ self.atom2+'  '+ self.atom3+'  '+ self.atom4+'  '+self.phi+'  '+ self.k+'  '+ self.multiplicity




