from ctools.Decorators import *
from ctools.Force.AbstractBond import *

class Bond(AbstractBond):

    @accepts_compatible_units(None, 
            None, 
            units.nanometers, 
            units.kilojoules_per_mole * units.nanometers**(-2))
    def __init__(self, atom1, atom2, length, k):
        """
        """
        AbstractBond.__init__(self, atom1, atom2)
        self.length = length
        self.k = k

    def getForceParameters(self):
        return (self.atom1, self.atom2, self.length, self.k)

    def __repr__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.length) +'  '+  str(self.k)
    
    def __str__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.length) +'  '+  str(self.k)

