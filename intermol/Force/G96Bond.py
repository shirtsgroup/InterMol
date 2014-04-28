from intermol.Decorators import *
from AbstractBond import *

class G96Bond(AbstractBond):

    @accepts_compatible_units(None, None, units.nanometers, units.kilojoules_per_mole * units.nanometers**(-4))
    def __init__(self, atom1, atom2, length, k):
        """
        """
        AbstractBond.__init__(self, atom1, atom2)
        self.length = length
        self.k = k

    def get_parameters(self):
        return (self.atom1, self.atom2, self.length, self.k)

    def __repr__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.length) +'  '+  str(self.k)

    def __str__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.length) +'  '+  str(self.k)

