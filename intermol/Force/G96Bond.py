from intermol.Decorators import *
from AbstractBond import *

class G96Bond(AbstractBond):
    __slots__ = ['length', 'k', 'order', 'c']
    @accepts_compatible_units(None, None, units.nanometers, units.kilojoules_per_mole * units.nanometers**(-4), None, None)
    def __init__(self, atom1, atom2, length, k, order=1, c=False):  # default bond order is 1
        """
        """
        AbstractBond.__init__(self, atom1, atom2)
        self.length = length
        self.k = k
    	self.order = order
    	self.c = c #constrained or not, Desmond only

    def get_parameters(self):
        return (self.atom1, self.atom2, self.length, self.k)

    def __repr__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.length) +'  '+  str(self.k)

    def __str__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.length) +'  '+  str(self.k)

