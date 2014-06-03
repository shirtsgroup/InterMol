from intermol.Decorators import *
from AbstractBond import *

class CubicBond(AbstractBond):
    __slots__ = ['length', 'C2', 'C3', 'order','c'] 
    @accepts_compatible_units(None,
            None,
            units.nanometers,
            units.kilojoules_per_mole * units.nanometers**(-2),
            units.kilojoules_per_mole * units.nanometers**(-3), None, None)
    def __init__(self, atom1, atom2, length, C2, C3, order=1, c=False): # default bond order is 1
        """
        """
        AbstractBond.__init__(self, atom1, atom2)
        self.length = length
        self.C2 = C2
        self.C3 = C3
    	self.order = order
    	self.c = c #constrained or not, Desmond only

    def get_parameters(self):
        return (self.atom1, self.atom2, self.length, self.C2, self.C3)

    def __repr__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.length) +'  '+  str(self.C2) +'  '+ str(self.C3)

    def __str__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.length) +'  '+  str(self.C2) +'  '+ str(self.C3)

