from intermol.decorators import *
from abstract_pair import *

class LJ1PairCR23(AbstractPair):
    __slots__ = ['V', 'W']
    @accepts_compatible_units(None, None,
            units.nanometers,
            units.kilojoules_per_mole)
    def __init__(self, atom1, atom2, V, W):
        """
        """
        AbstractPair.__init__(self, atom1, atom2)
        self.V = V
        self.W = W

    def get_parameters(self):
        return (self.atom1, self.atom2, self.V, self.W)

    def __repr__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.V) +'  '+  str(self.W)

    def __str__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.V) +'  '+  str(self.W)

