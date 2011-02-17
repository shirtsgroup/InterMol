from cctools.Decorators import *
from cctools.Force.AbstractPair import *


class LJNBPairCR23(AbstractPair):

    @accepts_compatible_units(None, None, units.elementary_charge, units.elementary_charge, units.nanometers, units.kilojoules_per_mole)
    def __init__(self, atom1, atom2, qi, qj, V, W):
        """
        """
        AbstractPair.__init__(self, atom1, atom2)
        self.qi = qi
        self.qj = qj
        self.V = V
        self.W = W 

    def getForceParameters(self):
        return (self.atom1, self.atom2, self.qi, self.qj, self.V, self.W) 

    def __repr__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.qi)+'  '+ str(self.qj) +'   '+str(self.V)+'   '+str(self.W)

    def __str__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.qi)+'  '+ str(self.qj) +'   '+str(self.V)+'   '+str(self.W)

