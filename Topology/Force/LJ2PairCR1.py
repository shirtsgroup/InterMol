from Topology.Decorators import *
from Topology.Force.AbstractPair import *

class LJ2PairCR1(AbstractPair):

    @accepts_compatible_units(None, None, None, units.elementary_charge, units.elementary_charge, units.kilojoules_per_mole * units.nanometers**(6), units.kilojoules_per_mole * units.nanometers**(12))
    def __init__(self, atom1, atom2, fudgeQQ, qi, qj, V, W):
        """
        """
        AbstractPair.__init__(self, atom1, atom2)
        self.fudgeQQ = fudgeQQ
        self.qi = qi
        self.qj = qj
        self.V = V
        self.W = W                

    def getForceParameters(self):
        return (self.atom1, self.atom2, self.fudgeQQ, self.qi, self.qj, self.V, self.W) 

    def __repr__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.fudgeQQ) +'  '+  str(self.qi)+'  '+ str(self.qj) +'   '+str(self.V)+'   '+str(self.W)

    def __str__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.fudgeQQ) +'  '+  str(self.qi)+'  '+ str(self.qj) +'   '+str(self.V)+'   '+str(self.W)
