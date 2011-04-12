from Topology.Decorators import *
from Topology.Force.AbstractBond import *

class CubicBond(AbstractBond):

    @accepts_compatible_units(None, 
            None, 
            units.nanometers, 
            units.kilojoules_per_mole * units.nanometers**(-2), 
            units.kilojoules_per_mole * units.nanometers**(-3))
    def __init__(self, atom1, atom2, length, C2, C3):
        """
        """
        AbstractBond.__init__(self, atom1, atom2)
        self.length = length
        self.C2 = C2
        self.C3 = C3

    def getForceParameters(self):
        return (self.atom1, self.atom2, self.length, self.C2, self.C3)    

    def __repr__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.length) +'  '+  str(self.C2) +'  '+ str(self.C3)

    def __str__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.length) +'  '+  str(self.C2) +'  '+ str(self.C3)

