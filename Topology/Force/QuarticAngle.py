from Topology.Decorators import *
from Topology.Force.AbstractAngle import *

class QuarticAngle(AbstractAngle):
    
    @accepts_compatible_units(None, 
            None, 
            None, 
            units.degrees, 
            units.kilojoules_per_mole, 
            units.kilojoules_per_mole * units.radians**(-1), 
            units.kilojoules_per_mole * units.radians**(-2), 
            units.kilojoules_per_mole * units.radians**(-3), 
            units.kilojoules_per_mole * units.radians**(-4)) 
    def __init__(self, atom1, atom2, atom3, theta, C0, C1, C2, C3, C4):
        """
        """
        AbstractAngle.__init__(self, atom1, atom2, atom3)
        self.theta = theta
        self.C0 = C0
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4

    def getForceParameters(self):
        return (self.atom1, self.atom2, self.atom3, self.theta, self.C0, self.C1, self.C2, self.C3, self.C4)

    def __repr__(self):
        print self.atom1+'  '+self.atom2+'  '+ self.atom3+'  '+self.theta+'  '+self.C0+'  '+self.C1+'  '+self.C2+'  '+self.C3+'  '+self.C4+'  '+self.k


    def __str__(self):
        print self.atom1+'  '+self.atom2+'  '+ self.atom3+'  '+self.theta+'  '+self.C0+'  '+self.C1+'  '+self.C2+'  '+self.C3+'  '+self.C4+'  '+self.k

           
