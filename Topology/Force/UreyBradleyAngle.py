from cctools.Decorators import *
from cctools.Force.AbstractAngle import *

class UreyBradleyAngle(AbstractAngle):

    @accepts_compatible_units(None, None, None, units.degrees, units.kilojoules_per_mole, units.nanometers, units.kilojoules_per_mole)
    def __init__(self, atom1, atom2, atom3, theta, k, r, kUB):
        """
        """
        AbstractAngle.__init__(self, atom1, atom2, atom3)
        self.theta = theta
        self.k = k
        self.r = r
        self.kUB = kUB

    def getForceParameters(self):
        return (self.atom1, self.atom2, self.atom3, self.theta, self.k, self.r, self.kUB)   


    def __repr__(self):
        print self.atom1+'  '+self.atom2+'  '+ self.atom3+'  '+self.theta+'  '+self.k+'  '+self.r+'  '+self.kUB


    def __str__(self):
        print self.atom1+'  '+self.atom2+'  '+ self.atom3+'  '+self.theta+'   '+ self.k+'  '+self.r+'   '+self.kUB

