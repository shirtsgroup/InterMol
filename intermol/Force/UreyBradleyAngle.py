from intermol.Decorators import *
from AbstractAngle import *

class UreyBradleyAngle(AbstractAngle):

    @accepts_compatible_units(None, None, None, units.degrees, units.kilojoules_per_mole * units.radians**(-2), units.nanometers, units.kilojoules_per_mole,None)
    def __init__(self, atom1, atom2, atom3, theta, k, r, kUB, c=False):
        """
        """
        AbstractAngle.__init__(self, atom1, atom2, atom3)
        self.theta = theta
        self.k = k
        self.r = r
        self.kUB = kUB
        self.c = c

    def get_parameters(self):
        return (self.atom1, self.atom2, self.atom3, self.theta, self.k, self.r, self.kUB, self.c)   


    def __repr__(self):
        print self.atom1+'  '+self.atom2+'  '+ self.atom3+'  '+self.theta+'  '+self.k+'  '+self.r+'  '+self.kUB+'  '+self.c  


    def __str__(self):
        print self.atom1+'  '+self.atom2+'  '+ self.atom3+'  '+self.theta+'  '+self.k+'  '+self.r+'  '+self.kUB+'  '+self.c

