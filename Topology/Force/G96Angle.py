from cctools.Decorators import *
from cctools.Force.AbstractAngle import *
      
class G96Angle(AbstractAngle):

    @accepts_compatible_units(None, None, None, units.degrees, units.kilojoules_per_mole)
    def __init__(self, atom1, atom2, atom3, theta, k):
        """
        """
        AbstractAngle.__init__(self, atom1, atom2, atom3)
        self.theta = theta
        self.k = k
    def getForceParameters(self):
        return (self.atom1, self.atom2, self.atom3, self.theta, self.k)

    def __repr__(self):
        print self.atom1+'  '+self.atom2+'  '+ self.atom3+'  '+ self.theta+'  '+self.k


    def __str__(self):
        print self.atom1+'  '+self.atom2+'  '+ self.atom3+'  '+ self.theta+'  '+self.k

