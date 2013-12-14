from ctools.Decorators import *
from Force.AbstractBond import *

class Morse(AbstractBond):

    @accepts_compatible_units(None, None, units.nanometers, units.kilojoules_per_mole, units.nanometers**(-1))
    def __init__(self, atom1, atom2, length, D, beta):
        """
        """
        AbstractBond.__init__(self, atom1, atom2)
        self.length = length
        self.D = D
        self.beta = beta

    def getForceParameters(self):
        return (self.atom1, self.atom2, self.length, self.D, self.beta)

    def __repr__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.length) +'  '+  str(self.D) +'   '+ str(self.beta)

    def __str__(self):
        return str(self.atom1) +'  '+ str(self.atom2) +'  '+  str(self.length) +'  '+  str(self.D) + '  '+ str(beta)

