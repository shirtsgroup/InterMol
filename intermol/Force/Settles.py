from intermol.Decorators import *

class Settles(object):
    
    @accepts_compatible_units(None, units.nanometers, units.nanometers)
    def __init__(self, atom1, dOH, dHH):
        """
        """
        if atom1:
            self.atom1 = atom1
        self.dOH = dOH
        self.dHH = dHH

    def getarameters(self):
        return (self.atom1, dOH, dHH)

    def __repr__(self):
        print self.atom1+'  '+self.dOH+'  '+self.dHH         

    def __str__(self):
        print self.atom1+'  '+self.dOH+'  '+self.dHH

