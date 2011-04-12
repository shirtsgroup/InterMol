from Topology.Decorators import *

class Constraint(object):

    @accepts_compatible_units(None, 
            None, 
            units.nanometers)
    def __init__(self, atom1, atom2, length):
        """
        """
        if atom1:
            self.atom1 = atom1
        if atom2:
            self.atom2 = atom2
        self.length = length

    def getForceParameters(self):
        return (self.atom1, self.atom2, self.length)

