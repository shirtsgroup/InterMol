from ctools.Decorators import *

class Constraint2(object):

    @accepts_compatible_units(None, 
            None,
            units.nanometers,
            None)
    def __init__(self, atom1, atom2, length1, type = 'AH1'):
        """
        """
	    
        if atom1:
            self.atom1 = atom1
        if atom2:
            self.atom2 = atom2
        self.length1 = length1
        self.type = type

    def getForceParameters(self):

        return (self.atom1, self.atom2, self.length1)

