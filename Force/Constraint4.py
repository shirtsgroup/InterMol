from ctools.Decorators import *

class Constraint4(object):

    @accepts_compatible_units(None, 
            None,
	    None,
	    None,
            units.nanometers,
	    units.nanometers,
	    units.nanometers,
            None)
    def __init__(self, atom1, atom2, atom3, atom4, length1, length2, length3, type = 'AH3'):
        """
        """
	    
        self.atom1 = atom1
        self.atom2 = atom2
	self.atom3 = atom3
	self.atom4 = atom4
        self.length1 = length1
	self.length2 = length2
	self.length3 = length3
        self.type = type
	
    def getForceParameters(self):

        return (self.atom1, self.atom2, self.atom3, self.atom4, self.length1, self.length2, self.length3)

