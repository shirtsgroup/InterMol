from ctools.Decorators import *

class Constraint3(object):

    @accepts_compatible_units(None, 
            None,
	    None,
            units.nanometers,
	    units.nanometers,
	    units.nanometers,
	    None)
    def __init__(self, atom1, atom2, atom3, length1, length2, length3, type):
        """
        """
	    
        self.atom1 = atom1
        self.atom2 = atom2
	self.atom3 = atom3
        self.length1 = length1
	self.length2 = length2
	self.length3 = length3
	self.type = type

    def getForceParameters(self):

        return (self.atom1, self.atom2, self.atom3, self.length1, self.length2, self.length3, self.type)

