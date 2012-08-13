from ctools.Decorators import *

class Constraint6(object):

    @accepts_compatible_units(None, 
            None,
	    None,
	    None,
	    None,
	    None,
            units.nanometers,
	    units.nanometers,
	    units.nanometers,
	    units.nanometers,
	    units.nanometers,
            None)
    def __init__(self, atom1, atom2, atom3, atom4, atom5, atom6, length1, length2, length3, length4, length5, type = 'AH5'):
        """
        """
	    
        self.atom1 = atom1
        self.atom2 = atom2
	self.atom3 = atom3
	self.atom4 = atom4
	self.atom5 = atom5
	self.atom6 = atom6
        self.length1 = length1
	self.length2 = length2
	self.length3 = length3
	self.length4 = length4
	self.length5 = length5
        self.type = type
	
    def getForceParameters(self):

        return (self.atom1, self.atom2, self.atom3, self.atom4, self.atom5, self.atom6, self.length1, self.length2, self.length3, self.length4, self.length5)

