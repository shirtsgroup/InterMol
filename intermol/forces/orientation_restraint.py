from intermol.decorators import *

class OrientationRestraint(object):

    @accepts_compatible_units(None, None, None, None, None, None, units.amu, units.amu**(-1))
    ### units for 'c' should actually be units.amu * units.nanometers**(alpha) - unsure of method for implementation
    def __init__(self, atom1, atom2, exp, label, alpha, c, obs, weight):
        """
        """
        if atom1:
            self.atom1 = atom1
        if atom2:
            self.atom2 = atom2
        self.exp = exp
        self.label = label
        self.alpha = alpha
        self.c = c
        self.obs = obs
        self.weight = weight

    def get_parameters(self):
        return (self.atom1, self.atom2, self.exp, self.label, self.alpha, self.c, self.obs, self.weight)

