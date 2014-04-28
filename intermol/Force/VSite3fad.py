from intermol.Decorators import *

class VSite3fad(object):

    @accepts_compatible_units(None, None, None, None, units.degrees, units.nanometers)
    def __init__(self, atom1, atom2, atom3, atom4, theta, d):
        if atom1:
            self.atom1 = atom1
        if atom2:
            self.atom2 = atom2
        if atom3:
            self.atom3 = atom3
        if atom4:
            self.atom4 = atom4
        self.theta = theta
        self.d = d

    def get_parameters(self):
        return(self.atom1, self.atom2, self.atom3, self.atom4, self.theta, self.d)

