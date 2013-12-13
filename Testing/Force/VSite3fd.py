import sys
sys.path.append('..')

from Decorators import *

class VSite3fd(object):

    @accepts_compatible_units(None, None, None, None, None, units.nanometers)
    def __init__(self, atom1, atom2, atom3, atom4, a, d):
        if atom1:
            self.atom1 = atom1
        if atom2:
            self.atom2 = atom2
        if atom3:
            self.atom3 = atom3
        if atom4:
            self.atom4 = atom4
        self.a = a
        self.d = d

    def getForceParameters(self):
        return(self.atom1, self.atom2, self.atom3, self.atom4, self.a, self.d)

