from intermol.decorators import *

class AngleRestraint(object):
    @accepts_compatible_units(None, None, None, None,
            units.degrees, units.kilojoules_per_mole, None)
    def __init__(self, atom1, atom2, atom3, atom4, theta, k, multiplicity):
        """
        """
        if atom1:
            self.atom1 = atom1
        if atom2:
            self.atom2 = atom2
        if atom3:
            self.atom3 = atom3
        if atom4:
            self.atom4 = atom4
        self.theta = theta
        self.k = k
        self.multiplicity = multiplicity

    def get_parameters(self):
        return (self. atom1, self.atom2, self.atom3, self.atom4, self.theta, self.k, self.multiplicity)
