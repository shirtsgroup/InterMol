from intermol.decorators import *


class TorsionTorsionCMAP(object):
    @accepts_compatible_units(None, None, None, None, None, None, None, None, None, None)
    def __init__(self, atom1, atom2, atom3, atom4, atom5, atom6, atom7, atom8, type, chart):
        """
        """
        self.type = type
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.atom5 = atom5
        self.atom6 = atom6
        self.atom7 = atom7
        self.atom8 = atom8
        self.chart = chart

    def getparameters(self):
        return self.atom1, self.atom2, self.atom3, self.atom4, self.atom5, self.atom6, self.atom7, self.atom8, self.type, self.chart