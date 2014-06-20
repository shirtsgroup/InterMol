from intermol.decorators import *
from abstract_dihedral import *
import pdb


class RBDihedral(AbstractDihedral):
    @accepts_compatible_units(None,
            None,
            None,
            None,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            None)
    def __init__(self, atom1, atom2, atom3, atom4, C0, C1, C2, C3, C4, C5, C6, i=0):
        """
        """
        AbstractDihedral.__init__(self, atom1, atom2, atom3, atom4, -1)
        self.C0 = C0
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4
        self.C5 = C5
        self.C6 = C6
        self.i = i #improper or not--Desmond only

    def get_parameters(self):
        return (self.atom1, self.atom2, self.atom3, self.atom4, self.C0, self.C1, self.C2, self.C3, self.C4, self.C5, self.C6)

    def __repr__(self):
        return "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}".format(self.atom1,
                self.atom2,
                self.atom3,
                self.atom4,
                self.C0,
                self.C1,
                self.C2,
                self.C3,
                self.C4,
                self.C5,
                self.C6)

    def __str__(self):
        return "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}".format(self.atom1,
                self.atom2,
                self.atom3,
                self.atom4,
                self.C0,
                self.C1,
                self.C2,
                self.C3,
                self.C4,
                self.C5,
                self.C6)


