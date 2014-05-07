import sys
sys.path.append('..')

from intermol.Decorators import *
from AbstractDihedralType import *

class DihedralTrigType(AbstractDihedralType):

    @accepts_compatible_units(None,
            None,
            None,
            None,
            units.degrees,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            None                  
            )

    def __init__(self, atom1, atom2, atom3, atom4, phi, fc0, fc1, fc2, fc3, fc4, fc5, fc6, improper = False):
        """
        """
        AbstractDihedralType.__init__(self, atom1, atom2, atom3, atom4)
        self.phi = phi
        self.fc0 = fc0
        self.fc1 = fc1
        self.fc2 = fc2
        self.fc3 = fc3
        self.fc4 = fc4
        self.fc5 = fc5
        self.fc6 = fc6
        self.improper = improper

    def sum_parameters(self, addterms):
        # addterms is another dihedral.

        compatible = False
        sign = 1
        if self.phi == addterms.phi:
            compatible = True
        else:
            # might not be compatible.  Check to see if they are offset by 180 degrees
            if abs((self.phi - addterms.phi).value_in_unit(units.degrees)) == 180:
                compatible = True
                sign = -1

        if (compatible):
            self.fc0 += addterms.fc0   # zero term never get sign flipped.
            self.fc1 += sign*addterms.fc1
            self.fc2 += sign*addterms.fc2
            self.fc3 += sign*addterms.fc3
            self.fc4 += sign*addterms.fc4
            self.fc5 += sign*addterms.fc5
            self.fc6 += sign*addterms.fc6
        else:
            print "Can't add terms for dihedral_trigs with different phase angles (%4f vs. %4f)" % (self.phi, addterms.phi)
