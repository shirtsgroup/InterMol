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
