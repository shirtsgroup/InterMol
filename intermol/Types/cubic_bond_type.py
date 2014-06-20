import sys
sys.path.append('..')

from intermol.decorators import *
from abstract_bond_type import *

class CubicBondType(AbstractBondType):
    @accepts_compatible_units(None, 
            None, 
            units.nanometers, 
            units.kilojoules_per_mole * units.nanometers**(-2), 
            units.kilojoules_per_mole * units.nanometers**(-3), None, None)
    def __init__(self, atom1, atom2, length, C2, C3, order=1, c=False):
        AbstractBondType.__init__(self, atom1, atom2)
        self.length = length
        self.C2 = C2
        self.C3 = C3
        self.order = order
        self.c = c #constrained or not, Desmond
