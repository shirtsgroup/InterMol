import sys
sys.path.append('..')

from intermol.decorators import *
from abstract_bond_type import *

class G96BondType(AbstractBondType):
    @accepts_compatible_units(None, None, units.nanometers, units.kilojoules_per_mole * units.nanometers**(-4), None, None)
    def __init__(self, atom1, atom2, length, k, order=1, c=False):
        AbstractBondType.__init__(self, atom1, atom2)
        self.length = length
        self.k = k
        self.order = order
        self.c = c #constrained or not, Desmond
