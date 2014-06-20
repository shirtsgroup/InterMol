
import sys
sys.path.append('..')
from intermol.decorators import *
from abstract_bond_type import *

class HarmonicBondType(AbstractBondType):
    @accepts_compatible_units(None, None, units.nanometers, units.kilojoules_per_mole * units.nanometers**(-2))
    def __init__(self, atom1, atom2, length, k):
        AbstractBondType.__init__(self, atom1, atom2)
        self.length = length
        self.k = k
