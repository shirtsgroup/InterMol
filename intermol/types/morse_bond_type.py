from intermol.decorators import *
from abstract_bond_type import *

class MorseBondType(AbstractBondType):
    @accepts_compatible_units(None, None, units.nanometers, units.kilojoules_per_mole,  units.nanometers**(-1))
    def __init__(self, atom1, atom2, length, D, beta):
        AbstractBondType.__init__(self, atom1, atom2)
        self.length = length
        self.D = D
        self.beta = beta
