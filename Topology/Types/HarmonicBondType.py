from Topology.Decorators import *
from Topology.Types.AbstractBondType import *

class HarmonicBondType(AbstractBondType):
    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole * units.nanometers**(-2))
    def __init__(self, atom1, atom2, type, length, k):
        AbstractBondType.__init__(self, atom1, atom2, type)
        self.length = length
        self.k = k
