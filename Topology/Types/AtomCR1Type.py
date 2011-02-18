from cctools.Decorators import *
from AbstractAtomType import AbstractAtomType

class AtomCR1Type(AbstractAtomType):
    @accepts_compatible_units(None, 
            None,
            None, 
            units.amu, 
            units.elementary_charge, 
            None, 
            units.kilojoules_per_mole * units.nanometers**(6), 
            units.kilojoules_per_mole * units.nanometers**(12))
    def __init__(self, 
            atomtype, 
            bondtype,
            Z, 
            mass, 
            charge, 
            ptype, 
            sigma, 
            epsilon):
        AbstractAtomType.__init__(self, 
            atomtype, 
            Z, 
            mass, 
            charge, 
            ptype)
        self.sigma = sigma
        self.epsilon = epsilon

