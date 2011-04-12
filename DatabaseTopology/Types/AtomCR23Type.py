from Topology.Decorators import *
from AbstractAtomType import AbstractAtomType


class AtomCR23Type(AbstractAtomType):
    @accepts_compatible_units(None,
            None, 
            None,
            units.amu, 
            units.elementary_charge, 
            None, 
            units.nanometers, 
            units.kilojoules_per_mole)

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
                bondtype,
                Z, 
                mass, 
                charge, 
                ptype)
        self.sigma = sigma
        self.epsilon = epsilon


