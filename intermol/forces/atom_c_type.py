from intermol.decorators import *
from abstract_atom_type import *

class AtomCType(AbstractAtomType):
    @accepts_compatible_units(None,  None, None,  units.amu,
            units.elementary_charge,  None,
            units.kilojoules_per_mole * units.nanometers**(6),
            units.kilojoules_per_mole * units.nanometers**(12))
    def __init__(self,  atomtype,  bondtype, atomic_number,  mass,  charge,
            ptype,  C6,  C12):
        AbstractAtomType.__init__(self,  atomtype,  atomic_number,  mass,
               charge,  ptype)
        self.C6 = C6
        self.C12 = C12

