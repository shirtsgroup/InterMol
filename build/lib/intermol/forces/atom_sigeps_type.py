import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_atom_type import AbstractAtomType


class AtomSigepsType(AbstractAtomType):
    @accepts_compatible_units(None, None,  None, units.amu,
            units.elementary_charge,  None,  units.nanometers,
            units.kilojoules_per_mole)
    def __init__(self, atomtype, bondtype, atomic_number, mass, charge, ptype,
                 sigma, epsilon):
        AbstractAtomType.__init__(self, atomtype, bondtype, atomic_number,
                mass, charge, ptype)
        self.sigma = sigma
        self.epsilon = epsilon


