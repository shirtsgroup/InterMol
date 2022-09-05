import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_type import AbstractType


class RigidWater(AbstractType):
    @accepts_compatible_units(None, None, None, units.nanometers, units.nanometers)
    def __init__(self, atom1, atom2, atom3, dOH, dHH):
        """
        """
        super(RigidWater, self).__init__()
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.dOH = dOH
        self.dHH = dHH
