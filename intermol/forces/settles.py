import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from abstract_type import AbstractType


class Settles(AbstractType):
    @accepts_compatible_units(None, units.nanometers, units.nanometers)
    def __init__(self, atom1, dOH, dHH):
        """
        """
        super(Settles, self).__init__()
        if atom1:
            self.atom1 = atom1
        self.dOH = dOH
        self.dHH = dHH
