import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_bond_type import AbstractBondType


class NonlinearBondType(AbstractBondType):
    __slots__ = ['epsilon', 'r0', 'lamda', 'order', 'c']

    @accepts_compatible_units(None, None, 
                              epsilon=units.kilojoules_per_mole,
                              r0=units.nanometers,
                              lamda=units.nanometers,
                              order=None,
                              c=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 epsilon=0.0 * units.kilojoules_per_mole,
                 r0=0.0 * units.nanometers,
                 lamda=0.0 * units.nanometers,
                 order=1, c=False):
        AbstractBondType.__init__(self, bondingtype1, bondingtype2, order, c)
        self.epsilon = epsilon
        self.r0 = r0
        self.lamda = lamda

