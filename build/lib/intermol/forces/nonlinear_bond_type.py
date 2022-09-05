import parmed.unit as units

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


class NonlinearBond(NonlinearBondType):
    """
    http://lammps.sandia.gov/doc/bond_nonlinear.html
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 epsilon=0.0 * units.kilojoules_per_mole,
                 r0=0.0 * units.nanometers,
                 lamda=0.0 * units.nanometers,
                 order=1, c=False):
        self.atom1 = atom1
        self.atom2 = atom2
        NonlinearBondType.__init__(self, bondingtype1, bondingtype2, 
                epsilon=epsilon,
                r0=r0,
                lamda=lamda,
                order=order, c=c)