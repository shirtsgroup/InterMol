import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_pair_type import AbstractPairType


class LjqCPairType(AbstractPairType):
    __slots__ = ['qi', 'qj', 'C6', 'C12', 'scaleLJ', 'scaleQQ', 'long']

    @accepts_compatible_units(None, None, 
                              qi=units.elementary_charge,
                              qj=units.elementary_charge,
                              C6=units.kilojoules_per_mole * units.nanometers ** (6),
                              C12=units.kilojoules_per_mole * units.nanometers ** (12),
                              scaleLJ=None,
                              scaleQQ=None,
                              long=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 qi=0.0 * units.elementary_charge,
                 qj=0.0 * units.elementary_charge,
                 C6=0.0 * units.kilojoules_per_mole * units.nanometers ** (6),
                 C12=0.0 * units.kilojoules_per_mole * units.nanometers ** (12),
                 scaleLJ=None, scaleQQ=None, long=False):
        AbstractPairType.__init__(self, bondingtype1, bondingtype2, scaleLJ, scaleQQ, long)
        self.qi = qi
        self.qj = qj
        self.C6 = C6
        self.C12 = C12
