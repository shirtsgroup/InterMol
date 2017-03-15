import parmed.unit as units

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


class LjqCPair(LjqCPairType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 qi=0.0 * units.elementary_charge,
                 qj=0.0 * units.elementary_charge,
                 C6=0.0 * units.kilojoules_per_mole * units.nanometers ** (6),
                 C12=0.0 * units.kilojoules_per_mole * units.nanometers ** (12),
                 scaleLJ=None, scaleQQ=None, long=False):
        self.atom1 = atom1
        self.atom2 = atom2
        LjqCPairType.__init__(self, bondingtype1, bondingtype2, 
                qi=qi,
                qj=qj,
                C6=C6,
                C12=C12,
                scaleLJ=scaleLJ, scaleQQ=scaleQQ, long=long)