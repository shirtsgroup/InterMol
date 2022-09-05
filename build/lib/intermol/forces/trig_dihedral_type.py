import parmed.unit as units

from intermol.decorators import accepts_compatible_units
from intermol.forces.abstract_dihedral_type import AbstractDihedralType


class TrigDihedralType(AbstractDihedralType):
    __slots__ = ['phi', 'fc0', 'fc1', 'fc2', 'fc3', 'fc4', 'fc5', 'fc6', 'improper']

    @accepts_compatible_units(None, None, None, None, 
                              phi=units.degrees,
                              fc0=units.kilojoules_per_mole,
                              fc1=units.kilojoules_per_mole,
                              fc2=units.kilojoules_per_mole,
                              fc3=units.kilojoules_per_mole,
                              fc4=units.kilojoules_per_mole,
                              fc5=units.kilojoules_per_mole,
                              fc6=units.kilojoules_per_mole,
                              improper=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                 phi=0.0 * units.degrees,
                 fc0=0.0 * units.kilojoules_per_mole,
                 fc1=0.0 * units.kilojoules_per_mole,
                 fc2=0.0 * units.kilojoules_per_mole,
                 fc3=0.0 * units.kilojoules_per_mole,
                 fc4=0.0 * units.kilojoules_per_mole,
                 fc5=0.0 * units.kilojoules_per_mole,
                 fc6=0.0 * units.kilojoules_per_mole,
                 improper=False):
        AbstractDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, improper)
        self.phi = phi
        self.fc0 = fc0
        self.fc1 = fc1
        self.fc2 = fc2
        self.fc3 = fc3
        self.fc4 = fc4
        self.fc5 = fc5
        self.fc6 = fc6


class TrigDihedral(TrigDihedralType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, atom4, bondingtype1=None, bondingtype2=None, bondingtype3=None, bondingtype4=None, 
                 phi=0.0 * units.degrees,
                 fc0=0.0 * units.kilojoules_per_mole,
                 fc1=0.0 * units.kilojoules_per_mole,
                 fc2=0.0 * units.kilojoules_per_mole,
                 fc3=0.0 * units.kilojoules_per_mole,
                 fc4=0.0 * units.kilojoules_per_mole,
                 fc5=0.0 * units.kilojoules_per_mole,
                 fc6=0.0 * units.kilojoules_per_mole,
                 improper=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        TrigDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                phi=phi,
                fc0=fc0,
                fc1=fc1,
                fc2=fc2,
                fc3=fc3,
                fc4=fc4,
                fc5=fc5,
                fc6=fc6,
                improper=improper)