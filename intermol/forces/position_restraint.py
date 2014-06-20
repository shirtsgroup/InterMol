from intermol.decorators import *

class PositionRestraint(object):

    @accepts_compatible_units(None, units.kilojoules_per_mole * units.nanometers**(-2),  units.kilojoules_per_mole * units.nanometers**(-2), units.kilojoules_per_mole * units.nanometers**(-2))
    def __init__(self, atom1, kx, ky, kz):
        """
        """
        if atom1:
            self.atom1 = atom1
        self.kx = kx
        self.ky = ky
        self.kz = kz

    def get_parameters(self):
        return (self.atom1, self.kx, self.ky, self.kz)
        
