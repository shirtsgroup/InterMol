from intermol.decorators import *

class DistanceRestraint(object): 

    @accepts_compatible_units(None, None, None, None, units.nanometers,  units.nanometers, units.nanometers, None)
    def __init__(self, atom1, atom2, type, label, low, up1, up2, weight):
        """
        """
        if atom1:
            self.atom1 = atom1
        if atom2:
            self.atom2 = atom2
        self.type = type
        self.label = label
        self.low = low
        self.up1 = up1
        self.up2 = up2
        self.weight = weight

    def get_parameters(self):
        return (self.atom1, self.atom2, self.type, self.label, self.low, self.up1, self.up2, self.weight)

