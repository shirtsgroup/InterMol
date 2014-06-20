from intermol.decorators import *

class DihedralRestraint(object):

    @accepts_compatible_units(None, 
            None, 
            None, 
            None, 
            None, 
            None, 
            units.degrees, 
            units.degrees, 
            None)
    def __init__(self, atom1, atom2, atom3, atom4, type, label, phi, delphi, weight):
        """
        """
        if atom1:
            self.atom1 = atom1
        if atom2:
            self.atom2 = atom2
        if atom3:
            self.atom3 = atom3
        if atom4:
            self.atom4 = atom4
        self.type = type
        self.label = label
        self.phi = phi
        self.delphi = delphi
        self.weight = weight
    
    def get_parameters(self):
        return (self.atom1, self.atom2, self.atom3, self.atom4, self.type, self.label, self.phi, self.delphi, self.weight)

