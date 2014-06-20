from intermol.decorators import *
from abstract_dihedral import *

class FourierDihedral(AbstractDihedral):

    @accepts_compatible_units(None, None, None, None,
        units.kilojoules_per_mole, units.kilojoules_per_mole,
        units.kilojoules_per_mole, units.kilojoules_per_mole)
    def __init__(self, atom1, atom2, atom3, atom4, c1, c2, c3, c4):
        """
        """
        AbstractDihedral.__init__(self, atom1, atom2, atom3, atom4)
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4

    def get_parameters(self):
        return (self.atom1, self.atom2, self.atom3, self.atom4, self.phi, self.c1, self.c2, self.c3, self.c4)  
    
    def __repr__(self):
        return str(self.atom1)+'  '+str(self.atom2)+'  '+ str(self.atom3)+'  '+str(self.atom4)+'  '+str(self.c1)+'  '+str(self.c2)+'  '+str(self.c3)+'  '+str(self.c4)

    def __str__(self):
        return str(self.atom1)+'  '+str(self.atom2)+'  '+ str(self.atom3)+'  '+str(self.atom4)+'  '+str(self.c1)+'  '+str(self.c2)+'  '+str(self.c3)+'  '+str(self.c4)
        


