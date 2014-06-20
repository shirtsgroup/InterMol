from intermol.decorators import *
from abstract_dihedral import *

class ProperDihedral9(AbstractDihedral):

    @accepts_compatible_units(None, None, None, None, units.degrees, units.kilojoules_per_mole, None)
    def __init__(self, atom1, atom2, atom3, atom4, phi, k, multiplicity):
        """
        """
        AbstractDihedral.__init__(self, atom1, atom2, atom3, atom4, multiplicity)

        self.phi = phi
        self.k = k

    def get_parameters(self):
        return (self.atom1, self.atom2, self.atom3, self.atom4, self.phi, self.k, self.multiplicity)

    def __repr__(self):
        return str(self.atom1)+'  '+str(self.atom2)+'  '+ str(self.atom3)+'  '+str(self.atom4)+'  '+str(self.phi)+'  '+str(self.k)+'  '+str(self.multiplicity)

    def __str__(self):
        return str(self.atom1)+'  '+str(self.atom2)+'  '+ str(self.atom3)+'  '+str(self.atom4)+'  '+str(self.phi)+'  '+str(self.k)+'  '+str(self.multiplicity)




