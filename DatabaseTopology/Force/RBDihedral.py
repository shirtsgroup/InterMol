from Topology.Decorators import *
from Topology.Force.AbstractDihedral import *


class RBDihedral(AbstractDihedral):

    @accepts_compatible_units(None, None, None, None, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole)
    def __init__(self, atom1, atom2, atom3, atom4, C0, C1, C2, C3, C4, C5):
        """
        """
        AbstractDihedral.__init__(self, atom1, atom2, atom3, atom4, -1)
        self.C0 = C0
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4
        self.C5 = C5

    def getForceParameters(self):
        return (self.atom1, self.atom2, self.atom3, self.atom4, self.C0, self.C1, self.C2, self.C3, self.C4, self.C5)

    def __repr__(self):
        print self.atom1+'  '+self.atom2+'  '+ self.atom3+'  '+self.atom4+'  '+self.C0+'  '+self.C1+'  '+self.C2+'  '+self.C3+'  '+self.C4+'  '+self.C5


    def __str__(self):
        print self.atom1+'  '+self.atom2+'  '+ self.atom3+'  '+self.atom4+'  '+self.C0+'  '+self.C1+'  '+self.C2+'  '+self.C3+'  '+self.C4+'  '+self.C5
