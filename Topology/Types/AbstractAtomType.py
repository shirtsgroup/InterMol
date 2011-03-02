from Topology.Decorators import *

class AbstractAtomType(object):

    def __init__(self, 
            atomtype, 
            Z =  None, 
            mass = None, 
            charge = None, 
            ptype = None):
        self.atomtype = atomtype
        self.Z = Z
        self.mass = mass
        self.charge = charge
        self.ptype = ptype

    def __eq__(self, object):
        return self.atomtype == object.atomtype

    def __hash__(self):
        return hash(self.atomtype)

