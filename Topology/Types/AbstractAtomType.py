from cctools.Decorators import *

class AbstractAtomType(object):

    def __init__(self, atomtype, Z=-1, m=-1, q=None, ptype=None):
        self.atomtype = atomtype
        self.Z = Z
        self.m = m
        self.q = q
        self.ptype = ptype

    def __eq__(self, object):
        return self.atomtype == object.atomtype

    def __hash__(self):
        return hash(self.atomtype)

