from Topology.OrderedSet import OrderedSet
from Topology.HashMap import HashMap
class MoleculeType(object):
    def __init__(self, name, nrexcl, settles):
        self.name = name 
        self.settles = None
        self.nrexcl = None
