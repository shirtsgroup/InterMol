from Topology.OrderedSet import OrderedSet
from Topology.HashMap import HashMap
class MoleculeType(object):
    def __init__(self, name):
        self.name = name 
        self.moleculeSet = OrderedSet()
        self.bondForceSet = HashMap()
        self.pairForceSet = HashMap()
        self.angleForceSet = HashMap()
        self.dihedralForceSet = HashMap()
        self.constraints = HashMap()
        self.exclusions = HashMap()
        self.settles = None
        self.nrexcl = None
    
    def addMolecule(self, molecule):
        self.moleculeSet.add(molecule)
    
    def removeMolecule(self, molecule):
        self.moleculeSet.remove(molecule)
    
    def getMolecule(self, molecule):
        return get_equivalent(self.moleculeSet, molecule, False)

    def addForce(self, force):
        self.forceSet.add(force)

    def removeForce(self, force):
        self.forceSet.remove(force)    
    
    def getForce(self, force):
        return get_equivalent(self.forceSet, force, False)

    def setNrexcl(self, nrexcl):
        self.nrexcl = nrexcl

    def getNrexcl(self):
        return self.nrexcl
