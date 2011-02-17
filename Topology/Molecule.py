import sys
from cctools.Atom import *
from cctools.OrderedSet import OrderedSet

class Molecule(object):
    def __init__(self, name = None):
        """
        Create a new Structure object.
        """
        if name != None:
            self.name = name
        else:
            # TODO Fix the naming resolution
            self.name = "Untitled"
        self.atoms = OrderedSet()
        
    def setNrexcl(self, nrexcl):
        self.nrexcl = nrexcl
    
    def getNrexcl(self):
        return self.nrexcl

    def addAtom(self, atom): 
        self.atoms.add(atom)   

    def removeAtom(self, atom):
        self.atoms.remove(atom)

    def getAtoms(self):
        return self.atoms    

    def __repr__(self):
        return self.name 
