"""
.. module:: Molecule
    :platform: Unix

.. moduleuthor:: Christoph Klein <ctk3b@virginia.edu>, Christopher
Lee <ctl4f@virginia.edu>
"""

import sys
from Topology.Atom import *
from Topology.OrderedSet import OrderedSet

class Molecule(object):
    """The abstract molecule object.
        
        blah blah blah
    """
    def __init__(self, name = None):
        """Create a new Structure object.

        Args:
            
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
