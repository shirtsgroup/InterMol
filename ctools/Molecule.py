"""
.. module:: Molecule
   :platform: Unix

.. moduleuthor:: Christoph Klein <ctk3b@virginia.edu>, Christopher
Lee <ctl4f@virginia.edu>
"""

from ctools.Atom import *
from ctools.OrderedSet import OrderedSet

class Molecule(object):
    """An abstract molecule object.
    """
    def __init__(self, name = None):
        """Initialize the molecule
        """
        if name != None:
            self.name = name
        else:
            # TODO Fix the naming resolution
            self.name = "Untitled"
        self._atoms = OrderedSet()
        
    def addAtom(self, atom): 
        self._atoms.add(atom)   

    def removeAtom(self, atom):
        self._atoms.remove(atom)

    def getAtoms(self):
        return self._atoms    

    def __repr__(self):
        return self.name 
