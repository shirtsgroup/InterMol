"""
.. module:: Molecule
   :platform: Unix

.. moduleuthor:: Christoph Klein <ctk3b@virginia.edu>, Christopher
Lee <ctl4f@virginia.edu>
"""

from Atom import *
from OrderedSet import OrderedSet

class Molecule(object):
    """An abstract molecule object.
    """
    def __init__(self, name = None):
        """Initialize the molecule
        
        Args:
            name (str): name of the molecule
        """
        if name != None:
            self.name = name
        else:
            # TODO Fix the naming resolution
            self.name = "Untitled"
        self._atoms = OrderedSet()
        
    def addAtom(self, atom): 
        """Add and atom 
        
        Args:
            atom (atom): the atom to add into the molecule
        """
        self._atoms.add(atom)   

    def removeAtom(self, atom):
        """Remove Atom
        
        Args:
            atom (atom): the atom to remove from the molecule
        """
        self._atoms.remove(atom)

    def getAtoms(self):
        """Return an orderedset of atoms
        """
        return self._atoms    

    def __repr__(self):
        return self.name 
