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
    def __init__(self, molType, molIndex):
        self.molType = molType
        self.molIndex = molIndex 
    
    def __repr__(self):
        return self.molType 
