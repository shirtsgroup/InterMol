import simtk.unit as units
from Topology.Converter import *
import warnings

class Bond(object):
    # TODO implement reordering of atoms
    def __init__(self, atom1, atom2, func):
        self.atom1 = atom1
        self.atom2 = atom2
        self.func = func
        self._k = units.Quantity(0)
        self._length = units.Quantity(0)
        self._alpha = units.Quantity(0)
        self._beta = units.Quantity(0)

    def setK(self, k):
        unit = None
        if self.func == 1 or self.func == 5 or self.func == 6 or self.func == 7:
            unit = units.kilojoules_per_mole * units.nanometers**(-2)
        elif self.func == 2:
            unit = units.kilojoules_per_mole * units.nanometers**(-4)
        else:
            warnings.warn("Bond with function number = %s does not have value k" %(self.func), UserWarning)
            return
        self._k = convert_units(k, unit) 
    
    def setLength(self, length):
        unit = units.nanometers
        self.length = convert_units(length, unit)

    def setAlpha(self, alpha):
        unit = None
        if self.func == 3:
            unit = units.kilojoules_per_mole
        elif self.func == 4:
            unit = units.kilojoules_per_mole * units.nanometers**(-2)
        else:
            warnings.warn("Bond with function number = %s does not use alpha" %(self.func), UserWarning)
            return
        self._alpha = convert_units(alpha, unit)

    def setBeta(self, beta):
        unit = None
        if self.func == 3:
            unit = units.nanometers**(-1)
        elif self.func == 4:
            unit = units.kilojoules_per_mole * units.nanometers**(-3)
        else:
            warnings.warn("Bond with function number = %s does not use beta" %(self.func), UserWarning)
            return
        self._beta = convert_units(beta, unit)     

    def __repr__(self):
        return str(self.atom1) +'  '+ str(self.atom2)
    
    def __str__(self):
        return str(self.atom1) +'  '+ str(self.atom2)

