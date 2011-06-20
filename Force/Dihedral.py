import warnings
import simtk.unit as units
from Topology.Converter import *

class Dihedral(object):
    def __init__(self, atom1, atom2, atom3, atom4, func):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.func = func
        self._phi = units.Quantity(0)
        self._k = units.Quantity(0)
        self.multiplicity = 0
        self._c1 = units.Quantity(0)
        self._c2 = units.Quantity(0)
        self._c3 = units.Quantity(0)
        self._c4 = units.Quantity(0)
        self._c5 = units.Quantity(0)
        self._c6 = units.Quantity(0)
        
    def setK(self, k):
        if self.func == 1 or self.func == 9 or self.func == 4 or self.func == 8:
            unit = units.kilojoules_per_mole
        elif self.func == 2:
            unit = units.kilojoules_per_mole * units.radians**(-2)
        else:
            warnings.warn("Dihedral with function number %s does not have value k" %(self.func), UserWarning)
            return
        self._k = convert_units(k, unit)

    def setPhi(self, phi):
        if self.func == 1 or self.func == 9 or self.func == 2 or self.func == 4:
            unit = units.degrees
        else:
            warnings.warn("Dihedral with function number %s does not have value phi" %(self.func), UserWarning)
            return
        self._phi = convert_units(phi, unit)

    def setC1(self, c1):
        if self.func == 3 or self.func == 5:
            unit = units.kilojoules_per_mole
        else:
            warnings.warn("Dihedral with function number %s does not have value c1" %(self.func), UserWarning)
            return
        self._c1 = convert_units(c1, unit)
 
    def setC2(self, c2):
        if self.func == 3 or self.func == 5:
            unit = units.kilojoules_per_mole
        else:
            warnings.warn("Dihedral with function number %s does not have value phi" %(self.func), UserWarning)
            return
        self._c2 = convert_units(c2, unit)

    def setC3(self, c3):
        if self.func == 3 or self.func == 5:
            unit = units.kilojoules_per_mole
        else:
            warnings.warn("Dihedral with function number %s does not have value phi" %(self.func), UserWarning)
            return
        self._c3 = convert_units(c3, unit)

    def setC4(self, c4):
        if self.func == 3 or self.func == 5:
            unit = units.kilojoules_per_mole
        else:
            warnings.warn("Dihedral with function number %s does not have value phi" %(self.func), UserWarning)
            return
        self._c4 = convert_units(c4, unit)

    def setC5(self, c5):
        if self.func ==3:
            unit = units.kilojoules_per_mole
        else:
            warnings.warn("Dihedral with function number %s does not have value phi" %(self.func), UserWarning)
            return
        self._c5 = convert_units(c5, unit)

    def setC6(self, c6):
        if self.func == 3:
            unit = units.kilojoules_per_mole
        else:
            warnings.warn("Dihedral with function number %s does not have value phi" %(self.func), UserWarning)
            return
        self._c6 = convert_units(c6, unit)
           
