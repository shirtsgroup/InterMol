import simtk.unit as units
from Topology.Converter import *
import warnings

class Angle(object):
    def __init__(self, atom1, atom2, atom3, func):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.func = func
        self._theta = units.Quantity(0)
        self._k = units.Quantity(0)
        self._k2 = units.Quantity(0)
        self._c1 = units.Quantity(0)
        self._c2 = units.Quantity(0)
        self._c3 = units.Quantity(0)
        self._c4 = units.Quantity(0)
        self._c5 = units.Quantity(0)
            
    def setTheta(self, theta):
        unit = units.degrees
        self._theta = convert_units(theta, unit) 
 
    def setK(self, k):
        if self.func == 1:
            unit = units.kilojoules_per_mole * units.radians**(-2)
            self._k = convert_units(k, unit)
        elif self.func == 2 or self.func == 5 or self.func == 8:
            unit = units.kilojoules_per_mole
            self._k = convert_units(k, unit)
        elif self.func == 3 or self.func == 4:
            unit = units.kilojoules_per_mole * unit.nanometers**(-2)
            self._k = convert_units(k, unit)
        else:
            warnings.warn("Angle with function number %s does not have value k" %(self.func), UserWarning)
    
    def setK2(self, k2):
        if self.func == 5:
            unit = units.kilojoules_per_mole
            self._k2 = convert_units(k2, unit)
        else:
            warnings.warn("Angle with function number %s does not have value k2" %(self.func), UserWarning)

    def setC1(self, c1):
        if self.func == 3 or self.func == 4 or self.func == 5:
            unit = units.nanometers
            self._c1 = convert_units(c1, unit)
        elif self.func == 6:
            unit = units.kilojoules_per_mole
            self._c1 = convert_units(c1, unit)
        else:
            warnings.warn("Angle with function number %s does not have value c1" %(self.func), UserWarning)

    
    def setC2(self, c2):
        if self.func == 3 or self.func == 4:
            unit = units.nanometers
            self._c2 = convert_units(c2, unit)
        elif self.func == 6:
            unit = units.kilojoules_per_mole * units.radians**(-1)
            self._c2 = convert_units(c2, unit)
        else:
            warnings.warn("Angle with function number %s does not have value c2" %(self.func), UserWarning)

    def setC3(self, c3):
        if self.func == 4:
            unit = units.nanometers
            self._c3 = convert_units(c3, unit)
        elif self.func == 6:
            unit = units.kilojoules_per_mole * units.radians**(-2)
            self._c3 = convert_units(c3, unit)
        else:
            warnings.warn("Angle with function number %s does not have value c3" %(self.func), UserWarning) 

    def setC4(self, c4):
        if self.func == 6:
            unit = units.kilojoules_per_mole * units.radians**(-3)
            self._c4 = covert_units(c4, unit)
        else:
            warnings.warn("Angle with function number %s does not have value c4" %(self.func), UserWarning)
    
    def setC5(self, c5):
        if self.func == 6:
            unit = units.kilojoules_per_mole * units.raidans**(-4)
            self._c5 = convert_units(c5, unit)
        else:
            warnings.warn("Angle with function number %s does not have value c5" %(self.func), UserWarning)

    def __repr__(self):
        return str(self.atom1)+'  '+str(self.atom2)+'  '+ str(self.atom3)+'  '+ str(self.theta)+'  '+str(self.k)


    def __str__(self):
        return str(self.atom1)+'  '+str(self.atom2)+'  '+ str(self.atom3)+'  '+ str(self.theta)+'  '+str(self.k)

