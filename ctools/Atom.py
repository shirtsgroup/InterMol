"""
... module:: Atom
    :platform: Unix

.. moduleuthor:: Christoph Klein <ctk3b@virginia.edu>, Christopher
Lee <ctl4f@virginia.edu>
"""
import simtk.unit as units
from ctools.Converter import *

class Atom(object):
    def __init__(self, 
                atomIndex,
                atomName = None,
                residueIndex = -1, 
                residueName = None):
        """Create an Atom object
        
        Args:
            atomIndex (int): index of atom in the molecule
            atomName (str): name of the atom (eg., N, C, H, O) 
            residueIndex (int): index of residue in the molecule
            residueName (str): name of the residue (eg., THR, CYS)
        """
        self._atomIndex = atomIndex
        self._atomName = atomName
        self._residueIndex = residueIndex
        self._residueName = residueName
        self._position = [0 * units.nanometers,
                            0 * units.nanometers,
                            0 * units.nanometers]
        self._velocity = [0 * units.nanometers * units.picosecond**(-1), 
                            0 * units.nanometers * units.picosecond**(-1), 
                            0 * units.nanometers * units.picosecond**(-1)]
        self._force = [0 * units.kilojoules_per_mole * units.nanometers**(-1),
                        0 * units.kilojoules_per_mole * units.nanometers**(-1),
                        0 * units.kilojoules_per_mole * units.nanometers**(-1)]
        
        # These will be added after the ctools is read in and come from [ atomtypes ]
        self._atomtype = list()
        self._bondtype = None
        self._Z = None
        self._cgnr = None
        self._mass = list()
        self._charge = list()
        self._ptype = "A"
        self._sigma = list() 
        self._epsilon = list()

    def setPosition(self, x, y, z):
        """Sets the position of the atom
        
        Args:
            x (float): x position
            y (float): y position
            z (float): z position
        """
        unit = units.nanometers
        x = convert_units(x, unit)
        y = convert_units(y, unit)
        z = convert_units(z, unit)
        self._position = [x, y, z]

    def getPosition(self):
        """Gets the position fo the atom
        
        Returns:
            Tuple [x, y, z]
        """
        return self.position
    
    def setVelocity(self, vx, vy, vz):
        """Sets the velocity of the atom
        
        Args:
            vx (float): x velocity
            vy (float): y velocity
            vz (float): z velocity
        """
        unit = units.nanometers * units.picoseconds**(-1)
        vx = convert_units(vx, unit)
        vy = convert_units(vy, unit)
        vz = convert_units(vz, unit) 
        self._velocity = [vx, vy, vz]
    
    def getVelocity(self):
        """Gets the velocity of the atom
        
        Returns:
            Tuple [vx, vy, vz]
        """
        return self._velocity
    
    def setForce(self, fx, fy, fz):
        """Sets the force of the atom
        
        Args:
            fx (float): x force
            fy (float): y force
            fz (float): z force 
        """
        unit = units.kilojoules_per_mole * units.nanometers**(-1)
        fx = convert_units(fx, unit)
        fy = convert_units(fy, unit)
        fz = convert_units(fz, unit)
        self._force = [fx, fy, fz]    
        
    def getForce(self):
        """Gets the force of the atom
            
        Returns:
            Tuple [fx, fy, fz]
        """
        return self._force    

    def setMass(self, mass):
        """Sets the mass of the atom
        
        Args:
            mass (float): mass of the atom
        """
        unit = units.amu
        self._mass = convert_units(mass, unit)

    def getMass(self):
        """Gets the mass of the atom
            
        Returns:
            mass (float): mass of the atom
        """
        return self._mass

    def setCharge(self, charge):
        """Sets the charge of the atom
        
        Args:
            charge (float): Charge of the atom
        """
        unit = units.elementary_charge
        self._charge = convert_units(charge)

    def getCharge(self):
        """Gets the charge of the atom
        
        Returns:
            charge (float): Charge of the atom
        """
        return self._charge 
    
    def __repr__(self):
        """
        """
        return ('Atom(' + str(self._atomIndex) + ", " +  str(self._atomName) + ")")

    def __cmp__(self, other):
        """
        """
        return self._atomIndex - other._atomIndex
    
    def __eq__(self, other):
        """
        """
        return self._atomIndex == other._atomIndex
    
    def __hash__(self):
        return hash(self._atomIndex)
