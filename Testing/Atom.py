"""
... module:: Atom
    :platform: Unix

.. moduleuthor:: Christoph Klein <ctk3b@virginia.edu>, Christopher
Lee <ctl4f@virginia.edu>
"""
import simtk.unit as units
from Converter import *

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
        self.atomIndex = atomIndex
        self.atomName = atomName
        self.residueIndex = residueIndex
        self.residueName = residueName
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
        self._atomtype = dict()
        self.bondtype = None
        self.Z = None
        self.cgnr = None
        self._mass = dict()
        self._charge = dict()
        self.ptype = "A"
        self._sigma = dict() 
        self._epsilon = dict()
    
    def setSigma(self, index, sigma):
        """Sets the sigma
        
        Args:
            sigma (float): sigma of the atom
            index (int): index to insert at
        """
        self._sigma[index] = sigma   

    def getSigma(self, index = None):
        """
        """
        if index:
            return self._sigma[index]
        return self._sigma 

    def setEpsilon(self, index, epsilon):
        """Sets the epsilon
    
        Args:
            epsilon (float): epsilon of the atom
            index(int): index corresponding to epsilon
        """
        self._epsilon[index] = epsilon

    def getEpsilon(self, index = None):
        """
        """
        if index:
            return self._epsilon[index]
        return self._epsilon    

    def setCgnr(self, index, cgnr):
        """Sets the Cgnr
        
        Args:
            cgnr (int): The charge group number
            index (int): the value corresponding with cgnr precedence

        """
        self._cgnr[index] = cgnr

    def getCgnr(self, index = None):
        """Gets the Cgnr
        
        Args:
            index (int): the index to retrieve, defaults to None
        
        Returns:
            cngr (dict, int): returns the index or the dictionary depending on if index is set
        """
        if index:
            return self._cgnr[index]
        return self._cgnr 

    def setAtomType(self, index, atomtype):
        """Sets the atomtype

        Args:
            atomtype (str): the atomtype of the atom
            index (str): the value corresponding with type precedence (A Type, B Type)
        """
        self._atomtype[index] = atomtype
    
    def getAtomType(self, index = None):
        """Gets the atomtype
        
        Args:
            index (str): the value corresponding with type precedence (A Type, B Type)
        
        Returns:
            atomtype (list, str): Returns the atomtype list or the value at index if index is specified
        """
        if index:
            return self._atomtype[index]
        return self._atomtype
    
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

    def setMass(self, index, mass):
        """Sets the mass of the atom
        
        Args:
            mass (float): mass of the atom
            index (str): the index corresponding with mass precedence (A Mass, B Mass)
        """
        unit = units.amu
        self._mass[index] = convert_units(mass, unit)

    def getMass(self, index = None):
        """Gets the mass of the atom
            
        Returns:
            mass (float): mass of the atom
            index (str): index to retrieve
        """
        if index:
            return self._mass[index]
        return self._mass

    def setCharge(self, index, charge):
        """Sets the charge of the atom
        
        Args:
            charge (float): Charge of the atom
            index (int): the index corresponding with charge precedence 
        """
        unit = units.elementary_charge
        self._charge[index] = convert_units(charge, unit)

    def getCharge(self, index = None):
        """Gets the charge of the atom
        
        Args:
            index (int): index of the charge to retrieve defaults to None   
     
        Returns:
            charge (float): Charge of the atom
        """
        if index:
            return self._charge[index]
        return self._charge 
    
    def __repr__(self):
        """
        """
        return ('Atom(' + str(self.atomIndex) + ", " +  str(self.atomName) + ")")

    def __cmp__(self, other):
        """
        """
        return self.atomIndex - other.atomIndex
    
    def __eq__(self, other):
        """
        """
        return self.atomIndex == other.atomIndex
    
    def __hash__(self):
        """
        """
        return hash(self.atomIndex)
