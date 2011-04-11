from Topology.Decorators import *

class Atom(object):
    @accepts_compatible_units(None, 
                None, 
                None,
                None,
                units.nanometers, 
                units.nanometers, 
                units.nanometers, 
                units.nanometers * units.picoseconds**(-1), 
                units.nanometers * units.picoseconds**(-1), 
                units.nanometers * units.picoseconds**(-1),
                units.kilojoules_per_mole * units.nanometers**(-1), 
                units.kilojoules_per_mole * units.nanometers**(-1), 
                units.kilojoules_per_mole * units.nanometers**(-1))

    def __init__(self, 
                atomNum,
                resNum = -1, 
                resName = None,
                atomName = None,
                x = None,
                y = None,
                z = None,
                vx = 0.0,
                vy = 0.0,
                vz = 0.0,
                fx = 0.0,
                fy = 0.0,
                fz = 0.0):
        self.atomNum = atomNum
        self.resNum = resNum
        self.resName = resName
        self.atomName = atomName
        self.position = [x, y, z]
        self.velocity = [vx, vy, vz]
        self.force = [fx, fy, fz]
        
        # These will be added after the Topology is read in and come from [ atomtypes ]
        self.atomtype = list()
        self.bondtype = None
        self.Z = None
        self.cgnr = None
        self.mass = list()
        self.charge = list()
        self.ptype = "A"
        self.sigma = list() 
        self.epsilon = list()

    @accepts_compatible_units(
        units.nanometers,
        units.nanometers,
        units.nanometers) 
    def setPosition(self, x, y , z):
        self.position = [x, y, z]

    def getPosition(self):
        return self.position
    
    @accepts_compatible_units(
        units.nanometers * units.picoseconds**(-1),
        units.nanometers * units.picoseconds**(-1),
        units.nanometers * units.picoseconds**(-1))
    def setVelocity(self, vx, vy, vz):
        self.velocity = [vx, vy, vz]
    
    def getVelocity(self):
        return self.velocity
    
    @accepts_compatible_units(
        units.kilojoules_per_mole * units.nanometers**(-1),
        units.kilojoules_per_mole * units.nanometers**(-1),
        units.kilojoules_per_mole * units.nanometers**(-1))
    def setForce(self, fx, fy, fz):
        self.force = [fx, fy, fz]    
        
    def getForce(self):
        return self.force

    def setId(self, id):
        self.id = id

    def getId(self):
        return self.id

    def setZ(self, Z):
        self.Z = Z

    def getZ(self):
        return self.Z

    @accepts_compatible_units(units.amu)
    def setMass(self, mass):
        self.mass = mass

    def getMass(self):
        return self.mass

    @accepts_compatible_units(units.elementary_charge)
    def setCharge(self, charge):
        self.charge = charge

    def getCharge(self):
        return self.charge 
    
    def setPtype(self, ptype):
        self.ptype = ptype

    def getPtype(self):
        return self.ptype

    def setSigma(self, sigma):
        self.sigma = sigma
    
    def getSigma(self):
        return self.sigma
    
    def setEpsilon(self, epsilon):
        self.epsilon = epsilon        
 
    def getEpsilon(self):
        return self.epsilon
    
    def addPotential(self, potential):
        self.potentials.add(potential)

    def getPotentials(self):
        return self.potentials

    def removePotential(self, potential):
        self.potentials.remove(potential)

    def __repr__(self):
        # TODO fill out this repr
        return str(self.atomNum) + " " +  self.atomName

    def __cmp__(self, other):
        return self.atomNum - other.atomNum
    
    def __eq__(self, other):
        return self.atomNum == other.atomNum
    #TODO need to update hash
    def __hash__(self):
        return hash(self.atomNum)
