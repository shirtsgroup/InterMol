from System import *


#=============================================================================================
# Force base class
#=============================================================================================

class Force(object):
    """
    Force objects apply forces to the particles in a System, or alter their behavior in other
    ways.  This is an abstract class.  Subclasses define particular forces.

    More specifically, a Force object can do any or all of the following:

    * Add a contribution to the force on each particle
    * Add a contribution to the potential energy of the System
    * Modify the positions and velocities of particles at the start of each time step
    * Define parameters which are stored in the Context and can be modified by the user
    * Change the values of parameters defined by other Force objects at the start of each time step

    """
    pass

#=============================================================================================
# NonbondedForce
#=============================================================================================

class NonbondedForce(Force):
    """
    This class implements nonbonded interactions between particles, including a Coulomb force to represent
    electrostatics and a Lennard-Jones force to represent van der Waals interactions.  It optionally supports
    periodic boundary conditions and cutoffs for long range interactions.  Lennard-Jones interactions are
    calculated with the Lorentz-Bertelot combining rule: it uses the arithmetic mean of the sigmas and the
    geometric mean of the epsilons for the two interacting particles.

    To use this class, create a NonbondedForce object, then call addParticle() once for each particle in the
    System to define its parameters.  The number of particles for which you define nonbonded parameters must
    be exactly equal to the number of particles in the System, or else an exception will be thrown when you
    try to create a Context.  After a particle has been added, you can modify its force field parameters
    by calling setParticleParameters().

    NonbondedForce also lets you specify "exceptions", particular pairs of particles whose interactions should be
    computed based on different parameters than those defined for the individual particles.  This can be used to
    completely exclude certain interactions from the force calculation, or to alter how they interact with each other.

    Many molecular force fields omit Coulomb and Lennard-Jones interactions between particles separated by one
    or two bonds, while using modified parameters for those separated by three bonds (known as "1-4 interactions").
    This class provides a convenience method for this case called createExceptionsFromBonds().  You pass to it
    a list of bonds and the scale factors to use for 1-4 interactions.  It identifies all pairs of particles which
    are separated by 1, 2, or 3 bonds, then automatically creates exceptions for them.

    EXAMPLES

    Create a force object.

    >>> force = NonbondedForce()

    Add a particle.

    >>> charge = 1.0 * units.elementary_charge 
    >>> sigma = 1.0 * units.angstrom
    >>> epsilon = 0.001 * units.kilocalories_per_mole
    >>> force.addParticle(charge, sigma, epsilon)
    0

    Set the nonbonded method.

    >>> nonbondedMethod = NonbondedForce.NoCutoff
    >>> force.setNonbondedMethod(nonbondedMethod)

    Return a Swig proxy.

    >>> force_proxy = force.asSwig()

    Create a new force from proxy.

    >>> force_from_proxy = NonbondedForce(force_proxy)

    Create a deep copy.

    >>> import copy
    >>> force_copy = copy.deepcopy(force)

    """

    NoCutoff = 0 #: No cutoff is applied to nonbonded interactions.  The full set of N^2 interactions is computed exactly. This necessarily means that periodic boundary conditions cannot be used.  This is the default.
    CutoffNonPeriodic = 1 #: Interactions beyond the cutoff distance are ignored.  Coulomb interactions closer than the cutoff distance are modified using the reaction field method.
    CutoffPeriodic = 2 #: Periodic boundary conditions are used, so that each particle interacts only with the nearest periodic copy of each other particle.  Interactions beyond the cutoff distance are ignored.  Coulomb interactions closer than the cutoff distance are modified using the reaction field method.
    Ewald = 3 #: Periodic boundary conditions are used, and Ewald summation is used to compute the interaction of each particle with all periodic copies of every other particle.
    PME = 4 #: Periodic boundary conditions are used, and Particle-Mesh Ewald (PME) summation is used to compute the interaction of each particle with all periodic copies of every other particle.

    def __init__(self, force=None):
        """
        Create a NonbondedForce.

        """
        # Initialize with defaults.
        self.particles = list() # particles[i] is the ParticleInfo object for particle i
        self.exceptions = list() # exceptions[i] is the ith ExceptionInfo entry
        self.exceptionMap = dict()
        self.nonbondedMethod = NonbondedForce.NoCutoff # the method of computing nonbonded interactions to use
        self.cutoffDistance = units.Quantity(1.0, units.nanometer) # cutoff used for cutoff-based nonbondedMethod choices, in units of distance
        self.rfDielectric = units.Quantity(78.3, units.dimensionless) # reaction-field dielectric
        self.ewaldErrorTol = 1.0e-4 # relative Ewald error tolerance

        # Populate data structures from swig object, if specified
        if force is not None:
            self._copyDataUsingInterface(self, force)

        return

    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        #TODO: Do we need to call __init__() first, in case class has already been initialized to something else?

        dest.setNonbondedMethod( src.getNonbondedMethod() )
        dest.setCutoffDistance( src.getCutoffDistance() )
        dest.setReactionFieldDielectric( src.getReactionFieldDielectric() )
        dest.setEwaldErrorTolerance( src.getEwaldErrorTolerance() )
        for index in range(src.getNumParticles()):
            args = src.getParticleParameters(index)
            dest.addParticle(*args)
        for index in range(src.getNumExceptions()):
            args = src.getExceptionParameters(index)
            dest.addException(*args)

        return

    def asSwig(self):
        """
        Construct a corresponding Swig object.

        """
        force = openmm.NonbondedForce()
        self._copyDataUsingInterface(force, self)
        return force

    def getNumParticles(self):
        """
        Get the number of particles for which force field parameters have been defined.

        """
        return len(self.particles)

    def getNumExceptions(self):
        """
        Get the number of special interactions that should be calculated differently from other interactions.

        """
        return len(self.exceptions)

    def getNonbondedMethod(self):
        """
        Get the method used for handling long range nonbonded interactions.

        """
        return self.nonbondedMethod

    def setNonbondedMethod(self, nonbondedMethod):
        """
        Set the method used for handling long range nonbonded interactions.

        """
        # TODO: Argument checking.
        self.nonbondedMethod = nonbondedMethod
        return

    def getCutoffDistance(self):
        """
        Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
        is NoCutoff, this value will have no effect.

        """
        return self.cutoffDistance

    @accepts_compatible_units(units.nanometer)
    def setCutoffDistance(self, cutoffDistance):
        """
        Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
        is NoCutoff, this value will have no effect.

        """
        self.cutoffDistance = cutoffDistance
        return

    def getReactionFieldDielectric(self):
        """
        Get the dielectric constant to use for the solvent in the reaction field approximation.

        """
        return self.rfDielectric

    def setReactionFieldDielectric(self, reactionFieldDielectric):
        """
        Set the dielectric constant to use for the solvent in the reaction field approximation.

        @param reactionFieldDielectric     the reaction field dielectric (unitless Quantity)

        """
        self.rfDielectric = reactionFieldDielectric
        return

    def getEwaldErrorTolerance(self):
        """
        Get the error tolerance for Ewald summation.  This corresponds to the fractional error in the forces
        which is acceptable.  This value is used to select the reciprocal space cutoff and separation
        parameter so that the average error level will be less than the tolerance.  There is not a
        rigorous guarantee that all forces on all atoms will be less than the tolerance, however.

        """
        return self.ewaldErrorTol

    @accepts(float)
    def setEwaldErrorTolerance(self, tol):
        """
        Set the error tolerance for Ewald summation.  This corresponds to the fractional error in the forces
        which is acceptable.  This value is used to select the reciprocal space cutoff and separation
        parameter so that the average error level will be less than the tolerance.  There is not a
        rigorous guarantee that all forces on all atoms will be less than the tolerance, however.

        """
        self.ewaldErrorTol = tol
        return 

    @accepts_compatible_units(units.elementary_charge, units.nanometers, units.kilojoules_per_mole)
    def addParticle(self, charge, sigma, epsilon):
        """
        Add the nonbonded force parameters for a particle.  This should be called once for each particle
        in the System.  When it is called for the i'th time, it specifies the parameters for the i'th particle.
        For calculating the Lennard-Jones interaction between two particles, the arithmetic mean of the sigmas
        and the geometric mean of the epsilons for the two interacting particles is used (the Lorentz-Bertelot
        combining rule).

        @param charge    the charge of the particle, measured in units of the proton charge
        @param sigma     the sigma parameter of the Lennard-Jones potential (corresponding to the van der Waals radius of the particle), measured in nm
        @param epsilon   the epsilon parameter of the Lennard-Jones potential (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
        @return the index of the particle that was added

        """
        particle = NonbondedForceParticleInfo(charge, sigma, epsilon)
        self.particles.append(particle)
        return (len(self.particles) - 1)

    def getParticleParameters(self, index):
        """
        Set the nonbonded force parameters for a particle.  When calculating the Lennard-Jones interaction between two particles,
        it uses the arithmetic mean of the sigmas and the geometric mean of the epsilons for the two interacting particles
        (the Lorentz-Bertelot combining rule).

        @param index     the index of the particle for which to set parameters
        @param charge    the charge of the particle, measured in units of the proton charge
        @param sigma     the sigma parameter of the Lennard-Jones potential (corresponding to the van der Waals radius of the particle), measured in nm
        @param epsilon   the epsilon parameter of the Lennard-Jones potential (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol

        """
        particle = self.particles[index]
        return (particle.charge, particle.sigma, particle.epsilon)

    @accepts_compatible_units(None, units.elementary_charge, units.nanometers, units.kilojoules_per_mole)
    def setParticleParameters(self, index, charge, sigma, epsilon):
        """
        Set the nonbonded force parameters for a particle.  When calculating the Lennard-Jones interaction between two particles,
        it uses the arithmetic mean of the sigmas and the geometric mean of the epsilons for the two interacting particles
        (the Lorentz-Bertelot combining rule).

        @param index     the index of the particle for which to set parameters
        @param charge    the charge of the particle, measured in units of the proton charge
        @param sigma     the sigma parameter of the Lennard-Jones potential (corresponding to the van der Waals radius of the particle), measured in nm
        @param epsilon   the epsilon parameter of the Lennard-Jones potential (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol

        """
        particle = NonbondedForceParticleInfo(charge, sigma, epsilon)
        self.particles[index] = particle
        return

    @accepts_compatible_units(None, None, units.elementary_charge**2, units.nanometers, units.kilojoules_per_mole, None)
    def addException(self, particle1, particle2, chargeProd, sigma, epsilon, replace=False):
        """
        Add an interaction to the list of exceptions that should be calculated differently from other interactions.
        If chargeProd and epsilon are both equal to 0, this will cause the interaction to be completely omitted from
        force and energy calculations.

        In many cases, you can use createExceptionsFromBonds() rather than adding each exception explicitly.

        @param particle1  the index of the first particle involved in the interaction
        @param particle2  the index of the second particle involved in the interaction
        @param chargeProd the scaled product of the atomic charges (i.e. the strength of the Coulomb interaction), measured in units of the proton charge squared
        @param sigma      the sigma parameter of the Lennard-Jones potential (corresponding to the van der Waals radius of the particle), measured in nm
        @param epsilon    the epsilon parameter of the Lennard-Jones potential (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
        @param replace    determines the behavior if there is already an exception for the same two particles.  If true, the existing one is replaced.  If false,
                          an exception is thrown.
        @return the index of the exception that was added

        """

        # TODO: Eliminate exceptionMap so that exceptions can be directly manipulated by user in a more pythonic manner without ill consequence?

        if particle1 not in range(self.getNumParticles()):
            raise ValueError("particle1 must be in range(0, getNumParticles())")
        if particle2 not in range(self.getNumParticles()):
            raise ValueError("particle1 must be in range(0, getNumParticles())")

        # Ensure particle1 < particle2.
        if (particle2 < particle1):
            (particle1, particle2) = (particle2, particle1)

        # Create exception entry.
        exception = NonbondedForceExceptionInfo(particle1, particle2, chargeProd, sigma, epsilon)

        # Store it.
        try:
            # Only one set of exceptions are allowed per pair of particles.
            index = self.exceptionMap[(particle1,particle2)]
            if replace:
                self.exceptions[index] = exception
            else:
                raise ValueError("Exception already exists for particle pair (%d,%d); set replace=True optional argument to allow replacement." % (particle1,particle2))
        except:
            # Add the exception and note the pair of particles added.
            index = len(self.exceptions)
            self.exceptions.append(exception)
            self.exceptionMap[(particle1,particle2)] = index

        return (len(self.exceptions) - 1)

    def getExceptionParameters(self, index):
        """
        Get the force field parameters for an interaction that should be calculated differently from others.

        @param index      the index of the interaction for which to get parameters
        @param particle1  the index of the first particle involved in the interaction
        @param particle2  the index of the second particle involved in the interaction
        @param chargeProd the scaled product of the atomic charges (i.e. the strength of the Coulomb interaction), measured in units of the proton charge squared
        @param sigma      the sigma parameter of the Lennard-Jones potential (corresponding to the van der Waals radius of the particle), measured in nm
        @param epsilon    the epsilon parameter of the Lennard-Jones potential (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
        """
        exception = self.exceptions[index]
        return (exception.particle1, exception.particle2, exception.chargeProd, exception.sigma, exception.epsilon)

    def setExceptionParameters(self, index, particle1, particle2, chargeProd, sigma, epsilon):
        """
        Set the force field parameters for an interaction that should be calculated differently from others.
        If chargeProd and epsilon are both equal to 0, this will cause the interaction to be completely omitted from
        force and energy calculations.

        @param index      the index of the interaction for which to get parameters
        @param particle1  the index of the first particle involved in the interaction
        @param particle2  the index of the second particle involved in the interaction
        @param chargeProd the scaled product of the atomic charges (i.e. the strength of the Coulomb interaction), measured in units of the proton charge squared
        @param sigma      the sigma parameter of the Lennard-Jones potential (corresponding to the van der Waals radius of the particle), measured in nm
        @param epsilon    the epsilon parameter of the Lennard-Jones potential (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
        """
        nparticles = self.getNumParticles()
        if particle1 not in range(0,nparticles):
            raise ValueError("particle1 must be in range(0, getNumParticles())")
        if particle2 not in range(0,nparticles):
            raise ValueError("particle1 must be in range(0, getNumParticles())")

        exception = NonbondedForceExceptionInfo(particle1, particle2, chargeProd, sigma, epsilon)
        self.exceptions[index] = exception

        return

    def createExceptionsFromBonds(self, bonds, coulomb14Scale, lj14Scale):
        """
        Identify exceptions based on the molecular topology.  Particles which are separated by one or two bonds are set
        to not interact at all, while pairs of particles separated by three bonds (known as "1-4 interactions") have
        their Coulomb and Lennard-Jones interactions reduced by a fixed factor.

        @param bonds           the set of bonds based on which to construct exceptions.  Each element specifies the indices of
                               two particles that are bonded to each other.
        @param coulomb14Scale  pairs of particles separated by three bonds will have the strength of their Coulomb interaction
                               multiplied by this factor
        @param lj14Scale       pairs of particles separated by three bonds will have the strength of their Lennard-Jones interaction
                               multiplied by this factor
        """

        # Create a Swig object.
        swig_self = self.asSwig()
        # Create exceptions from bonds in Swig object.
        swig_self.createExceptionsFromBonds(bonds, coulomb14Scale, lj14Scale)
        # Reconstitute Python object.
        self._copyDataUsingInterface(self, swig_self)

        return

    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def nparticles(self):
        return len(self.particles)

    @property
    def nexceptions(self):
        return len(self.exceptions)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the System object.

        """

        r = ""

        # Show settings.
        r += "nonbondedMethod: %s" % str(self.nonbondedMethod)
        r += "cutoffDistance: %s" % str(self.cutoffDistance)
        r += "rfDielectric: %s" % str(self.rfDielectric)
        r += "ewaldErrorTolerance: %s" % str(self.ewaldErrorTol)

        # Show particles.
        if (self.nparticles > 0):
            r += "Particles:\n"
            r += "%8s %24s %24s %24s %24s\n" % ("particle", "charge", "sigma", "epsilon", "lambda")
            for (index, particle) in enumerate(self.particles):
                r += "%8d %24s %24s %24s %24f\n" % (index, str(particle.charge), str(particle.sigma), str(particle.epsilon), particle.lambda_)
            r += "\n"

        # Show exceptions
        if (self.nexceptions > 0):
            r += "Exceptions:\n"
            r += "%8s %8s %24s %24s %24s\n" % ("particle1", "particle2", "chargeProd", "sigma", "epsilon")
            for exception in self.exceptions:
                r += "%8d %8d %24s %24s %24s" % (particle1, particle2, str(exception.chargeProd), str(exception.sigma), str(exception.epsilon))
            r += "\n"

        return r

    def _appendForce(self, force, offset):
        """
        Append atoms defined in another force of the same type.

        @param force      the force to be appended
        @param offset     integral offset for atom number

        EXAMPLES
        Create two force objects and append the second to the first.

        >>> force1 = NonbondedForce()
        >>> force2 = NonbondedForce()
        >>> charge = 1.0 * units.elementary_charge
        >>> sigma = 1.0 * units.angstrom
        >>> epsilon = 0.001 * units.kilocalories_per_mole
        >>> force1.addParticle(charge, sigma, epsilon)
        0
        >>> force2.addParticle(charge, sigma, epsilon)
        0
        >>> offset = 1
        >>> force1._appendForce(force2, offset)

        """
        # Check to make sure both forces are compatible.
        if type(self) != type(force):
            raise ValueError("other force object must be of identical Force subclass")

        # Check to make sure both forces have same settings.
        # TODO: This checkins is probably more paranoid than we need.
        if (self.nonbondedMethod != force.nonbondedMethod):
            raise ValueError("other force has incompatible nonbondedMethod")
        if (self.cutoffDistance != force.cutoffDistance):
            raise ValueError("other force has incompatible cutoffDistance")
        if (self.rfDielectric != force.rfDielectric):
            raise ValueError("other force has incompatible rfDielectric")
        if (self.ewaldErrorTol != force.ewaldErrorTol):
            raise ValueError("other force has incompatible ewaldErrorTol")

        # Combine systems.
        for particle in force.particles:
            self.particles.append(particle)
        for exception in force.exceptions:
            exception.particle1 += offset
            exception.particle2 += offset
            self.exceptions.append(exception)

        return

#==========================================================================
# CONTAINER CLASSES
#==========================================================================

class NonbondedForceParticleInfo(object):
    @accepts_compatible_units(units.elementary_charge, units.nanometers, units.kilojoules_per_mole, None)
    def __init__(self, charge, sigma, epsilon, lambda_ = 1.0):
        """
        """
        self.charge = charge
        self.sigma = sigma
        self.epsilon = epsilon
        self.lambda_ = lambda_

        return

class NonbondedForceExceptionInfo(object):
    @accepts_compatible_units(None, None, units.elementary_charge**2, units.nanometers, units.kilojoules_per_mole)
    def __init__(self, particle1, particle2, chargeProd, sigma, epsilon):
        """
        """
        self.particle1 = particle1
        self.particle2 = particle2
        self.chargeProd = chargeProd
        self.sigma = sigma
        self.epsilon = epsilon

        return
#=============================================================================================
# BondForce
#=============================================================================================

class BondForce(Force):


    def __init__(self, force=None):
        """
        Create a BondForce.

        """
        # Initialize defaults.
        self.bonds = list()

        # Populate data structures from swig object, if specified
        if force is not None:
            self._copyDataUsingInterface(self, force)

        return

    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        dest.__init__()
        for index in range(src.getNumBonds()):
            args = src.getBondParameters(index)
            dest.addBond(*args)

        return

    def asSwig(self):
        """
        Construct a Swig object.

        """
        force = openmm.BondForce()
        self._copyDataUsingInterface(force, self)
        return force

    def getNumBonds(self):
        """
        Get the number of bond stretch terms in the potential function

        """
        return len(self.bonds)

    @accepts_compatible_units(None, None, units.nanometers, units.kilojoules_per_mole / units.nanometers**2)
    def addBond(self, particle1, particle2, length, k):
        """
        Add a bond term to the force field.
        *
        @param particle1 the index of the first particle connected by the bond
        @param particle2 the index of the second particle connected by the bond
        @param length    the equilibrium length of the bond
        @param k         the harmonic force constant for the bond
        @return the index of the bond that was added

        """
        bond = BondForceBondInfo(particle1, particle2, length, k)
        self.bonds.append(bond)
        return

    def getBondParameters(self, index):
        """
        Get the force field parameters for a bond term.

        @param index     the index of the bond for which to get parameters
        @returns particle1 the index of the first particle connected by the bond
        @returns particle2 the index of the second particle connected by the bond
        @returns length    the equilibrium length of the bond
        @returns k         the force constant for the bond

        """
        # TODO: Return deep copy?
        bond = self.bonds[index]
        return (bond.particle1, bond.particle2, bond.length, bond.k)

    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole / units.nanometers**2)
    def setBondParameters(self, index, particle1, particle2, length, k):
        """
        Set the force field parameters for a bond term.

        @param index     the index of the bond for which to set parameters
        @param particle1 the index of the first particle connected by the bond
        @param particle2 the index of the second particle connected by the bond
        @param length    the equilibrium length of the bond
        @param k         the harmonic force constant for the bond
        """

        bond = BondForceBondInfo(particle1, particle2, length, k)
        self.bonds[index] = bond
        return

    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def nbonds(self):
        return len(self.bonds)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the System object.

        """

        r = ""

        # Show bonds.
        if (self.nbonds > 0):
            r += "Bonds:\n"
            r += "%8s %10s %10s %24s %24s\n" % ("index", "particle1", "particle2", "length", "k")
            for (index, bond) in enumerate(self.bonds):
                r += "%8d %10d %10d %24s %24s\n" % (index, bond.particle1, bond.particle2, str(bond.length), str(bond.k))
            r += "\n"
            r += "\n"

        return r

    def _appendForce(self, force, offset):
        """
        Append atoms defined in another force of the same type.

        @param force      the force to be appended
        @param offset     integral offset for atom number

        """
        # Check to make sure both forces are compatible.
        if type(self) != type(force):
            raise ValueError("other force object must be of identical Force subclass")

        # Combine systems.
        for bond in force.bonds:
            bond.particle1 += offset
            bond.particle2 += offset
            self.bonds.append(bond)

        return

class BondForceBondInfo(object):
    """
    Information about a bond.
    """
    @accepts_compatible_units(None, None, units.nanometers, units.kilojoules_per_mole / units.nanometers**2)
    def __init__(self, particle1, particle2, length, k):
        # Store data.
        self.particle1 = particle1
        self.particle2 = particle2
        self.length = length
        self.k = k
        return

#=============================================================================================
# HarmonicBondForce
#=============================================================================================

class HarmonicBondForce(Force):
    """
    This class implements an interaction between pairs of particles that varies harmonically with the distance
    between them.  To use it, create a HarmonicBondForce object then call addBond() once for each bond.  After
    a bond has been added, you can modify its force field parameters by calling setBondParameters().

    EXAMPLE

    Create and populate a HarmonicBondForce object.

    >>> bondforce = HarmonicBondForce()
    >>> bondforce.addBond(0, 1, 1.0*units.angstrom, 1.0*units.kilocalories_per_mole/units.angstrom**2)
    >>> bondforce.addBond(0, 2, 1.0*units.angstrom, 1.0*units.kilocalories_per_mole/units.angstrom**2)
    >>> bondforce.addBond(0, 3, 1.0*units.angstrom, 1.0*units.kilocalories_per_mole/units.angstrom**2)

    Create a Swig object.

    >>> bondforce_proxy = bondforce.asSwig()

    Create a deep copy.

    >>> bondforce_copy = copy.deepcopy(bondforce)

    Append a set of bonds.

    >>> bondforce_copy._appendForce(bondforce, 3)

    """

    def __init__(self, force=None):
        """
        Create a HarmonicBondForce.

        """
        # Initialize defaults.
        self.bonds = list()

        # Populate data structures from swig object, if specified
        if force is not None:
            self._copyDataUsingInterface(self, force)

        return

    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        dest.__init__()
        for index in range(src.getNumBonds()):
            args = src.getBondParameters(index)
            dest.addBond(*args)

        return

    def asSwig(self):
        """
        Construct a Swig object.

        """
        force = openmm.HarmonicBondForce()
        self._copyDataUsingInterface(force, self)
        return force
 
    def getNumBonds(self):
        """
        Get the number of harmonic bond stretch terms in the potential function

        """
        return len(self.bonds)

    @accepts_compatible_units(None, None, units.nanometers, units.kilojoules_per_mole / units.nanometers**2)
    def addBond(self, particle1, particle2, length, k):
        """
        Add a bond term to the force field.
        *
        @param particle1 the index of the first particle connected by the bond
        @param particle2 the index of the second particle connected by the bond
        @param length    the equilibrium length of the bond
        @param k         the harmonic force constant for the bond
        @return the index of the bond that was added

        """
        bond = HarmonicBondForceBondInfo(particle1, particle2, length, k)
        self.bonds.append(bond)
        return

    def getBondParameters(self, index):
        """
        Get the force field parameters for a bond term.

        @param index     the index of the bond for which to get parameters
        @returns particle1 the index of the first particle connected by the bond
        @returns particle2 the index of the second particle connected by the bond
        @returns length    the equilibrium length of the bond
        @returns k         the harmonic force constant for the bond

        """
        # TODO: Return deep copy?
        bond = self.bonds[index]
        return (bond.particle1, bond.particle2, bond.length, bond.k)

    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole / units.nanometers**2)
    def setBondParameters(self, index, particle1, particle2, length, k):
        """
        Set the force field parameters for a bond term.

        @param index     the index of the bond for which to set parameters
        @param particle1 the index of the first particle connected by the bond
        @param particle2 the index of the second particle connected by the bond
        @param length    the equilibrium length of the bond
        @param k         the harmonic force constant for the bond
        """

        bond = HarmonicBondForceBondInfo(particle1, particle2, length, k)
        self.bonds[index] = bond
        return

    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def nbonds(self):
        return len(self.bonds)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the System object.

        """

        r = ""

        # Show bonds.
        if (self.nbonds > 0):
            r += "Bonds:\n"
            r += "%8s %10s %10s %24s %24s\n" % ("index", "particle1", "particle2", "length", "k")
            for (index, bond) in enumerate(self.bonds):
                r += "%8d %10d %10d %24s %24s\n" % (index, bond.particle1, bond.particle2, str(bond.length), str(bond.k))
            r += "\n"
            r += "\n"

        return r

    def _appendForce(self, force, offset):
        """
        Append atoms defined in another force of the same type.

        @param force      the force to be appended
        @param offset     integral offset for atom number

        """
        # Check to make sure both forces are compatible.
        if type(self) != type(force):
            raise ValueError("other force object must be of identical Force subclass")

        # Combine systems.
        for bond in force.bonds:
            bond.particle1 += offset
            bond.particle2 += offset
            self.bonds.append(bond)

        return

class HarmonicBondForceBondInfo(object):
    """
    Information about a harmonic bond.
    """
    @accepts_compatible_units(None, None, units.nanometers, units.kilojoules_per_mole / units.nanometers**2)
    def __init__(self, particle1, particle2, length, k):
        # Store data.
        self.particle1 = particle1
        self.particle2 = particle2
        self.length = length
        self.k = k
        return

#=============================================================================================
# HarmonicAngleForce
#=============================================================================================

class HarmonicAngleForce(Force):
    """
    This class implements an interaction between groups of three particles that varies harmonically with the angle
    between them.  To use it, create a HarmonicAngleForce object then call addAngle() once for each angle.  After
    an angle has been added, you can modify its force field parameters by calling setAngleParameters().

    EXAMPLE

    Create and populate a HarmonicAngleForce object.

    >>> angle = 120.0 * units.degree
    >>> k = 0.01 * units.kilocalories_per_mole / units.degree**2
    >>> angleforce = HarmonicAngleForce()
    >>> angleforce.addAngle(0, 1, 2, angle, k)

    Create a Swig object.

    >>> angleforce_proxy = angleforce.asSwig()

    Create a deep copy.

    >>> angleforce_copy = copy.deepcopy(angleforce)

    Append a set of angles.

    >>> angleforce_copy._appendForce(angleforce, 3)

    """

    def __init__(self, force=None):
        """
        Create a HarmonicAngleForce.

        """
        # Initialize defaults.
        self.angles = list()

        # Populate data structures from swig object, if specified
        if force is not None:
            self._copyDataUsingInterface(self, force)

        return

    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        dest.__init__()
        for index in range(src.getNumAngles()):
            args = src.getAngleParameters(index)
            dest.addAngle(*args)

        return

    def asSwig(self):
        """
        Construct a Swig object.

        """
        force = openmm.HarmonicAngleForce()
        self._copyDataUsingInterface(force, self)
        return force

    def getNumAngles(self):
        """
        Get the number of harmonic angle stretch terms in the potential function

        """
        return len(self.angles)

    @accepts_compatible_units(None, None, None, units.radians, units.kilojoules_per_mole / units.radians**2)
    def addAngle(self, particle1, particle2, particle3, angle, k):
        """
        Add an angle term to the force field.

        @param particle1 the index of the first particle forming the angle
        @param particle2 the index of the second particle forming the angle
        @param particle3 the index of the third particle forming the angle
        @param angle     the equilibrium angle
        @param k         the harmonic force constant for the angle
        @return the index of the angle that was added

        """
        angle = HarmonicAngleForceAngleInfo(particle1, particle2, particle3, angle, k)
        self.angles.append(angle)
        return

    def getAngleParameters(self, index):
        """
        Get the force field parameters for an angle term.

        @param index     the index of the angle for which to get parameters
        @returns particle1 the index of the first particle forming the angle
        @returns particle2 the index of the second particle forming the angle
        @returns particle3 the index of the third particle forming the angle
        @returns angle     the equilibrium angle
        @returns k         the harmonic force constant for the angle

        """
        angle = self.angles[index]
        return (angle.particle1, angle.particle2, angle.particle3, angle.angle, angle.k)

    @accepts_compatible_units(None, None, None, None, units.radians, units.kilojoules_per_mole / units.radians**2)
    def setAngleParameters(self, index, particle1, particle2, particle3, angle, k):
        """
        Set the force field parameters for an angle term.

        @param index     the index of the angle for which to set parameters
        @param particle1 the index of the first particle forming the angle
        @param particle2 the index of the second particle forming the angle
        @param particle3 the index of the third particle forming the angle
        @param angle     the equilibrium angle
        @param k         the harmonic force constant for the angle
        """

        angle = HarmonicAngleForceAngleInfo(particle1, particle2, particle3, angle, k)
        self.angles[index] = angle
        return

    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def nangles(self):
        return len(self.angles)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the System object.

        """

        r = ""

        # Show angles.
        if (self.nangles > 0):
            r += "Angles:\n"
            r += "%8s %10s %10s %10s %24s %24s\n" % ("index", "particle1", "particle2", "particle3", "length", "k")
            for (index, angle) in enumerate(self.angles):
                r += "%8d %10d %10d %10d %24s %24s\n" % (index, angle.particle1, angle.particle2, angle.particle3, str(angle.angle), str(angle.k))
            r += "\n"
            r += "\n"

        return r

    def _appendForce(self, force, offset):
        """
        Append atoms defined in another force of the same type.

        @param force      the force to be appended
        @param offset     integral offset for atom number

        """
        # Check to make sure both forces are compatible.
        if type(self) != type(force):
            raise ValueError("other force object must be of identical Force subclass")

        # Combine systems.
        for angle in force.angles:
            angle.particle1 += offset
            angle.particle2 += offset
            angle.particle3 += offset
            self.angles.append(angle)

        return

class HarmonicAngleForceAngleInfo(object):
    """
    Information about a harmonic angle.
    """
    @accepts_compatible_units(None, None, None, units.radians, units.kilojoules_per_mole / units.radians**2)
    def __init__(self, particle1, particle2, particle3, angle, k):
        # Store data.
        self.particle1 = particle1
        self.particle2 = particle2
        self.particle3 = particle3
        self.angle = angle
        self.k = k
        return


#=============================================================================================
# LJ1Force
#=============================================================================================

class LJ1Force(Force):


    def __init__(self, force=None):
        """
        Create a LJ1Force.

        """
        # Initialize defaults.
        self.torsions = list()

        # Populate data structures from swig object, if specified
        if force is not None:
            self._copyDataUsingInterface(self, RBforce)

        return

    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        dest.__init__()
        for index in range(src.getNumTorsions()):
            args = src.getTorsionParameters(index)
            dest.addTorsion(*args)

        return

    def asSwig(self):
        """
        Construct a Swig object.
        """
        force = openmm.LJ1TorsionForce()
        self._copyDataUsingInterface(force, self)
        return force

    def getNumTorsions(self):
        """
        Get the number of LJ1 pair terms in the potential function

        """
        return len(self.torsions)

    @accepts_compatible_units(None, None, None, None)
    def addTorsion(self, particle1, particle2, V, W):
        """
        Add a LJ1 pair term to the force field.

        @param particle1    the index of the first particle forming the pair
        @param particle2    the index of the second particle forming the pair
        @param V
        @param W
        @return the index of the torsion that was added

        """
        torsion = LJ1ForceLJ1Info(particle1, particle2, V, W)
        self.torsions.append(torsion) 
        return

    def getTorsionParameters(self, index):
        """
        Get the force field parameters for a LJ1 pair term.

        @param index        the index of the torsion for which to get parameters
        @returns particle1    the index of the first particle forming the pair
        @returns particle2    the index of the second particle forming the pair
        @returns V
        @returns W

        """
        torsion = self.torsions[index]
        return (torsion.particle1, torsion.particle2, torsion.V, torsion.W)

    @accepts_compatible_units(None, None, None, None, None)
    def setTorsionParameters(self, index, particle1, particle2, V, W):
        """
        Set the force field parameters for a LJ1 pair term.

        @param index        the index of the torsion for which to set parameters
        @param particle1    the index of the first particle forming the pair
        @param particle2    the index of the second particle forming the pair 
        @param V
        @param W

        """
        torsion = LJ1ForceLJ1Info(particle1, particle2, V, W)
        self.torsions[index] = torsion
        return

    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def ntorsions(self):
        return len(self.torsions)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the System object.

        """

        r = ""

        # Show torsions.
        if (self.ntorsions > 0):
            r += "Torsions:\n"
            r += "%8s %10s %10s %24s %24s\n" % ("index", "particle1", "particle2", "V", "W")
            for (index, torsion) in enumerate(self.torsions):
                r += "%8d %10d %10d %24s %24s\n" % (index, torsion.particle1, torsion.particle2, str(torsion.V), str(torsion.W))
            r += "\n"
            r += "\n"

        return r

    def _appendForce(self, force, offset):
        """
        Append atoms defined in another force of the same type.

        @param force      the force to be appended
        @param offset     integral offset for atom number

        """
        # Check to make sure both forces are compatible.
        if type(self) != type(force):
            raise ValueError("other force object must be of identical Force subclass")

        # Combine systems.
        for torsion in force.torsions:
            torsion.particle1 += offset
            torsion.particle2 += offset
            torsion.particle3 += offset
            torsion.particle4 += offset
            self.torsions.append(torsion)

        return

class LJ1ForceLJ1Info(object):
    """
    Information about a LJ1 pair.

    """
    @accepts_compatible_units(None, None, None, None)
    def __init__(self, particle1, particle2, V, W):
        self.particle1 = particle1
        self.particle2 = particle2
        self.V = V
        self.W = W

        return


#=============================================================================================
# LJ2Force
#=============================================================================================

class LJ2Force(Force):


    def __init__(self, force=None):
        """
        Create a LJ2Force.

        """
        # Initialize defaults.
        self.torsions = list()

        # Populate data structures from swig object, if specified
        if force is not None:
            self._copyDataUsingInterface(self, RBforce)

        return

    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        dest.__init__()
        for index in range(src.getNumTorsions()):
            args = src.getTorsionParameters(index)
            dest.addTorsion(*args)

        return

    def asSwig(self):
        """
        Construct a Swig object.

        """
        force = openmm.LJ2TorsionForce()
        self._copyDataUsingInterface(force, self)
        return force

    def getNumTorsions(self):
        """
        Get the number of LJ2 pair terms in the potential function

        """
        return len(self.torsions)

    @accepts_compatible_units(None, None, None, units.elementary_charge, units.elementary_charge,  None, None)
    def addTorsion(self, particle1, particle2, fudgeQQ, qi, qj, V, W):
        """
        Add a LJ2 pair term to the force field.

        @param particle1    the index of the first particle forming the pair
        @param particle2    the index of the second particle forming the pair
        @param fudgeQQ
        @param qi
        @param qj
        @param V
        @param W
        @return the index of the torsion that was added

        """
        torsion = LJ2ForceLJ2Info(particle1, particle2, fudgeQQ, qi, qj, V, W)
        self.torsions.append(torsion) 
        return

    def getTorsionParameters(self, index):
        """
        Get the force field parameters for a LJ2 pair term.

        @param index        the index of the torsion for which to get parameters
        @returns particle1    the index of the first particle forming the pair
        @returns particle2    the index of the second particle forming the pair
        @returns fudgeQQ
        @returns qi
        @returns qj
        @returns V
        @returns W


        """
        torsion = self.torsions[index]
        return (torsion.particle1, torsion.particle2, torsion.fudgeQQ, torsion.qi, torsion.qj, torsion.V, torsion.W)

    @accepts_compatible_units(None, None, None, None, units.elementary_charge, units.elementary_charge, None, None)
    def setTorsionParameters(self, index, particle1, particle2, fudgeQQ, qi, qj, V, W):
        """
        Set the force field parameters for a LJ1 pair term.

        @param index        the index of the torsion for which to set parameters
        @param particle1    the index of the first particle forming the pair
        @param particle2    the index of the second particle forming the pair 
        @param fudgeQQ
        @param qi
        @param qj
        @param V
        @param W

        """
        torsion = LJ2ForceLJ2Info(particle1, particle2, fudgeQQ, qi, qj, V, W)
        self.torsions[index] = torsion
        return

    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def ntorsions(self):
        return len(self.torsions)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the System object.

        """

        r = ""

        # Show torsions.
        if (self.ntorsions > 0):
            r += "Torsions:\n"
            r += "%8s %10s %10s %24s %24s %24s %24s %24s\n" % ("index", "particle1", "particle2", "fudgeQQ", "qi", "qj", "V", "W")
            for (index, torsion) in enumerate(self.torsions):
                r += "%8d %10d %10d %24s %24s %24s %24s %24s\n" % (index, torsion.particle1, torsion.particle2, str(torsion.fudgeQQ), str(torsion.qi), str(torsion.qj), str(torsion.V), str(torsion.W))
            r += "\n"
            r += "\n"

        return r

    def _appendForce(self, force, offset):
        """
        Append atoms defined in another force of the same type.

        @param force      the force to be appended
        @param offset     integral offset for atom number

        """
        # Check to make sure both forces are compatible.
        if type(self) != type(force):
            raise ValueError("other force object must be of identical Force subclass")

        # Combine systems.
        for torsion in force.torsions:
            torsion.particle1 += offset
            torsion.particle2 += offset
            torsion.particle3 += offset
            torsion.particle4 += offset
            self.torsions.append(torsion)

        return

class LJ2ForceLJ1Info(object):
    """
    Information about a LJ2 pair.
    """
    @accepts_compatible_units(None, None, None, units.elementary_charge, units.elementary_charge, None, None)
    def __init__(self, particle1, particle2, fudgeQQ, qi, qj, V, W):
        self.particle1 = particle1
        self.particle2 = particle2
        self.fudgeQQ = fudgeQQ
        self.qi = qi
        self.qj = qj
        self.V = V
        self.W = W

        return


#=============================================================================================
# LJNBForce
#=============================================================================================

class LJNBForce(Force):


    def __init__(self, force=None):
        """
        Create a LJNBForce.

        """
        # Initialize defaults.
        self.torsions = list()

        # Populate data structures from swig object, if specified
        if force is not None:
            self._copyDataUsingInterface(self, RBforce)

        return

    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        dest.__init__()
        for index in range(src.getNumTorsions()):
            args = src.getTorsionParameters(index)
            dest.addTorsion(*args)

        return

    def asSwig(self):
        """
        Construct a Swig object.

        """
        force = openmm.LJNBTorsionForce()
        self._copyDataUsingInterface(force, self)
        return force

    def getNumTorsions(self):
        """
        Get the number of LJNB pair terms in the potential function

        """
        return len(self.torsions)

    @accepts_compatible_units(None, None, units.elementary_charge, units.elementary_charge, None,  None)
    def addTorsion(self, particle1, particle2, qi, qj, V, W):
        """
        Add a LJ2 pair term to the force field.

        @param particle1    the index of the first particle forming the pair
        @param particle2    the index of the second particle forming the pair
        @param qi
        @param qj
        @param V
        @param W
        @return the index of the torsion that was added

        """
        torsion = LJNBForceLJNBInfo(particle1, particle2, qi, qj, V, W)
        self.torsions.append(torsion) 
        return

    def getTorsionParameters(self, index):
        """
        Get the force field parameters for a LJNB pair term.

        @param index        the index of the torsion for which to get parameters
        @returns particle1    the index of the first particle forming the pair
        @returns particle2    the index of the second particle forming the pair
        @returns qi
        @returns qj
        @returns V
        @returns W


        """
        torsion = self.torsions[index]
        return (torsion.particle1, torsion.particle2, torsion.qi, torsion.qj, torsion.V, torsion.W)

    @accepts_compatible_units(None, None, None, units.elementary_charge, units.elementary_charge, None, None)
    def setTorsionParameters(self, index, particle1, particle2, qi, qj, V, W):
        """
        Set the force field parameters for a LJ1 pair term.

        @param index        the index of the torsion for which to set parameters
        @param particle1    the index of the first particle forming the pair
        @param particle2    the index of the second particle forming the pair 
        @param qi
        @param qj
        @param V
        @param W

        """
        torsion = LJNBForceLJNBInfo(particle1, particle2, qi, qj, V, W)
        self.torsions[index] = torsion
        return

    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def ntorsions(self):
        return len(self.torsions)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the System object.

        """

        r = ""

        # Show torsions.
        if (self.ntorsions > 0):
            r += "Torsions:\n"
            r += "%8s %10s %10s %24s %24s %24s %24s\n" % ("index", "particle1", "particle2", "qi", "qj", "V", "W")
            for (index, torsion) in enumerate(self.torsions):
                r += "%8d %10d %10d %24s %24s %24s %24s\n" % (index, torsion.particle1, torsion.particle2, str(torsion.qi), str(torsion.qj), str(torsion.V), str(torsion.W))
            r += "\n"
            r += "\n"

        return r

    def _appendForce(self, force, offset):
        """
        Append atoms defined in another force of the same type.

        @param force      the force to be appended
        @param offset     integral offset for atom number

        """
        # Check to make sure both forces are compatible.
        if type(self) != type(force):
            raise ValueError("other force object must be of identical Force subclass")

        # Combine systems.
        for torsion in force.torsions:
            torsion.particle1 += offset
            torsion.particle2 += offset
            torsion.particle3 += offset
            torsion.particle4 += offset
            self.torsions.append(torsion)

        return

class LJNBForceLJNBInfo(object):
    """
    Information about a LJNB pair.

    """
    @accepts_compatible_units(None, None, units.elementary_charge, units.elementary_charge, None,  None)
    def __init__(self, particle1, particle2, qi, qj, V, W):
        self.particle1 = particle1
        self.particle2 = particle2
        self.qi = qi
        self.qj = qj
        self.V = V
        self.W = W

        return

#=============================================================================================
# PeriodicTorsionForce
#=============================================================================================

class PeriodicTorsionForce(Force):
    """
    This class implements an interaction between groups of four particles that varies periodically with the torsion angle
    between them.  To use it, create a PeriodicTorsionForce object then call addTorsion() once for each torsion.  After
    a torsion has been added, you can modify its force field parameters by calling setTorsionParameters().

    EXAMPLE

    Create and populate a PeriodicTorsionForce object.

    >>> periodicity = 3
    >>> phase = 30 * units.degrees
    >>> k = 1.0 * units.kilocalories_per_mole
    >>> torsionforce = PeriodicTorsionForce()
    >>> torsionforce.addTorsion(0, 1, 2, 3, periodicity, phase, k)

    Create a Swig object.

    >>> torsionforce_proxy = torsionforce.asSwig()

    Create a deep copy.

    >>> torsionforce_copy = copy.deepcopy(torsionforce)

    Append a set of angles.

    >>> torsionforce_copy._appendForce(torsionforce, 4)

    """

    def __init__(self, force=None):
        """
        Create a PeriodicTorsionForce.

        """
        # Initialize defaults.
        self.torsions = list()

        # Populate data structures from swig object, if specified
        if force is not None:
            self._copyDataUsingInterface(self, force)

        return

    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        dest.__init__()
        for index in range(src.getNumTorsions()):
            args = src.getTorsionParameters(index)
            dest.addTorsion(*args)

        return

    def asSwig(self):
        """
        Construct a Swig object.

        """
        force = openmm.PeriodicTorsionForce()
        self._copyDataUsingInterface(force, self)
        return force

    def getNumTorsions(self):
        """
        Get the number of periodic torsion terms in the potential function

        """
        return len(self.torsions)

    @accepts_compatible_units(None, None, None, None, None, units.radians, units.kilojoules_per_mole)
    def addTorsion(self, particle1, particle2, particle3, particle4, periodicity, phase, k):
        """
        Add a periodic torsion term to the force field.

        @param particle1    the index of the first particle forming the torsion
        @param particle2    the index of the second particle forming the torsion
        @param particle3    the index of the third particle forming the torsion
        @param particle3    the index of the fourth particle forming the torsion
        @param periodicity  the periodicity of the torsion
        @param phase        the phase offset of the torsion
        @param k            the force constant for the torsion
        @return the index of the torsion that was added

        """
        torsion = PeriodicTorsionForcePeriodicTorsionInfo(particle1, particle2, particle3, particle4, periodicity, phase, k)
        self.torsions.append(torsion) 
        return

    def getTorsionParameters(self, index):
        """
        Get the force field parameters for a periodic torsion term.

        @param index        the index of the torsion for which to get parameters
        @returns particle1    the index of the first particle forming the torsion
        @returns particle2    the index of the second particle forming the torsion
        @returns particle3    the index of the third particle forming the torsion
        @returns particle3    the index of the fourth particle forming the torsion
        @returns periodicity  the periodicity of the torsion
        @returns phase        the phase offset of the torsion
        @returns k            the force constant for the torsion

        """
        torsion = self.torsions[index]
        return (torsion.particle1, torsion.particle2, torsion.particle3, torsion.particle4, torsion.periodicity, torsion.phase, torsion.k)

    @accepts_compatible_units(None, None, None, None, None, None, units.radians, units.kilojoules_per_mole)
    def setTorsionParameters(self, index, particle1, particle2, particle3, particle4, periodicity, phase, k):
        """
        Set the force field parameters for a periodic torsion term.

        @param index        the index of the torsion for which to set parameters
        @param particle1    the index of the first particle forming the torsion
        @param particle2    the index of the second particle forming the torsion
        @param particle3    the index of the third particle forming the torsion
        @param particle3    the index of the fourth particle forming the torsion
        @param periodicity  the periodicity of the torsion
        @param phase        the phase offset of the torsion
        @param k            the force constant for the torsion

        """
        torsion = PeriodicTorsionForcePeriodicTorsionInfo(particle1, particle2, particle3, particle4, periodicity, phase, k)
        self.torsions[index] = torsion
        return

    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def ntorsions(self):
        return len(self.torsions)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the System object.

        """

        r = ""

        # Show torsions.
        if (self.ntorsions > 0):
            r += "Torsions:\n"
            r += "%8s %10s %10s %10s %10s %10s %24s %24s\n" % ("index", "particle1", "particle2", "particle3", "particle4", "periodicity", "phase", "k")
            for (index, torsion) in enumerate(self.torsions):
                r += "%8d %10d %10d %10d %10d %10d %24s %24s\n" % (index, torsion.particle1, torsion.particle2, torsion.particle3, torsion.particle4,  torsion.periodicity, str(torsion.phase), str(torsion.k))
            r += "\n"
            r += "\n"

        return r

    def _appendForce(self, force, offset):
        """
        Append atoms defined in another force of the same type.

        @param force      the force to be appended
        @param offset     integral offset for atom number

        """
        # Check to make sure both forces are compatible.
        if type(self) != type(force):
            raise ValueError("other force object must be of identical Force subclass")

        # Combine systems.
        for torsion in force.torsions:
            torsion.particle1 += offset
            torsion.particle2 += offset
            torsion.particle3 += offset
            torsion.particle4 += offset
            self.torsions.append(torsion)

        return

class PeriodicTorsionForcePeriodicTorsionInfo(object):
    """
    Information about a periodic torsion.

    """
    @accepts_compatible_units(None, None, None, None, None, units.radians, units.kilojoules_per_mole)
    def __init__(self, particle1, particle2, particle3, particle4, periodicity, phase, k):
        self.particle1 = particle1
        self.particle2 = particle2
        self.particle3 = particle3
        self.particle4 = particle4
        self.periodicity = periodicity
        self.phase = phase
        self.k = k
        return

#=============================================================================================
# NonPeriodicTorsionForce
#=============================================================================================

class NonPeriodicTorsionForce(Force):


    def __init__(self, force=None):
        """
        Create a PeriodicTorsionForce.

        """
        # Initialize defaults.
        self.torsions = list()

        # Populate data structures from swig object, if specified
        if force is not None:
            self._copyDataUsingInterface(self, force)

        return

    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        dest.__init__()
        for index in range(src.getNumTorsions()):
            args = src.getTorsionParameters(index)
            dest.addTorsion(*args)

        return

    def asSwig(self):
        """
        Construct a Swig object.

        """
        force = openmm.PeriodicTorsionForce()
        self._copyDataUsingInterface(force, self)
        return force

    def getNumTorsions(self):
        """
        Get the number of periodic torsion terms in the potential function

        """
        return len(self.torsions)

    @accepts_compatible_units(None, None, None, None, units.radians, units.kilojoules_per_mole)
    def addTorsion(self, particle1, particle2, particle3, particle4, phase, k):
        """
        Add a periodic torsion term to the force field.

        @param particle1    the index of the first particle forming the torsion
        @param particle2    the index of the second particle forming the torsion
        @param particle3    the index of the third particle forming the torsion
        @param particle3    the index of the fourth particle forming the torsion
        @param phase        the phase offset of the torsion
        @param k            the force constant for the torsion
        @return the index of the torsion that was added

        """
        torsion = NonPeriodicTorsionForceNonPeriodicTorsionInfo(particle1, particle2, particle3, particle4, phase, k)
        self.torsions.append(torsion) 
        return

    def getTorsionParameters(self, index):
        """
        Get the force field parameters for a periodic torsion term.

        @param index        the index of the torsion for which to get parameters
        @returns particle1    the index of the first particle forming the torsion
        @returns particle2    the index of the second particle forming the torsion
        @returns particle3    the index of the third particle forming the torsion
        @returns particle3    the index of the fourth particle forming the torsion
        @returns phase        the phase offset of the torsion
        @returns k            the force constant for the torsion

        """
        torsion = self.torsions[index]
        return (torsion.particle1, torsion.particle2, torsion.particle3, torsion.particle4, torsion.phase, torsion.k)

    @accepts_compatible_units(None, None, None, None, None, units.radians, units.kilojoules_per_mole)
    def setTorsionParameters(self, index, particle1, particle2, particle3, particle4, phase, k):
        """
        Set the force field parameters for a periodic torsion term.

        @param index        the index of the torsion for which to set parameters
        @param particle1    the index of the first particle forming the torsion
        @param particle2    the index of the second particle forming the torsion
        @param particle3    the index of the third particle forming the torsion
        @param particle3    the index of the fourth particle forming the torsion
        @param phase        the phase offset of the torsion
        @param k            the force constant for the torsion

        """
        torsion = NonPeriodicTorsionForceNonPeriodicTorsionInfo(particle1, particle2, particle3, particle4, phase, k)
        self.torsions[index] = torsion
        return

    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def ntorsions(self):
        return len(self.torsions)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the System object.

        """

        r = ""

        # Show torsions.
        if (self.ntorsions > 0):
            r += "Torsions:\n"
            r += "%8s %10s %10s %10s %10s %24s %24s\n" % ("index", "particle1", "particle2", "particle3", "particle4", "phase", "k")
            for (index, torsion) in enumerate(self.torsions):
                r += "%8d %10d %10d %10d %10d %24s %24s\n" % (index, torsion.particle1, torsion.particle2, torsion.particle3, torsion.particle4, str(torsion.phase), str(torsion.k))
            r += "\n"
            r += "\n"

        return r

    def _appendForce(self, force, offset):
        """
        Append atoms defined in another force of the same type.

        @param force      the force to be appended
        @param offset     integral offset for atom number

        """
        # Check to make sure both forces are compatible.
        if type(self) != type(force):
            raise ValueError("other force object must be of identical Force subclass")

        # Combine systems.
        for torsion in force.torsions:
            torsion.particle1 += offset
            torsion.particle2 += offset
            torsion.particle3 += offset
            torsion.particle4 += offset
            self.torsions.append(torsion)

        return

class NonPeriodicTorsionForceNonPeriodicTorsionInfo(object):
    """
    Information about a periodic torsion.

    """
    @accepts_compatible_units(None, None, None, None, units.radians, units.kilojoules_per_mole)
    def __init__(self, particle1, particle2, particle3, particle4, phase, k):
        self.particle1 = particle1
        self.particle2 = particle2
        self.particle3 = particle3
        self.particle4 = particle4
        self.phase = phase
        self.k = k
        return

#=============================================================================================
# RBTorsionForce
#=============================================================================================

class RBTorsionForce(Force):


    def __init__(self, force=None):
        """
        Create a RBTorsionForce.

        """
        # Initialize defaults.
        self.torsions = list()

        # Populate data structures from swig object, if specified
        if force is not None:
            self._copyDataUsingInterface(self, force)

        return

    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        dest.__init__()
        for index in range(src.getNumTorsions()):
            args = src.getTorsionParameters(index)
            dest.addTorsion(*args)

        return

    def asSwig(self):
        """
        Construct a Swig object.

        """
        force = openmm.RBTorsionForce()
        self._copyDataUsingInterface(force, self)
        return force

    def getNumTorsions(self):
        """
        Get the number of RB torsion terms in the potential function

        """
        return len(self.torsions)

    @accepts_compatible_units(None, None, None, None, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole)
    def addTorsion(self, particle1, particle2, particle3, particle4, C0, C1, C2, C3, C4, C5):
        """
        Add a RB torsion term to the force field.

        @param particle1    the index of the first particle forming the torsion
        @param particle2    the index of the second particle forming the torsion
        @param particle3    the index of the third particle forming the torsion
        @param particle4    the index of the fourth particle forming the torsion
        @param C0
        @param C1
        @param C2
        @param C3
        @param C4
        @param C5
        @return the index of the torsion that was added

        """
        torsion = RBTorsionForceRBTorsionInfo(particle1, particle2, particle3, particle4, C0, C1, C2, C3, C4, C5)
        self.torsions.append(torsion)
        return

    def getTorsionParameters(self, index):
        """
        Get the force field parameters for a RB torsion term.

        @param index        the index of the torsion for which to get parameters
        @returns particle1    the index of the first particle forming the torsion
        @returns particle2    the index of the second particle forming the torsion
        @returns particle3    the index of the third particle forming the torsion
        @returns particle3    the index of the fourth particle forming the torsion
        @returns C0
        @returns C1
        @returns C2
        @returns C3
        @returns C4
        @returns C5

        """
        torsion = self.torsions[index]
        return (torsion.particle1, torsion.particle2, torsion.particle3, torsion.particle4, torsion.C0, torsion.C1, torsion.C2, torsion.C3, torsion.C4, torsion.C5)

    @accepts_compatible_units(None, None, None, None, None, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole)
    def setTorsionParameters(self, index, particle1, particle2, particle3, particle4, C0, C1, C2, C3, C4, C5):
        """
        Set the force field parameters for a periodic torsion term.

        @param index        the index of the torsion for which to set parameters
        @param particle1    the index of the first particle forming the torsion
        @param particle2    the index of the second particle forming the torsion
        @param particle3    the index of the third particle forming the torsion
        @param particle3    the index of the fourth particle forming the torsion
        @param C0
        @param C1
        @param C2
        @param C3
        @param C4
        @param C5

        """
        torsion =  RBForceRBTorsionInfo(particle1, particle2, particle3, particle4, C0, C1, C2, C3, C4, C5)
        self.torsions[index] = torsion
        return

    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def ntorsions(self):
        return len(self.torsions)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the System object.

        """

        r = ""

        # Show torsions.
        if (self.ntorsions > 0):
            r += "Torsions:\n"
            r += "%8s %10s %10s %10s %10s %24s %24s %24s %24s %24s %24s\n" % ("index", "particle1", "particle2", "particle3", "particle4", "C0", "C1", "C2", "C3", "C4", "C5")
            for (index, torsion) in enumerate(self.torsions):
                r += "%8d %10d %10d %10d %10d %24s %24s %24s %24s %24s\n" % (index, torsion.particle1, torsion.particle2, torsion.particle3, torsion.particle4, str(torsion.C0), str(torsion.C1), str(torsion.C2), str(torsion.C3), str(torsion.C4), str(torsion.C5))
            r += "\n"
            r += "\n"

        return r

    def _appendForce(self, force, offset):
        """
        Append atoms defined in another force of the same type.

        @param force      the force to be appended
        @param offset     integral offset for atom number

        """
        # Check to make sure both forces are compatible.
        if type(self) != type(force):
            raise ValueError("other force object must be of identical Force subclass")

        # Combine systems.
        for torsion in force.torsions:
            torsion.particle1 += offset
            torsion.particle2 += offset
            torsion.particle3 += offset
            torsion.particle4 += offset
            self.torsions.append(torsion)

        return

class RBTorsionForceRBTorsionInfo(object):
    """
    Information about a  RB torsion.

    """
    @accepts_compatible_units(None, None, None, None, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole)
    def __init__(self, particle1, particle2, particle3, particle4, C0, C1, C2, C3, C4, C5):
        self.particle1 = particle1
        self.particle2 = particle2
        self.particle3 = particle3
        self.particle4 = particle4
        self.C0 = C0
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4
        self.C5 = C5
        return

#=============================================================================================
# FourierTorsionForce
#=============================================================================================

class FourierTorsionForce(Force):


    def __init__(self, force=None):
        """
        Create a FourierTorsionForce.

        """
        # Initialize defaults.
        self.torsions = list()

        # Populate data structures from swig object, if specified
        if force is not None:
            self._copyDataUsingInterface(self, RBforce)

        return

    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        dest.__init__()
        for index in range(src.getNumTorsions()):
            args = src.getTorsionParameters(index)
            dest.addTorsion(*args)

        return

    def asSwig(self):
        """
        Construct a Swig object.

        """
        force = openmm.FourierTorsionForce()
        self._copyDataUsingInterface(force, self)
        return force

    def getNumTorsions(self):
        """
        Get the number of Fourier torsion terms in the potential function

        """
        return len(self.torsions)

    @accepts_compatible_units(None, None, None, None, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole)
    def addTorsion(self, particle1, particle2, particle3, particle4, C1, C2, C3, C4):
        """
        Add a RB torsion term to the force field.

        @param particle1    the index of the first particle forming the torsion
        @param particle2    the index of the second particle forming the torsion
        @param particle3    the index of the third particle forming the torsion
        @param particle4    the index of the fourth particle forming the torsion
        @param C1
        @param C2
        @param C3
        @param C4
        @return the index of the torsion that was added

        """
        torsion = FourierTorsionForceFourierTorsionInfo(particle1, particle2, particle3, particle4, C1, C2, C3, C4)
        self.torsions.append(torsion) 
        return

    def getTorsionParameters(self, index):
        """
        Get the force field parameters for a fourier torsion term.

        @param index        the index of the torsion for which to get parameters
        @returns particle1    the index of the first particle forming the torsion
        @returns particle2    the index of the second particle forming the torsion
        @returns particle3    the index of the third particle forming the torsion
        @returns particle3    the index of the fourth particle forming the torsion
        @returns C1
        @returns C2
        @returns C3RB
        @returns C4

        """
        torsion = self.torsions[index]
        return (torsion.particle1, torsion.particle2, torsion.particle3, torsion.particle4, torsion.C1, torsion.C2, torsion.C3, torsion.C4)

    @accepts_compatible_units(None, None, None, None, None, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole)
    def setTorsionParameters(self, index, particle1, particle2, particle3, particle4, C1, C2, C3, C4):
        """
        Set the force field parameters for a fourier torsion term.

        @param index        the index of the torsion for which to set parameters
        @param particle1    the index of the first particle forming the torsion
        @param particle2    the index of the second particle forming the torsion
        @param particle3    the index of the third particle forming the torsion
        @param particle3    the index of the fourth particle forming the torsion
        @param C1
        @param C2
        @param C3
        @param C4

        """
        torsion = FourierTorsionForceFourierTorsionInfo(particle1, particle2, particle3, particle4, C1, C2, C3, C4)
        self.torsions[index] = torsion
        return

    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def ntorsions(self):
        return len(self.torsions)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the System object.

        """

        r = ""

        # Show torsions.
        if (self.ntorsions > 0):
            r += "Torsions:\n"
            r += "%8s %10s %10s %10s %10s %24s %24s %24s %24s\n" % ("index", "particle1", "particle2", "particle3", "particle4", "C1", "C2", "C3", "C4")
            for (index, torsion) in enumerate(self.torsions):
                r += "%8d %10d %10d %10d %10d %24s %24s %24s %24s\n" % (index, torsion.particle1, torsion.particle2, torsion.particle3, str(torsion.C1), str(torsion.C2), str(torsion.C3), str(torsion.C4))
            r += "\n"
            r += "\n"

        return r

    def _appendForce(self, force, offset):
        """
        Append atoms defined in another force of the same type.

        @param force      the force to be appended
        @param offset     integral offset for atom number

        """
        # Check to make sure both forces are compatible.
        if type(self) != type(force):
            raise ValueError("other force object must be of identical Force subclass")

        # Combine systems.
        for torsion in force.torsions:
            torsion.particle1 += offset
            torsion.particle2 += offset
            torsion.particle3 += offset
            torsion.particle4 += offset
            self.torsions.append(torsion)

        return

class FourierTorsionForceFourierTorsionInfo(object):
    """
    Information about a fourier torsion.

    """
    @accepts_compatible_units(None, None, None, None, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole)
    def __init__(self, particle1, particle2, particle3, particle4, C1, C2, C3, C4):
        self.particle1 = particle1
        self.particle2 = particle2
        self.particle3 = particle3
        self.particle4 = particle4
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4
        return



#=============================================================================================
# GBSAOBCForce
#=============================================================================================

class GBSAOBCForce(Force):
    """
    This class implements an implicit solvation force using the GBSA-OBC model.

    To use this class, create a GBSAOBCForce object, then call addParticle() once for each particle in the
    System to define its parameters.  The number of particles for which you define GBSA parameters must
    be exactly equal to the number of particles in the System, or else an exception will be thrown when you
    try to create a Context.  After a particle has been added, you can modify its force field parameters
    by calling setParticleParameters().

    EXAMPLE

    Create and populate a GBSAOBCForce object.

    >>> charge = 1.0 * units.elementary_charge
    >>> radius = 1.5 * units.angstrom
    >>> scalingFactor = 1.0
    >>> gbsaforce = GBSAOBCForce()
    >>> gbsaforce.addParticle(charge, radius, scalingFactor)
    0
    >>> gbsaforce.addParticle(charge, radius, scalingFactor)
    1

    Create a Swig object.

    >>> gbsaforce_proxy = gbsaforce.asSwig()

    Create a deep copy.

    >>> gbsaforce_copy = copy.deepcopy(gbsaforce)

    Append a set of particles.

    >>> gbsaforce_copy._appendForce(gbsaforce, 2)

    """

    NoCutoff = 0 #: No cutoff is applied to nonbonded interactions.  The full set of N^2 interactions is computed exactly. This necessarily means that periodic boundary conditions cannot be used.  This is the default.
    CutoffNonPeriodic = 1 #: Interactions beyond the cutoff distance are ignored.
    CutoffPeriodic = 2 #: Periodic boundary conditions are used, so that each particle interacts only with the nearest periodic copy of each other particle.  Interactions beyond the cutoff distance are ignored.

    def __init__(self, force=None):
        """
        Create a GBSAOBCForce.

        """

        # Initialize with defaults.
        self.particles = list()
        self.nonbondedMethod = GBSAOBCForce.NoCutoff
        self.cutoffDistance = units.Quantity(1.0, units.nanometer)
        self.solventDielectric = units.Quantity(78.3, units.dimensionless)
        self.soluteDielectric = units.Quantity(78.3, units.dimensionless)

        # Populate data structures from swig object, if specified
        if force is not None:
            self._copyDataUsingInterface(self, force)

        return

    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        dest.__init__()
        dest.setNonbondedMethod( src.getNonbondedMethod() )
        dest.setCutoffDistance( src.getCutoffDistance() )
        dest.setSoluteDielectric( src.getSoluteDielectric() )
        dest.setSolventDielectric( src.getSolventDielectric() )
        for index in range(src.getNumParticles()):
            args = src.getParticleParameters(index)
            dest.addParticle(*args)

        return

    def asSwig(self):
        """
        Construct a Swig object.

        """
        force = openmm.GBSAOBCForce()
        self._copyDataUsingInterface(force, self)
        return force

    def getNumParticles(self):
        """
        Get the number of particles in the system.

        """
        return len(self.particles)

    def addParticle(self, charge, radius, scalingFactor, nonPolarScalingFactor = 1.0):
        """
        Add the GBSA parameters for a particle.  This should be called once for each particle
        in the System.  When it is called for the i'th time, it specifies the parameters for the i'th particle.

        @param charge         the charge of the particle
        @param radius         the GBSA radius of the particle
        @param scalingFactor  the OBC scaling factor for the particle
        @return the index of the particle that was added

        """
        particle = GBSAOBCForceParticleInfo(charge, radius, scalingFactor)
        self.particles.append(particle)
        return (len(self.particles) - 1)

    def getParticleParameters(self, index):
        """
        Get the force field parameters for a particle.

        @param index          the index of the particle for which to get parameters
        @returns charge         the charge of the particle
        @returns radius         the GBSA radius of the particle
        @returns scalingFactor  the OBC scaling factor for the particle

        """
        particle = self.particles[index]
        return (particle.charge, particle.radius, particle.scalingFactor)

    def setParticleParameters(self, index, charge, radius, scalingFactor):
        """
        Set the force field parameters for a particle.

        @param index          the index of the particle for which to set parameters
        @param charge         the charge of the particle
        @param radius         the GBSA radius of the particle
        @param scalingFactor  the OBC scaling factor for the particle

        """
        particle = GBSAOBCForceParticleInfo(charge, radius, scalingFactor, nonPolarScalingFactor)
        self.particles[index] = particle
        return

    def getSolventDielectric(self):
        """
        Get the dielectric constant for the solvent.

        """
        return self.solventDielectric

    @accepts_compatible_units(units.dimensionless)
    def setSolventDielectric(self, solventDielectric):
        """
        Set the dielectric constant for the solvent.

        """
        self.solventDielectric = solventDielectric
        return

    def getSoluteDielectric(self):
        """
        Get the dielectric constant for the solute.

        """
        return self.soluteDielectric

    @accepts(float)
    def setSoluteDielectric(self, soluteDielectric):
        """
        Set the dielectric constant for the solute.

        """
        self.soluteDielectric = soluteDielectric
        return

    def getNonbondedMethod(self):
        """
        Get the method used for handling long range nonbonded interactions.

        """
        return self.nonbondedMethod

    def setNonbondedMethod(self, nonbondedMethod):
        """
        Set the method used for handling long range nonbonded interactions.

        """
        # TODO: Argument checking.
        self.nonbondedMethod = nonbondedMethod
        return

    def getCutoffDistance(self):
        """
        Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
        is NoCutoff, this value will have no effect.

        """
        return self.cutoffDistance

    @accepts_compatible_units(units.nanometer)
    def setCutoffDistance(self, cutoffDistance):
        """
        Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
        is NoCutoff, this value will have no effect.

        """
        self.cutoffDistance = cutoffDistance
        return

    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def nparticles(self):
        return len(self.particles)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the System object.

        """

        r = ""

        # Show settings.
        r += "nonbondedMethod: %s" % str(self.nonbondedMethod)
        r += "cutoffDistance: %s" % str(self.cutoffDistance)
        r += "solventDielectric: %s" % str(self.solventDielectric)
        r += "soluteDielectric: %s" % str(self.soluteDielectric)

        # Show particles.
        if (self.nparticles > 0):
            r += "Particles:\n"
            r += "%8s %24s %24s %24s %24s\n" % ("particle", "charge", "radius", "scalingFactor", "lambda")
            for (index, particle) in enumerate(self.particles):
                r += "%8d %24s %24s %24s %24f\n" % (index, str(particle.charge), str(particle.radius), str(particle.scalingFactor), 1.0)
            r += "\n"

        return r

    def _appendForce(self, force, offset):
        """
        Append atoms defined in another force of the same type.

        @param force      the force to be appended
        @param offset     integral offset for atom number

        EXAMPLES
        Create two force objects and append the second to the first.

        >>> force1 = GBSAOBCForce()
        >>> force2 = GBSAOBCForce()
        >>> charge = 1.0 * units.elementary_charge 
        >>> radius = 1.0 * units.angstrom
        >>> scalingFactor = 1.0
        >>> force1.addParticle(charge, radius, scalingFactor)
        0
        >>> force2.addParticle(charge, radius, scalingFactor)
        0
        >>> offset = 1
        >>> force1._appendForce(force2, offset)

        """
        # Check to make sure both forces are compatible.
        if type(self) != type(force):
            raise ValueError("other force object must be of identical Force subclass")

        # Check to make sure both forces have same settings.
        # TODO: This checking is probably more paranoid than we need.  Can this be relaxed?
        if (self.nonbondedMethod != force.nonbondedMethod):
            raise ValueError("other force has incompatible nonbondedMethod")
        if (self.cutoffDistance != force.cutoffDistance):
            raise ValueError("other force has incompatible cutoffDistance")
        if (self.solventDielectric != force.solventDielectric):
            raise ValueError("other force has incompatible solventDielectric")
        if (self.soluteDielectric != force.soluteDielectric):
            raise ValueError("other force has incompatible soluteDielectric")

        # Combine systems.
        for particle in force.particles:
            self.particles.append(particle)

        return

class GBSAOBCForceParticleInfo(object):
    @accepts_compatible_units(units.elementary_charge, units.nanometers, None)
    def __init__(self, charge, radius, scalingFactor):
        """
        """
        self.charge = charge
        self.radius = radius
        self.scalingFactor = scalingFactor
        return

#=============================================================================================
# CustomNonbondedForce
#=============================================================================================

class CustomNonbondedForce(Force):
    """
    This class implements nonbonded interactions between particles.  Unlike NonbondedForce, the functional form
    of the interaction is completely customizable, and may involve arbitrary algebraic expressions and tabulated
    functions.  It may depend on the distance between particles, as well as on arbitrary global and
    per-particle parameters.  It also optionally supports periodic boundary conditions and cutoffs for long range interactions.

    To use this class, create a CustomNonbondedForce object, passing an algebraic expression to the constructor
    that defines the interaction energy between each pair of particles.  The expression may depend on r, the distance
    between the particles, as well as on any parameters you choose.  Then call addPerParticleParameter() to define per-particle
    parameters, and addGlobalParameter() to define global parameters.  The values of per-particle parameters are specified as
    part of the system definition, while values of global parameters may be modified during a simulation by calling Context::setParameter().

    Next, call addParticle() once for each particle in the System to set the values of its per-particle parameters.
    The number of particles for which you set parameters must be exactly equal to the number of particles in the
    System, or else an exception will be thrown when you try to create a Context.  After a particle has been added,
    you can modify its parameters by calling setParticleParameters().

    CustomNonbondedForce also lets you specify 'exclusions', particular pairs of particles whose interactions should be
    omitted from force and energy calculations.  This is most often used for particles that are bonded to each other.

    As an example, the following code creates a CustomNonbondedForce that implements a 12-6 Lennard-Jones potential:

    >>> force = CustomNonbondedForce('4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)')

    This force depends on two parameters: sigma and epsilon.  The following code defines these parameters, and
    specifies combining rules for them which correspond to the standard Lorentz-Bertelot combining rules:

    >>> force.addPerParticleParameter('sigma')
    0
    >>> force.addPerParticleParameter('epsilon')
    1

    Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
    functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, sinh, cosh, tanh.  All trigonometric functions
    are defined in radians, and log is the natural logarithm.  The names of per-particle parameters have the suffix '1' or '2'
    appended to them to indicate the values for the two interacting particles.  As seen in the above example, the expression may
    also involve intermediate quantities that are defined following the main expression, using ';' as a separator.

    In addition, you can call addFunction() to define a new function based on tabulated values.  You specify a vector of
    values, and an interpolating or approximating spline is created from them.  That function can then appear in the expression.

    EXAMPLES


    """

    NoCutoff = 0 #: No cutoff is applied to nonbonded interactions.  The full set of N^2 interactions is computed exactly. This necessarily means that periodic boundary conditions cannot be used.  This is the default.
    CutoffNonPeriodic = 1 #: Interactions beyond the cutoff distance are ignored.
    CutoffPeriodic = 2 #: Periodic boundary conditions are used, so that each particle interacts only with the nearest periodic copy of each other particle.  Interactions beyond the cutoff distance are ignored.

    def __init__(self, energy=None, parameter_units=None, force=None):
        """
        Create a CustomNonbondedForce.

        @param energy    an algebraic expression giving the interaction energy between two particles as a function of r, the distance between them, as well as any global and per-particle parameters



        """

        # Initialize with defaults.
        self.energyExpression = energy #: Energy expression
        self.cutoffDistance = units.Quantity(1.0, units.nanometer)
        self.nonbondedMethod = CustomNonbondedForce.NoCutoff

        # Ensure that either the energy or the force object is specified.
        if (energy is None) and (force is None):
            raise ValueException("An energy function expression must be specified.")

        self.parameters = list() #: parameters[i] is the ParticleInfo object for particle i
        self.globalParameters = list() #: globalParameters[i] is the GlobalParameterInfo object for particle i
        self.particles = list()
        self.exclusions = list()
        self.functions = list()

        # Populate data structures from swig object, if specified
        if force is not None:
            self._copyDataUsingInterface(self, force)

        # If parameters are specified, check that the expression is correct.
        if (parameter_units is not None):
            self.checkUnits(energy, parameter_units)

        return

    def _checkUnits(self, energy_expression, parameter_units):
        """
        Check to ensure that the specified energy_expression and parameter_units have compatible dimensions.

        @param energy_expression    a string specifying the energy expression for Lepton
        @param parameter_units      a dict specifying the units for each parameter (e.g. parameter_units['epsilon'] = 'kilocalories_per_mole')

        """

        # TODO

        return

    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        dest.setNonbondedMethod( src.getNonbondedMethod() )
        dest.setCutoffDistance( src.getCutoffDistance() )
        for index in range(src.getNumPerParticleParameters()):
            name = src.getPerParticleParameterName(index)
            dest.addPerParticleParameter(name)
        for index in range(src.getNumGlobalParameters()):
            name = src.getGlobalParameterName(index)
            defaultValue = src.getGlobalParameterDefaultValue(index)
            dest.addGlobalParameter(name, defaultValue)
        for index in range(src.getNumFunctions()):
            args = src.getFunctionParameters(index)
            dest.addFunction(*args)
        for index in range(src.getNumParticles()):
            parameters = src.getParticleParameters(index)
            dest.addParticle(*parameters)
        for index in range(src.getNumExclusions()):
            (particle1, particle2) = src.getExclusionParticles(index)
            dest.addExclusion(particle1, particle2)

        return

    def asSwig(self):
        """
        Construct a Swig object.
        """
        force = openmm.CustomNonbondedForce(self.energyExpression)
        self._copyDataUsingInterface(force, self)
        return force

    def getNumParticles(self):
        """
        Get the number of particles for which force field parameters have been defined.

        """
        return len(self.particles)

    def getNumExclusions(self):
        """
        Get the number of particle pairs whose interactions should be excluded. 

        """
        return len(self.exclusions)

    def getNumPerParticleParameters(self):
        """
        Get the number of per-particle parameters that the interaction depends on.

        """
        return len(self.parameters)

    def getNumGlobalParameters(self):
        """
        Get the number of global parameters that the interaction depends on.

        """
        return len(self.globalParameters)

    def getNumFunctions(self):
        """
        Get the number of tabulated functions that have been defined.

        """
        return len(self.functions)

    def getEnergyFunction(self):
        """
        Get the number of special interactions that should be calculated differently from other interactions.

        """
        return self.energyExpression

    def setEnergyFunction(self, energy):
        """
        Set the algebraic expression that gives the interaction energy between two particles.

        @params energy  an algebraic expression giving the interaction energy between two particles as a function of r, the distance between them, as well as any global and per-particle parameters 

        """
        self.energyExpression = energy
        return

    def getNonbondedMethod(self):
        """
        Get the method used for handling long range nonbonded interactions.

        """
        return self.nonbondedMethod

    def setNonbondedMethod(self, nonbondedMethod):
        """
        Set the method used for handling long range nonbonded interactions.

        """
        # TODO: Argument checking.
        self.nonbondedMethod = nonbondedMethod
        return

    def getCutoffDistance(self):
        """
        Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
        is NoCutoff, this value will have no effect.

        """
        return self.cutoffDistance

    @accepts_compatible_units(units.nanometers)
    def setCutoffDistance(self, cutoffDistance):
        """
        Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
        is NoCutoff, this value will have no effect.

        """
        self.cutoffDistance = cutoffDistance
        return

    def addPerParticleParameter(self, name):
        """
        Add a new per-particle parameter that the interaction may depend on.

        @param name     the name of the parameter
        @return the index of the parameter that was added

        """
        self.parameters.append(name)
        return (len(self.parameters) - 1)

    def getPerParticleParameterName(self, index):
        """
        Get the name of a per-particle parameter.

        @param index     the index of the parameter for which to get the name
        @return the parameter name
        """
        return self.parameters[index]

    def setPerParticleParameterName(self, index, name):
        """
        Set the name of a per-particle parameter.

        @param index          the index of the parameter for which to set the name
        @param name           the name of the parameter

        """
        self.parameters[index] = name
        return

    def addGlobalParameter(self, name, defaultValue):
        """
        Add a new global parameter that the interaction may depend on.

        @param name             the name of the parameter
        @param defaultValue     the default value of the parameter
        @return the index of the parameter that was added

        """
        parameter = CustomNonbondedForceGlobalParameterInfo(name, defaultValue)
        self.globalParameters.append(parameter)
        return (len(self.globalParameters) - 1)

    def getGlobalParameterName(self, index):
        """
        Get the name of a global parameter.

        @param index     the index of the parameter for which to get the name
        @return the parameter name

        """
        return self.globalParameters[index].name

    def setGlobalParameterName(self, index, name):
        """
        Set the name of a global parameter.

        @param index          the index of the parameter for which to set the name
        @param name           the name of the parameter

        """
        self.globalParameters[index].name = name
        return

    def getGlobalParameterDefaultValue(self, index):
        """
        Get the default value of a global parameter.

        @param index     the index of the parameter for which to get the default value
        @return the parameter default value

        """
        return self.globalParameters[index].defaultValue

    def setGlobalParameterDefaultValue(self, index, defaultValue):
        """
        Set the default value of a global parameter.

        @param index          the index of the parameter for which to set the default value
        @param name           the default value of the parameter

        """
        self.globalParameters[index].defaultValue = defaultValue
        return

    def addParticle(self, *parameters):
        """
        Add the nonbonded force parameters for a particle.  This should be called once for each particle
        in the System.  When it is called for the i'th time, it specifies the parameters for the i'th particle.

        @param parameters    the list of parameters for the new particle
        @return the index of the particle that was added

        """
        # Unpack until we have just a list or a tuple of parameters.
        while (len(parameters) == 1) and type(parameters[0]) in (list, tuple):
            parameters = parameters[0]
        particle = CustomNonbondedForceParticleInfo(parameters)
        self.particles.append(particle)
        return (len(self.particles) - 1)

    def getParticleParameters(self, index):
        """
        Get the nonbonded force parameters for a particle.

        @param index       the index of the particle for which to get parameters
        @param parameters  the list of parameters for the specified particle

        """
        particle = self.particles[index]
        return particle.parameters

    def setParticleParameters(self, index, *parameters):
        """
        Set the nonbonded force parameters for a particle.

        @param index       the index of the particle for which to set parameters
        @param parameters  the list of parameters for the specified particle 

        """
        # Unpack until we have just a list or a tuple of parameters.
        while (len(parameters) == 1) and type(parameters[0]) in (list, tuple):
            parameters = parameters[0]
        particle = CustomNonbondedForceParticleInfo(parameters)
        self.particles[index] = particle
        return

    def addExclusion(self, particle1, particle2):
        """
        Add a particle pair to the list of interactions that should be excluded.

        @param particle1  the index of the first particle in the pair
        @param particle2  the index of the second particle in the pair
        @return the index of the exclusion that was added

        """
        exclusion = CustomNonbondedForceExclusionInfo(particle1, particle2)
        self.exclusions.append(exclusion)
        return (len(self.exclusions) - 1)

    def getExclusionParticles(self, index):
        """
        Get the particles in a pair whose interaction should be excluded.

        @param index      the index of the exclusion for which to get particle indices
        @return particle1  the index of the first particle in the pair
        @return particle2  the index of the second particle in the pair

        """
        exclusion = self.exclusions[index]
        return (exclusion.particle1, exclusion.particle2)

    def setExclusionParticles(self, index, particle1, particle2):
        """
        Set the particles in a pair whose interaction should be excluded.

        @param index      the index of the exclusion for which to set particle indices
        @param particle1  the index of the first particle in the pair
        @param particle2  the index of the second particle in the pair

        """
        self.exclusions[index] = CustomNonbondedForceExclusionInfo(particle1, particle2)
        return

    def addFunction(self, name, values, min, max, interpolating):
        """
        Add a tabulated function that may appear in the energy expression.

        @param name           the name of the function as it appears in expressions
        @param values         the tabulated values of the function f(x) at uniformly spaced values of x between min and max.
                              The function is assumed to be zero for x &lt; min or x &gt; max.
        @param min            the value of the independent variable corresponding to the first element of values
        @param max            the value of the independent variable corresponding to the last element of values
        @param interpolating  if true, an interpolating (Catmull-Rom) spline will be used to represent the function.
                              If false, an approximating spline (B-spline) will be used.
        @return the index of the function that was added

        """
        function = CustomNonbondedForceFunctionInfo(name, values, min, max, interpolating)
        self.functions.append(function)
        return (len(self.functions) - 1)

    def getFunctionParameters(self, index):
        """
        Get the parameters for a tabulated function that may appear in the energy expression.

        @param index          the index of the function for which to get parameters
        @return name          the name of the function as it appears in expressions
        @return values        the tabulated values of the function f(x) at uniformly spaced values of x between min and max.
                              The function is assumed to be zero for x &lt; min or x &gt; max.
        @return min           the value of the independent variable corresponding to the first element of values
        @return max           the value of the independent variable corresponding to the last element of values
        @return interpolating if true, an interpolating (Catmull-Rom) spline will be used to represent the function.
                              If false, an approximating spline (B-spline) will be used.

        """
        function = self.functions[index]
        return (function.name, function.values, function.min, function.max, function.interpolating)

    def setFunctionParameters(self, index, name, values, min, max, interpolating):
        """
        Set the parameters for a tabulated function that may appear in algebraic expressions.

        @param index          the index of the function for which to set parameters
        @param name           the name of the function as it appears in expressions
        @param values         the tabulated values of the function f(x) at uniformly spaced values of x between min and max.
                              The function is assumed to be zero for x &lt; min or x &gt; max.
        @param min            the value of the independent variable corresponding to the first element of values
        @param max            the value of the independent variable corresponding to the last element of values
        @param interpolating  if true, an interpolating (Catmull-Rom) spline will be used to represent the function.

        """
        function = CustomNonbondedForceFunctionInfo(name, values, min, max, interpolating)
        self.functions[index] = function
        return
    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def nparameters(self):
        return len(self.parameters)

    @property
    def nglobalparameters(self):
        return len(self.globalParameters)

    @property
    def nparticles(self):
        return len(self.particles)

    @property
    def nexclusions(self):
        return len(self.exclusions)

    @property
    def nfunctions(self):
        return len(self.functions)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the CustomNonbondedForce object.

        """

        r = ""

        # Show settings.
        r += "energyExpression: %s" % str(self.energyExpression)

        # TODO: Show nonbondedmethod, cutoff, particles, functions, exceptions, parameters, etc.

        return r

    def _appendForce(self, force, offset):
        """
        Append atoms defined in another force of the same type.

        @param force      the force to be appended
        @param offset     integral offset for atom number

        EXAMPLES

        Create two force objects and append the second to the first.

        >>> energy = 'k*r^2'
        >>> force1 = CustomNonbondedForce(energy)
        >>> force1.addGlobalParameter('k', 1.0 * units.kilocalories_per_mole / units.nanometers)
        0
        >>> force1.addParticle()
        0
        >>> force2 = CustomNonbondedForce(energy)
        >>> force2.addGlobalParameter('k', 1.0 * units.kilocalories_per_mole / units.nanometers)
        0
        >>> force2.addParticle()
        0
        >>> offset = 1
        >>> force1._appendForce(force2, offset)

        """

        # TODO: Allow coercion of Swig force objects into Python force objects.

        # Check to make sure both forces are compatible.
        if type(self) != type(force):
            raise ValueError("other force object must be of identical Force subclass")

        # Check to make sure force objects are compatible.
        if (self.energyExpression != force.energyExpression):
            raise ValueError("other force has incompatible energyExpression")
        if (self.nonbondedMethod != force.nonbondedMethod):
            raise ValueError("other force has incompatible nonbondedMethod")
        if (self.cutoffDistance != force.cutoffDistance):
            raise ValueError("other force has incompatible cutoffDistance")
        if (self.nfunctions != force.nfunctions):
            raise ValueError("other force has incompatible functions")
        for (function1, function2) in zip(self.functions, force.functions):
            if (function1 != function2):
                raise ValueError("functions are incompatible: (%s) and (%s)" % (function1, function2))
        # TODO: More checking?

        # Combine systems.
        for particle in force.particles:
            self.particles.append(particle)
        for exclusion in force.exclusions:
            offset_exclusion = CustomNonbondedForceExclusionInfo(exclusion.particle1 + offset, exclusion.particle2 + offset)
            self.exclusions.append(offset_exclusion)

        return


class CustomNonbondedForceParticleInfo(object):
    def __init__(self, *parameters):
        self.parameters = parameters
        return

class CustomNonbondedForcePerParticleParameterInfo(object):
    def __init__(self, name):
        self.name = name
        return

class CustomNonbondedForceGlobalParameterInfo(object):
    def __init__(self, name, defaultValue):
        self.name = name
        self.defaultValue = defaultValue
        return

class CustomNonbondedForceExclusionInfo(object):
    def __init__(self, particle1, particle2):
        self.particle1 = particle1
        self.particle2 = particle2
        return

class CustomNonbondedForceFunctionInfo(object):
    def __init__(self, name, values, min, max, interpolating):
        self.name = name
        # TODO: Check to ensure 'values' is numpy array or list.
        self.values = values
        self.min = min
        self.max = max
        self.interpolating = interpolating
        return

    def __eq__(self, other):
        if ((self.name != other.name) or
            (self.values != other.values) or
            (self.min != other.min) or
            (self.max != other.max) or
            (self.interpolating != other.interpolating)): 
            return false
        return true

#=============================================================================================
# CustomExternalForce
#=============================================================================================

class CustomExternalForce(Force):
    """
    This class implements an 'external' force on particles.  The force may be applied to any subset of the particles
    in the System.  The force on each particle is specified by an arbitrary algebraic expression, which may depend
    on the current position of the particle as well as on arbitrary global and per-particle parameters.

    To use this class, create a CustomExternalForce object, passing an algebraic expression to the constructor
    that defines the potential energy of each affected particle.  The expression may depend on the particle's x, y, and
    z coordinates, as well as on any parameters you choose.  Then call addPerParticleParameter() to define per-particle
    parameters, and addGlobalParameter() to define global parameters.  The values of per-particle parameters are specified as
    part of the system definition, while values of global parameters may be modified during a simulation by calling Context::setParameter().
    Finally, call addParticle() once for each particle that should be affected by the force.  After a particle has been added,
    you can modify its parameters by calling setParticleParameters().

    As an example, the following code creates a CustomExternalForce that attracts each particle to a target position (x0, y0, z0)
    via a harmonic potential:

    >>> force = CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2')

    This force depends on four parameters: the spring constant k and equilibrium coordinates x0, y0, and z0.  The following code defines these parameters:

    >>> force.addGlobalParameter('k', 1.0)
    0
    >>> force.addPerParticleParameter('x0')
    0
    >>> force.addPerParticleParameter('y0')
    1
    >>> force.addPerParticleParameter('z0')
    2

    Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
    functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, sinh, cosh, tanh, step.  All trigonometric functions
    are defined in radians, and log is the natural logarithm.  step(x) = 0 if x is less than 0, 1 otherwise.

    """
    def __init__(self, energy):
        """
        Create a CustomExternalForce.

        @param energy    an algebraic expression giving the potential energy of each particle as a function
                         of its x, y, and z coordinates

        """
        self.energyExpression = energy #: The algebraic expression giving the energy for each particle
        self.parameters = list() #: parameters[i] is the ith per-particle PerParticleParameterInfo object
        self.globalParameters = list() #: globalParameters[i] is the ith GlobalParameterInfo object
        self.particles = list() #: particles[i] is the ith ParticleInfo object containing per-particle parameters for particle i
        return


    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        #dest.__init__( src.getEnergyFunction() ) # TODO: Do we need this?
        for index in range(src.getNumPerParticleParameters()):
            name = src.getPerParticleParameterName(index)
            dest.addPerParticleParameter(name)
        for index in range(src.getNumGlobalParameters()):
            name = src.getGlobalParameterName(index)
            defaultValue = src.getGlobalParameterDefaultValue(index)
            dest.addGlobalParameter(name, defaultValue)
        for index in range(src.getNumParticles()):
            (particle, parameters) = src.getParticleParameters(index)
            dest.addParticle(particle, parameters)

        return

    def asSwig(self):
        """
        Construct a Swig object.
        """
        force = openmm.CustomExternalForce(self.energyExpression)
        self._copyDataUsingInterface(force, self)
        return force

    def getNumParticles(self):
        """
        Get the number of particles for which force field parameters have been defined.

        @returns the number of particles for which force field parameters have been defined

        """
        return len(self.particles)

    def getNumPerParticleParameters(self):
        """
        Get the number of per-particle parameters that the force depends on.

        @returns the number of per-particle parameters that have been specified

        """
        return len(self.parameters)

    def getNumGlobalParameters(self):
        """
        Get the number of global parameters that the force depends on.

        @returns the number of global parameters that have been specified

        """

        return len(self.globalParameters)

    def getEnergyFunction(self):
        """
        Get the algebraic expression that gives the potential energy of each particle

        @returns the algebraic energy expression

        """
        return self.energyExpression

    def setEnergyFunction(self, energy):
        """
        Set the algebraic expression that gives the potential energy of each particle

        @param energy   an algebraic expression specifying the potential energy of each particle

        """
        # TODO: Argument checking with Lepton.
        self.energyExpression = energy

    def addPerParticleParameter(self, name):
        """
        Add a new per-particle parameter that the force may depend on.

        @param name             the name of the parameter
        @return the index of the parameter that was added

        """
        parameter = CustomExternalForcePerParticleParameterInfo(name)
        self.parameters.append(parameter)
        return (len(self.parameters) - 1)

    def getPerParticleParameterName(self, index):
        """
        Get the name of a per-particle parameter.

        @param index     the index of the parameter for which to get the name
        @return the parameter name

        """
        return self.parameters[index].name

    def setPerParticleParameterName(self, index, name):
        """
        Set the name of a per-particle parameter.

        @param index          the index of the parameter for which to set the name
        @param name           the name of the parameter

        """
        self.parameters[index].name = name

    def addGlobalParameter(self, name, defaultValue):
        """
        Add a new global parameter that the force may depend on.

        @param name             the name of the parameter
        @param defaultValue     the default value of the parameter
        @return the index of the parameter that was added

        """
        globalParameter = CustomExternalForceGlobalParameterInfo(name, defaultValue)
        self.globalParameters.append(globalParameter)
        return (len(self.globalParameters) - 1)

    def getGlobalParameterName(self, index):
        """
        Get the name of a global parameter.

        @param index     the index of the parameter for which to get the name
        @return the parameter name

        """
        return self.globalParameters[index].name

    def setGlobalParameterName(self, index, name):
        """
        Set the name of a global parameter.

        @param index          the index of the parameter for which to set the name
        @param name           the name of the parameter

        """
        self.globalParameters[index].name = name

    def getGlobalParameterDefaultValue(self, index):
        """
        Get the default value of a global parameter.

        @param index     the index of the parameter for which to get the default value
        @return the parameter default value

        """
        return self.globalParameters[index].defaultValue

    def setGlobalParameterDefaultValue(self, index, defaultValue):
        """
        Set the default value of a global parameter.

        @param index          the index of the parameter for which to set the default value
        @param name           the default value of the parameter

        """
        self.globalParameters[index].defaultValue = defaultValue

    def addParticle(self, particle, *parameters):
        """
        Add a particle term to the force field.

        @param particle     the index of the particle this term is applied to
        @param parameters   the list of parameters for the new force term
        @return the index of the particle term that was added

        """
        # Unpack until we have just a list or a tuple of parameters.
        while (len(parameters) == 1) and type(parameters[0]) in (list, tuple):
            parameters = parameters[0]
        particle = CustomExternalForceParticleInfo(particle, parameters)
        self.particles.append(particle)
        return (len(self.particles) - 1)

    def getParticleParameters(self, index):
        """
        Get the force field parameters for a force field term.

        @param index         the index of the particle term for which to get parameters
        @returns particle      the index of the particle this term is applied to
        @returns parameters    the list of parameters for the force field term

        """
        particleinfo = self.particles[index]
        return (particleinfo.particle, particleinfo.parameters)

    def setParticleParameters(self, index, particle, *parameters):
        """
        Set the force field parameters for a force field term.

        @param index         the index of the particle term for which to set parameters
        @param particle      the index of the particle this term is applied to
        @param parameters    the list of parameters for the force field term

        """
        # Unpack until we have just a list or a tuple of parameters.
        while (len(parameters) == 1) and type(parameters[0]) in (list, tuple):
            parameters = parameters[0]
        self.particles[index] = self.ParticleInfo(particle, parameters)

    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def nparameters(self):
        return len(self.parameters)

    @property
    def nglobalparameters(self):
        return len(self.globalParameters)

    @property
    def nparticles(self):
        return len(self.particles)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the CustomNonbondedForce object.

        """

        r = ""

        # Show settings.
        r += "energyExpression: %s\n" % str(self.energyExpression)

        for (index, particleinfo) in enumerate(self.particles):
            r += "%8d %8d %s\n" % (index, particleinfo.particle, str(particleinfo.parameters))

        # TODO

        return r

    def _appendForce(self, force, offset):
        """
        Append atoms defined in another force of the same type.

        @param force      the force to be appended
        @param offset     integral offset for atom number

        EXAMPLES

        Create two force objects and append the second to the first.

        >>> energy = 'Ez*q*z'
        >>> force1 = CustomExternalForce(energy)
        >>> force1.addGlobalParameter('Ez', 1.0 * units.kilocalories_per_mole / units.nanometers / units.elementary_charge)
        0
        >>> force1.addPerParticleParameter('q')
        0
        >>> force1.addParticle(0, 1.0 * units.elementary_charge)
        0
        >>> force2 = CustomExternalForce(energy)
        >>> force2.addGlobalParameter('Ez', 1.0 * units.kilocalories_per_mole / units.nanometers / units.elementary_charge)
        0
        >>> force2.addPerParticleParameter('q')
        0
        >>> force2.addParticle(0, 1.0 * units.elementary_charge)
        0
        >>> offset = 1
        >>> force1._appendForce(force2, offset)

        """

        # TODO: Allow coercion of Swig force objects into Python force objects.

        # Check to make sure both forces are compatible.
        if type(self) != type(force):
            raise ValueError("other force object must be of identical Force subclass")

        # Check to make sure force objects are compatible.
        if (self.energyExpression != force.energyExpression):
            raise ValueError("other force has incompatible energyExpression")
        # TODO: More checking?

        # Combine systems.
        for particleinfo in force.particles:
            modified_particleinfo = CustomExternalForceParticleInfo(particleinfo.particle + offset, particleinfo.parameters)
            self.particles.append(modified_particleinfo)

        return

class CustomExternalForceParticleInfo(object):
    def __init__(self, particle, parameters):
        self.particle = particle
        self.parameters = parameters
        return

class CustomExternalForcePerParticleParameterInfo(object):
    def __init__(self, name):
        self.name = name
        return

class CustomExternalForceGlobalParameterInfo(object):
    def __init__(self, name, defaultValue):
        self.name = name
        self.defaultValue = defaultValue
        return
