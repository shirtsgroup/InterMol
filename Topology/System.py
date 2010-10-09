#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
System.py

Classes and methods for creating and manipulating molecular system topologies.

REQUIREMENTS

The SimTK python_units package must be installed See: https://simtk.org/home/python_units




HISTORY

This is based *heavily* on John Chodera's
Pure Python pyopenmm.py code, which is in turn based heavily on Chris Bruns' SimTK PyOpenMM code.

COPYRIGHT

@author Vincent Voelz <vvoelz@gmail.com>
@author John D. Chodera <jchodera@gmail.com>



EXAMPLES



TODO


"""


#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import re, os, sys, pdb
import copy
import numpy
from Force import *
from Structure import *

defaultdictionary={
"nbfunc":"NBFunc",
"cr":"CombinationRule",
"fudgeQQ":"CoulombCorrection",
"fudgeLJ":"LJCorrection"
}



#=============================================================================================
# System base class
#=============================================================================================

class System(object):
    """
    This class represents a molecular system topology.

    The System object is the outermost container, and stores a list of Topology objects along with  corresponding structure objects

    """

    def __init__(self, name = None, structureList=None, topologyList=None, inputFileType=None):
        """
        Create a new System object.

        If a System object is specified, it will be queried to construct the class.

        """

        self.name = name
        if self.name == None:
            self.name = "Untitled"

        self.structures = list()
        self.structureMolecules = list()
        self.topologyMolecules = list()

        self.loadStructure(structureList, inputFileType)
        self.loadTopology(topologyList, inputFileType)

        self.setSystemParameters()
        self.mapping = Mapping(self.structureMolecules, self.topologyMolecules)


    def loadStructure(self, structureFiles, inputFileType):
        """
        Load a structureFile into the system

        """

        if type(structureFiles).__name__ == 'str':
            temp = structureFiles
            structureFiles = list()
            structureFiles.append(temp)

        if inputFileType == "Gromacs":
            from GromacsTopology import *
            for structure in structureFiles:
                # Create GromacsStructure/readStructureFile
                gs = GromacsStructure(grofile=structure)
                gs.readStructureFile(structure)
                # Create generic Structure/convertFromGromacsStructure
                s = Structure()
                s.convertFromGromacsStructure(gs)
                self.structureMolecules.extend(s.molecules)
                self.structures.append(s)

        if inputFileType == "Amber":
            # Import Amber stuff?
            for structure in structureFiles:
                # Create AmberStructure/readStructureFile
                gs = AmberStructure(grofile=structure)
                gs.readStructureFile(structure)
                # Create generic Structure/convertFromGromacsStructure
                s = Structure()
                s.convertFromAmberStructure(gs)
                self.structureMolecules.extend(s.molecules)
                self.structures.append(s)

        if inputFileType == "Desmond":
            # Import Desmond
            for structure in structureFiles:
                # Create DesmondStructure/readStructureFile
                gs = DesmondStructure(grofile=structure)
                gs.readStructureFile(structure)
                # Create generic Structure/convertFromGromacsStructure
                s = Structure()
                s.convertFromDesmondStructure(gs)
                self.structureMolecules.extend(s.molecules)
                self.structures.append(s)


    def loadTopology(self, topologyFiles, inputFileType):
        """
        Load a topologyFile into the system

        """

        if type(topologyFiles).__name__ == 'str':
            temp = topologyFiles
            topologyFiles = list()
            topologyFiles.append(temp)

        if inputFileType == "Gromacs":
            from GromacsTopology import *
            for topology in topologyFiles:
                # Create object/readTopologyFile
                g = GromacsTopology(topfile=topology)
                self.topologyMolecules.extend(g.molecules)

        if inputFileType == "Amber":
            # Import Amber
            for topology in topologyFiles:
                # Create object/readTopologyFile
                g = AmberTopology(topfile=topology)
                self.topologyMolecules.extend(g.molecules)

        if inputFileType == "Desmond":
            # Import Desmond
            for topology in topologyFiles:
                # Create object/readTopologyFile
                g = DesmondTopology(topfile=topology)
                self.topologyMolecules.extend(g.molecules)


    def setSystemParameters(self):
        # From Structure
        if len(self.structures) == 1:
            self.name = self.structures[0].title
            self.boxvector = self.structures[0].boxvector
        else:
            # Need a method of determining which name and boxvector to use
            self.name = ''
            self.boxvector = ''

        # From Topology
        for topology in self.topologyMolecules:
	    # Gromacs defaults
		if len(topology.defaults) == 5:
                	self.NBFunc = topology.defaults[0]
                	self.CombinationRule = topology.defaults[1]
               		self.LJCorrection = topology.defaults[3]
                	self.CoulombCorrection = topology.defaults[4]
        return


    def deleteStructure(self, title):
        for structure in self.structures:
            if title == structure.title:
                self.structures.remove(structure)
        return


    def deleteTopology(self, name):
        for topology in topologyMolecules:
            if name == topology.name:
                self.topologyMolecules.remove(topology)
        return


    def getName(self):
        """
        Get the name of the System.

        """
        return self.name


    def setName(self, name):
        """
        Set the name of the System.

        """
        self.name = name
        return


    def getNumMolecules(self):
        """
        Get the number of molecules in the System.

        """
        return len(self.structureMolecules)


    def writeToFiles(self, printFileType, topologyName='system.top', structureName='system.gro', molecules='All'):
        """
        Write the system's contents out to useable files in the desired format.

        """

        # Change inputs to lists
        # If list then go through a list like in writeStructureFile

        if printFileType == "Gromacs":
            from GromacsTopology import *

            # Topology
            if molecules == 'All':
                gromacsTopology = GromacsTopology(self.topologyMolecules)
            else:
                # Function to connect input name to actual topology molecules
                gromacsTopology = GromacsTopology(molecules)
            gromacsTopology.writeTopologyFile(topologyName, ExpandIncludes=False, RebuildDirectives=True)

            # Structure
            if molecules == 'All':
                gromacsStructure = GromacsStructure(self.structures)
            else:
                # Function to connect input name to actual structure molecules
                gromacsStructure = GromacsStructure(molecules)
            gromacsStructure.writeStructureFile(structureName, RebuildStructure=True)

        if printFileType == "Amber":
            pass

        if printFileType == "Desmond":
            pass

        return


    #def addMolecule(self, molecule=None):
        #"""
        ##Add a molecule (i.e. a Topology object) to the System object.

        #"""
        #if molecule == None:
            #self.molecules.append( Topology() )
        #else:
            #self.molecules.append( Topology( system=molecule) )
        #return


    #def delMolecule(self, index):
        #"""
        ##Delete a molecule from the Topology object at the specified index.

        #"""
        #self.molecules.pop(index)
        #return


    #def insertMolecule(self, index, molecule=None):
        #"""
        ##Insert a molecule into the Topology object at the specified index.

        #"""
        #if molecule == None:
            #self.molecules.insert( index, TopologySystem() )
        #else:
            #self.molecules.insert( index, TopologySystem( system=molecule) )
        #return



class Mapping(object):


    def __init__(self, structureMolecules, topologyMolecules):
        """
        Set up 'proteinMap' with the following form for each entry:

        {key:list()}

        where - 'key' is the name of a topology
        and   - 'list()' contains pairwise entries of the form ['molecule', 'dict{strucAtom:topAtom}']

        Example:

         TopMol      StrucMol          Atoms

        {Protein : (['Protein', {strucAtom:topAtom])
                                 strucAtom:topAtom
                                    ...   :  ...  }
         solute: (['MOL', {strucAtom:topAtom])
                           strucAtom:topAtom
                              ...   :  ...  }
         WAT: (['WAT1', {strucAtom:topAtom  ], ['WAT2', {strucAtom:topAtom  ], ['WAT3', {strucAtom:topAtom  ])}
                         strucAtom:topAtom               strucAtom:topAtom               strucAtom:topAtom
                            ...   :  ...  }                 ...   :  ...  }                 ...   :  ...  }
        """

        self.proteinMap = {}

        # for every 'key'
        for topology in topologyMolecules:
            molnameList = list()
            # find every corresponding 'molecule'
            for molecule in structureMolecules:
                if molecule.name == topology.name:
                    atomMap = {}
                    # for every 'topAtom' in 'key'
                    for topAtom in topology.atoms:
                        # find every corresponding 'strucAtom' in 'molecule'
                        for strucAtom in molecule.atomList:
                            if strucAtom.atomname == topAtom.atomname and strucAtom.resname == topAtom.metadata.resname:
                                # add {'strucAtom' : 'topAtom'} pairs into an 'atomMap' dictionary
                                atomMap[strucAtom.atomnum] = topAtom.ID
                    # append current ['molecule', 'dict{strucAtom:topAtom}'] entry to 'molnameList()'
                    molnameList.append([str(molecule.atomList[0].resnum) + molecule.name, atomMap])
            if topology.name not in self.proteinMap:
                self.proteinMap[topology.name] = molnameList


    def find_key(dic, val):
        """
        Return the key of dictionary dic given the value

        """

        return [k for k, v in symbol_dic.iteritems() if v == val][0]


class Topology(object):
    """
    A Topology object consists of:

        <ol>
        <li>The set of particles in the system</li>
        <li>The forces acting on them</li>
        <li>Pairs of particles whose separation should be constrained to a fixed value</li>
        <li>For periodic systems, the dimensions of the periodic box</li>
        </ol>

        The particles and constraints are defined directly by the Topology object, while
        forces are defined by objects that extend the Force class.  After creating a
        Topology, call addParticle() once for each particle, addConstraint() for each constraint,
        and addForce() for each Force.


    Examples:

    >>> topology = Topology()

    Add a particle.

    >>> mass = 12.0 * units.amu
    >>> topology.addParticle(mass)
    0

    Add a NonbondedForce.

    >>> nonbondedForce = NonbondedForce()
    >>> topology.addForce(nonbondedForce)
    0

    Create a deep copy.

    >>> import copy
    >>> topology_copy = copy.deepcopy(topology)


    * A pyopenmm.System-like object stores the atoms in each molecule, their connectivity and forces:

      NOTE that unlike PyOpenMM's System.py, extra metadata information is stored here:
      Atom names, Atom types, along with functions to search through them.



    *** USAGE EXAMPLES ***


    Topology class should have methods for:
- insert and delete atoms
- insert and delete connections between atoms

- insert and delete atom groups
- merge and split atom groups
 -adding/deleting/joining groups.




- methods for selecting atom groups
   - select by residue
   - select by atomic type
   - select by atom number
   - select by physical location



    VAV: see  GromacsTopology.py for these:
    - outputting a gromacs top + gro file
    - reading in a topology from a gromacs .itp/.top/.gro




- constructing relative free energy topologies given n starting topologies
- building restraints between atoms
- taking a topology and a pdb/gro/etc. and outputting a new topology with the coordinates from the pdb
- creating a 'blank' topology from a pdb or mol2 file (will this be possible?)
- reset arbitrary key pairs in the input information (i.e., we load in the information once from a template file, and then can edit it keyword by keyword as desired).
- setting up A and B states -- that is, a gromacs topology that has a starting state (A) and an ending state (B). (Sort of implied by the relative free energy First: mention above)

    *** THINGS THAT WERE ALREADY IMPLEMENTED pyopenmm version***

    - merging (and splitting?) multiple topologies
    - generating angles from connectivities
    - generating torsions from connectivities
    - generating 1,4's from connectivities
      *  VAV: Force.NonbondedForce.createExceptionsFromBonds()

    VAV:  THIS IS CONTAINED IN THE Force.NonbondedForceParticleInfo class:
    - nonbonded parameters for all atoms
       - VdW parameters:
             epsilon (float) (kJ/mol)
             sigma (float)
       - partial charge (float)
       - mass (float)

    VAV:  THIS IS CONTAINED IN THE MetadataInfo class:
         - involved in free energy calculations (bool)
         - atom name (string)
         - residue name (string)

    VAV:  THIS IS CONTAINED IN THE Force.HarmonicBondForceBondInfo class:
     -list of bonds
         - list of parameters: two possibilities -- just a list that can be interpreted by whatever code we want, or exactly what's in the gromacs records

    VAV:  THIS IS CONTAINED IN THE Force.HarmonicAngleForceInfo class:
     -list of angles
         - list of parameters: see above

    VAV:  THIS IS CONTAINED IN THE Force.PeriodicTorsionForcePeriodicTorsionInfo class:
     -list of torsions
         - list of parameters: see above

    VAV:  THIS IS CONTAINED IN THE Force.NonbondedForceExceptionInfo class:
    -list of 1,4's
         - list of parameters: see above

    VAV:  THIS IS implemented via the AtomGroup class:
   - list of groups of atoms
      Groups are lists of atom indices, and a name string.

    VAV: THIS IS implemented via self.molecules()
      -can represent either single molecules, or systems of molecules

        Q: Can we allow for arbitrary per-atom data to be added later, such as implicit solvent parameters?
           You can either include this information with each atom, or you could think of including all GB data
           as a separate block of information separate from the nonbonded atom terms.I favor the OpenMM strategy
           of having a separate 'Force' object for each force component, so that NonbondedForce terms
           (charges, LJ sigma and epsilon) and metadata on how these nonbonded interactions are computed is separated
           from, say, a generalized Born section which also may have atomic parameters.  The reason for this is that
           the nonbonded parameters will also need associated with them a list of exclusions or exceptions for those
           interactions where the pairwise combining rule needs to be overridden. -John Chodera 5/26/10 7:47 PM

        VAV sez: For GB parameters, we can use the GB 'Force' objects provided in Force.py.  Otherwise, you can always add new
        attributes on the fly to the MetadataInfo() objects


        Lists of "molecules"?

        VAV sez:  One requirement was to store a list of "molecules" or objects (since an object may contain multiple molecules
        in the sense of containing multiple disconnected objects); each containing a list of atoms.  In this case, each "molecule"
        should be a separate Topology instance.  The GromacsTopology (and others, potentially) can serve as containers to store
        multiple "molecules"


        Combination rules to get vdW parameters for unspecified pairs?

        MRS sez:  another desired feature is to track nonbonded parameters for specific pairs of atoms to override the above
        JDC sez:  Agreed.  There needs to be the ability to exclude 1,2 and 1,3 pairs and attenuate 1,4 pairs, at the very least.
                  The simplest way is to just include a list of 'exceptions' for nonbonded interactions. -John Chodera 5/26/10 7:49 PM
        VAV:  This can be implemented via the Force.NonbondedForceExceptionInfo classes and associated methods.


    *** TO DO ***

    VAV: A task for MRS's students?
    - select by SMARTS string
     - Hand off the SMARTS string to open eye -- sdf/mol2 file


    VAV:  These are NOT implemented yet:  They would require their own derived Topology classes (like GromacsTopology)
    - reading in a topology from a hetgrpffgen file + pdb (Schrodinger OPLS-AA ?)  VAV: SchrodingerTopology.py ?
    - reading in a topology from Desmond templates(?)  VAV: DesmondTopology.py ?
    - reading in a topology from amber topology files (not necessary yet because of acpypi?) VAV: AmberTopology.py ?



    *** DISCUSSION - SHOULD THESE TO BE IMPLEMENTED? *** 

    -location: 3d coordinates: triplet of float (nm).
    -velocities?: triplet of float

    - includes the information necessary for the .gro and the .top in gromacs. ???

    - list of dummy atoms and parameters for construction of dummy atoms ???



    *** PURPOSELY NOT IMPLEMENTED ***

    MRS wrote:  "important note for free energy calculations:Each atom/angle/torsion should be able to have an arbitrary
                number of parameter sets!  For now, it will usually be one or two.  But eventually, we'd like to be able
                to support more.  So it should be a list of lists.

    VAV: This can be done for now using multiple objects with different sets of parms.

    JDC says:  "A critical question is how alchemical permutations are to be handled by this topology class.
                Should one topology object contain just a single topology, or should there be multiple sets
                of parameters for part or all of the molecule?  I favor just having *one* set of parameters
                defined, and using transformers that take one or two topology objects with identical ordering
                to generate alchemical transformation input files for gromacs, YANK, etc." -John Chodera 5/27/10 6:14 PM

    Q: Should we potentially include the runfile information, perhaps as a series of string pairs?
       Translators between formats might come later.
    A: I don't think we want runfile info mixed with the topology object. -David Mobley 5/25/10 5:29 PM



    *** TO DO ***

    - Very important:  need to lay out what changes need to be made to other files in MM tools
                   to be able to adapt to the new topology definition:


    """

    def __init__(self, Topology=None):
        """
        Create a new Topology.

        TODO:
        If a Topology object is specified, its attributes will be queried to construct the class.

        """

        self.name = '' # name of the topology
        self.defaults = list()
        self.nrexcl = ''
        self.atoms = list() # atoms[i] is the ith atom entry
        self.forces = list() # forces[i] is the ith force term
        self.constraints = list() # constraints[i] is the ith ConstraintInfo entry

        self.atomgroups = list() # atomgroups[i] is the ith AtomGroup object


        #self.periodicBoxVectors = [ units.Quantity((2.,0.,0.), units.nanometer), units.Quantity((0.,2.,0.), units.nanometer), units.Quantity((0.,0.,2.), units.nanometer) ] # periodic box vectors (only for periodic systems)
        # TODO: Store periodicBoxVectors as units.Quantity(numpy.array([3,3], numpy.float64), units.nanometer)?

        # Populate the system from a provided system, if given.
        if Topology is not None:
            self._copyDataUsingInterface(self, Topology)

        return


    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.

        """
        #dest.__init__()
        for index in range(src.getNumAtoms()):
            atoms = src.getAtomParameters(index)
            dest.addAtom(mass)
        for index in range(src.getNumConstraints()):
            args = src.getConstraintParameters(index)
            dest.addConstraint(*args)
        for index in range(src.getNumForces()):
            force = src.getForce(index)
            dest.addForce(force)
        box_vectors = src.getPeriodicBoxVectors()
        dest.setPeriodicBoxVectors(*box_vectors)

        return


    def getName(self):
        """
        Get the name of the Topology.

        """
        return self.name


    def setName(self, name):
        """
        Set the name of the Topology.

        """
        self.name = name
        return


    def getNumAtoms(self):
        """
        Get the number of atoms

        """
        return len(self.atoms)


    def addAtom(self, index):
        """
        Add a particle to the Topology.

        """
        return self.atoms[index]


    def delAtom(self, index, renumber=True):
        """
        Delete a particle from the Topology.
        Removes all references (forces, constraints, metadata), bonds) to this particle

        @return the index of the particle that was deleted

        """
        pass


    def insertAtom(self, index):
        """
        Insert a particle into the Topology at the specifed insertion index
        Renumbers all references (forces, constraints, metadata), bonds) to this particle accordingly

        @return the index of the particle that was inserted

        """
        pass


    def getAtomParameters(self, index):
        pass


    def getNumConstraints(self):
        """
        Get the number of distance constraints in this Topology.

        """
        return len(self.constraints)


    @accepts_compatible_units(None, None, units.nanometer)
    def addConstraint(self, particle1, particle2, distance):
        """
        Add a constraint to the Topology.

        @param particle1 the index of the first particle involved in the constraint
        @param particle2 the index of the second particle involved in the constraint
        @param distance  the required distance between the two particles, measured in nm
        @return the index of the constraint that was added

        """
        if particle1 not in range(self.getNumParticles()):
            raise ValueError("particle1 must be in range(0, getNumParticles())")
        if particle2 not in range(self.getNumParticles()):
            raise ValueError("particle1 must be in range(0, getNumParticles())")
        constraint = self.ConstraintInfo(particle1, particle2, distance)
        self.constraints.append(constraint)

        return


    def getConstraintParameters(self, index):
        """
        Get the parameters defining a distance constraint.

        @param index     the index of the constraint for which to get parameters
        @return a tuple of (particle1, particle2, distance) for the given constraint index

        """
        constraint = self.constraints[index]
        return (constraint.particle1, constraint.particle2, constraint.distance)


    @accepts_compatible_units(None, None, None, units.nanometer)
    def setConstraintParameters(self, index, particle1, particle2, distance):
        """
        Set the parameters defining a distance constraint.

        @param index     the index of the constraint for which to set parameters
        @param particle1 the index of the first particle involved in the constraint
        @param particle2 the index of the second particle involved in the constraint
        @param distance  the required distance between the two particles, measured in nm

        """
        if particle1 not in range(self.getNumParticles()):
            raise ValueError("particle1 must be in range(0, getNumParticles())")
        if particle2 not in range(self.getNumParticles()):
            raise ValueError("particle1 must be in range(0, getNumParticles())")
        constraint = self.ConstraintInfo(particle1,particle2,distance)
        self.constraints[index] = constraint

        return


    def addForce(self, force):
        """
        Add a Force to the Topology.

        @param force   the Force object to be added
        @return        the index within the System of the Force that was added

        NOTES

        """

        # Append the force.
        self.forces.append(force)
        return len(self.forces)-1


    def getNumForces(self):
        """
        Get the number of Force objects that have been added to the Topology.

        """
        return len(self.forces)


    def getForce(self, index):
        """
        Get a const reference to one of the Forces in this Topology.

        @param index  the index of the Force to get

        """
        return self.forces[index]


    def getPeriodicBoxVectors(self):
        """
        Get the vectors which define the axes of the periodic box (measured in nm).  These will affect
        any Force added to the System that uses periodic boundary conditions.

        Currently, only rectangular boxes are supported.  This means that a, b, and c must be aligned with the
        x, y, and z axes respectively.  Future releases may support arbitrary triclinic boxes.

        @returns a      the vector defining the first edge of the periodic box
        @returns b      the vector defining the second edge of the periodic box
        @returns c      the vector defining the third edge of the periodic box

        """
        return self.periodicBoxVectors


    def setPeriodicBoxVectors(self, a, b, c):
        """
        Set the vectors which define the axes of the periodic box (measured in nm).  These will affect
        any Force added to the System that uses periodic boundary conditions.

        Currently, only rectangular boxes are supported.  This means that a, b, and c must be aligned with the
        x, y, and z axes respectively.  Future releases may support arbitrary triclinic boxes.

        @param a      the vector defining the first edge of the periodic box
        @param b      the vector defining the second edge of the periodic box
        @param c      the vector defining the third edge of the periodic box

        """
        # TODO: Argument checking.
        self.periodicBoxVectors = [a,b,c]

    #==========================================================================
    # CONTAINERS
    #==========================================================================

    class ConstraintInfo(object):
        """
        Distance constraint information for particles in System object.

        """

        @accepts_compatible_units(None, None, units.nanometers)
        def __init__(self, particle1, particle2, distance):
            self.particle1 = particle1
            self.particle2 = particle2
            self.distance = distance
            return



    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================

    @property
    def nparticles(self):
        """
        The number of particles.

        """
        return len(self.masses)

    @property
    def nforces(self):
        """
        The number of force objects in the topology.

        """
        return len(self.forces)

    @property
    def nconstraints(self):
        """
        The number of interparticle distance constraints defined.

        """
        return len(self.constraints)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the Topology object.

        """

        r = ""
        r += "System object\n\n"

        # Show particles.
        r += "Particle masses:\n"
        r += "%8s %24s\n" % ("particle", "mass")
        for index in range(self.getNumParticles()):
            mass = self.getParticleMass(index)
            r += "%8d %24s\n" % (index, str(mass))
        r += "\n"

        # Show constraints.
        r += "Constraints:\n"
        r += "%8s %8s %16s\n" % ("particle1", "particle2", "distance")
        for index in range(self.getNumConstraints()):
            (particle1, particle2, distance) = self.getConstraintParameters(index)
            r += "%8d %8d %s" % (particle1, particle2, str(distance))
        r += "\n"

        # Show forces.
        r += "Forces:\n"
        for force in self.forces:
            r += str(force)

        return r

    def __add__(self, other):
        """
        Binary concatenation of two systems.

        The atoms of the second system appear, in order, after the atoms of the first system
        in the new combined system.

        USAGE

        combined_system = system1 + system2

        NOTES

        Both systems must have identical ordering of Force terms.
        Any non-particle settings from the first System override those of the second, if they differ.

        EXAMPLES

        Concatenate two systems in vacuum.

        >>> import testsystems
        >>> [system1, coordinates1] = testsystems.LennardJonesFluid()
        >>> system1 = System(system1) # convert to pure Python system
        >>> [system2, coordinates2] = testsystems.LennardJonesFluid()
        >>> system2 = System(system2) # convert to pure Python system
        >>> combined_system = system1 + system2

        """
        system = copy.deepcopy(self)
        system += other
        return system

    def __iadd__(self, other):
        """
        Append specified system.

        USAGE

        system += additional_system

        NOTES

        Both systems must have identical ordering of Force terms.
        Any non-particle settings from the first System override those of the second, if they differ.

        EXAMPLES

        Append atoms from a second system to the first.

        >>> import testsystems
        >>> [system1, coordinates1] = testsystems.LennardJonesFluid()
        >>> system1 = System(system1) # convert to pure Python system
        >>> [system2, coordinates2] = testsystems.LennardJonesFluid()
        >>> system2 = System(system2) # convert to pure Python system
        >>> system1 += system2

        """
        # Check to make sure both systems are compatible.
        if not isinstance(other, System):
            raise ValueError("both arguments must be System objects")
        if (self.nforces != other.nforces):
            raise ValueError("both System objects must have identical number of Force classes")
        for (force1,force2) in zip(self.forces, other.forces):
            if type(force1) != type(force2):
                raise ValueError("both System objects must have identical ordering of Force classes")

        # TODO: Make sure other system is Pythonic instance.

        # Combine systems.
        for mass in other.masses:
            mass_copy = copy.deepcopy(mass)
            self.masses.append(mass_copy)
        for constraint in other.constraints:
            constraint_copy = copy.deepcopy(constraint)
            self.constraints.append(constraint_copy)
        for (force1, force2) in zip(self.forces, other.forces):
            force2_copy = copy.deepcopy(force2)
            offset = self.nparticles
            force1._appendForce(force2_copy, offset)

        return

    #==========================================================================
    # Methods for selecting atom groups based on Metadata
    #==========================================================================

    def selectAtomGroup(self, atomname=None, atomtype=None, resname=None, resnum=None, \
                              atomnum=None, atomcharge=None, freeEnergyAtom=False):

        """Select particles according to specified metadata attributes.  Selections can be lists,
        or individual entries.

        Examples:

        Select all atomnames that are either N, CA or C, and in residue 4

        >>> selection = selectAtomGroup(atomname=['N','CA','C'], resnum=4)
        >>> selection = selectAtomGroup(atomname=['N','CA','C'], resnum=[4])

        RETURNS an AtomGroup object
        """

        # Convert the specification to lists (or None)
        if not isinstance(atomname,list) and (atomname!=None):
            AtomNameList = [ atomname ]
        if not isinstance(atomtype,list) and (atomtype!=None):
            AtomTypeList = [ atomtype ]
        if not isinstance(resname,list) and (resname!=None):
            ResidueNameList = [ resname ]
        if not isinstance(resnum,list) and (resnum!=None):
            ResidueNumList = [ resnum ]
        if not isinstance(atomcharge,list) and (atomcharge!=None):
            AtomChargeList = [ atomcharge ]
        if not isinstance(freeEnergyAtom,list) and (freeEnergyAtom!=None):
            freeEnergyList = [ freeEnergyAtom ]

        group = AtomGroup()
        for metadatum in self.metadata:
            Keep = True
            Keep *= ((metadatum.atomname in AtomNameList) or AtomNameList==None)
            Keep *= ((metadatum.atomtype in AtomTypeList) or AtomTypeList==None)
            Keep *= ((metadatum.resname in ResidueNameList) or ResidueNameList==None)
            Keep *= ((metadatum.resnum in ResidueNumList) or ResidueNumList==None)
            Keep *= ((metadatum.atomcharge in AtomChargeList) or AtomChargeList==None)
            Keep *= ((metadatum.freeEnergyAtom in freeEnergyList) or freeEnergyList==None)
            if Keep:
                group.addIndex(metadatum.particle)


    #==========================================================================
    # Methods for get and set of MetaData quantities
    #==========================================================================

    # by particle (index)

    def getAtomNameByParticle(self):
        pass

    def setAtomNameByParticle(self):
        pass

    def getAtomTypeByParticle(self):
        pass

    def setAtomTypeByParticle(self):
        pass

    def getResidueNameByParticle(self):
        pass

    def setResidueNameByParticle(self):
        pass

    def getResidueNumByParticle(self):
        pass

    def setResidueNumByParticle(self):
        pass

    def getAtomChargeByParticle(self):
        pass

    def setAtomChargeByParticle(self):
        pass

    def getCommentByParticle(self):
        pass

    def setCommentByParticle(self):
        pass

    def getFreeEnergyAtomByParticle(self):
        pass

    def setFreeEnergyAtomByParticle(self):
        pass


    # by atomgroup (index)

    def setAtomNameByAtomGroup(self):
        pass

    def getAtomNameByAtomGroup(self):
        pass

    def getAtomTypeByAtomGroup(self):
        pass

    def setAtomTypeByAtomGroup(self):
        pass

    def getResidueNameByAtomGroup(self):
        pass

    def setResidueNameByAtomGroup(self):
        pass

    def getResidueNumByAtomGroup(self):
        pass

    def setResidueNumByAtomGroup(self):
        pass

    def getAtomChargeByAtomGroup(self):
        pass

    def setAtomChargeByAtomGroup(self):
        pass

    def getCommentByAtomGroup(self):
        pass

    def setCommentByAtomGroup(self):
        pass

    def getFreeEnergyAtomByAtomGroup(self):
        pass

    def setFreeEnergyAtomByAtomGroup(self):
        pass


class TopAtom(object):
    """
    ;     nr type        resnr residue   atom   cgnr charge     mass       typeB chargeB massB
            1 amber99_39      1   NPRO      N      1 -0.202     14.01      ; qtot -0.202

    """
    @accepts_compatible_units(None, None, None, None, None, None, None, units.elementary_charge, units.amu, units.nanometers, units.kilojoules_per_mole, None, None)
    def __init__(self, ID, atomtype, resnum, resname, atomname, Z, cgnr, charge, mass, sigma, epsilon, exclusionList=list(), freeEnergyAtom=False):
        self.ID = ID
        self.atomname = atomname
        self.Z = Z
        self.metadata = Metadata(atomtype, resnum, resname, cgnr)
        self.mass = mass
        self.charge = charge
        self.sigma = sigma
        self.epsilon = epsilon
        self.exclusionList = exclusionList
        self.freeEnergyAtom = freeEnergyAtom

    def getID(self):
        return self.ID

    def setID(self, ID):
        self.ID = ID
        return


class Metadata(object):
    """
    A container class for storing atomtypes, residue names, residue numbers and charge numbers

    """

    def __init__(self, atomtype=None, resnum=None, resname=None, cgnr=None):

        self.atomtype = atomtype
        self.resnum = resnum
        self.resname = resname
        self.cgnr = cgnr
        return

    def getAtomType(self):
        return self.atomtype

    def setAtomType(self, atomtype):
        self.atomtype = atomtype
        return

    def getPtype(self):
        return self.ptype

    def setPtype(self, ptype):
        self.ptype = ptype
        return

    def getResidueName(self):
        return self.resname

    def setResidueName(self, resname):
        self.resname = resname
        return

    def getResidueNum(self):
        return self.resnum

    def setResidueNum(self, resnum):
        self.resnum = resnum
        return

    def getCGNR(self):
        return self.cgnr

    def setCGNR(self, cgnr):
        self.cgnr = cgnr
        return


class AtomGroup(object):
    """
    A class that stores a list of atom indices and an atomgroup name.
    Selection methods will return AtomGroup objects as the result of a selection query.

    Metadata properties can also be set using AtomGroup objects

    """

    def __init__(self, atomgroup=None):
        """Initialize an AtomGroup class.  If an input class is provided, instantiate using the information from it."""

        self.indices = list()
        self.name = None

        if atomgroup:
            self.indices = copy.deepcopy(atomgroup.indices)
            self.name = copy.deepcopy(atomgroup.name)

        return

    def addIndex(self, index):
        self.indices.append(index)
