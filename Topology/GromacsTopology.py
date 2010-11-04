#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import sys, os, tempfile, string, copy, pdb
from System import *
from Decorators import *
from GromacsParameter import *
from GromacsTopologyParser import *

ForceDictionary = {
"BondForce":"bonds",
"G96BondForce":"bonds",
"MorseBondForce":"bonds",
"CubicBondForce":"bonds",
"ConnectionBondForce":"bonds",
"HarmonicBondForce":"bonds",
"LJ1Force":"pairs",
"LJ2Force":"pairs",
"LJNBForce":"pairs_nb",
"AngleForce":"angles",
"G96AngleForce":"angles",
"CrossBondBondAngleForce":"angles",
"CrossBondAngleAngleForce":"angles",
"UreyBradleyAngleForce":"angles",
"QuarticAngleForce":"angles",
"PeriodicTorsionForce":"dihedrals",
"NonPeriodicTorsionForce":"dihedrals",
"RBTorsionForce":"dihedrals",
"FourierTorsionForce":"dihedrals",
"SettleForce":"settles"
}

# Classes

class GromacsTopology(object):

    """
    DESCRIPTION

    A derived class of Topology.py but with:

    1) extra attributes for gromacs topology parameters, #ifdefs, #defines
    2) extra methods for reading and writing to GROMACS *.top files


    REQUIREMENTS

    * Needs an installation of the simtk.units package.  See: https://simtk.org/home/python_units

    * The environment variable $GMXLIB needs to be defined, listing the directory where Gromacs forcefield and *.itp files
      are stored.  This can be added to .bash_profile, and sourced, for example:

      export GMXLIB="/home/my-account/local/gromacs-3.1.4/share/top"

    TO DO

    * Some parameters will be specified using the #included *.itp, while others are specified in the file.
    * Be able to write specified sets of parameters to *.itp files at will.

    NOTES

    """

    def __init__(self, topology=None, topfile=None, name=None):
        """
        Create a new GromacsTopology object (a derived Topology class)

        If a GromacsTopology object is specified, it will be queried to construct the class.
        If a *.top filename is specified, it will be read in.

        """

        # name of GromacsTopology (will also be used as name of [ system ] in the *.top file)
        self.name = name
        if self.name == None:
            self.name = "Untitled"

        self.defaults = list()

        self.parameters = list() # These are GromacsParameter[i] objects, akin to Force objects,
                                 # but holding Force parameters for various atomtype combinations 
        self.molecules = list()  # molecules[i] is the ith Topology object
        self.index = list()

        if (topfile):
            self.readTopologyFile(topfile=topfile)
        else:
            self.GromacsTopologyFileObject = None

        if (topology):
            if type(topology).__name__ == 'list':
                for top in topology:
                    self.molecules.append(top)
            else:
                self.molecules.append(topology)
        return

    def setName(self, name):
        self.name = name
        return

    def addMolecule(self, molecule):
        self.molecules.append(molecule)
        return

    def readTopologyFile(self, topfile):
        """Read in the contents of the topology file.
        NOTES

        """

        # Read in the *.top file, creating a new GromacsTopologyFile object
        self.GromacsTopologyFileObject = GromacsTopologyParser(topfile=topfile)

	parmIndex = list()
        # Translate all the parameter directives to GromacsParameter objects
        for parm in self.GromacsTopologyFileObject.parameters:
            self.parameters.append(self.ConvertParameterDirectiveToGromacsParameter(parm))
            parmIndex.append(parm.getName())
        self.index.append(parmIndex)
	
	# Create Molecules Parameter
	#for line in self.SystemDirectives[1].lines
	#	fields, comment  = line.splitLine(line)
	#	compounds.append(fields[0])
	#	molList.append(fields[1])
	#newGroParm = GromacsMoleculesParameterInfo(compound, molList, comment)
		

        #incredibly sloppy function that combines the contents of duplicate parameter entries
        set = {}
        copy = [set.setdefault(e,e) for e in self.index[0] if e not in set]

        reallist = list()
        for i in range(len(self.index[0])):
            walid = list()
            for j in range(i, len(self.index[0])):
                if self.index[0][i] == self.index[0][j]:
                    parmInfo = self.parameters[j]
                    walid.extend(parmInfo)
            reallist.append(walid)

        self.index[0] = copy
        while len(reallist) > len(self.index[0]):
            reallist.pop()
        self.parameters = reallist

        # Translate each set of molecule directives to Topology objects
        for mol in range(len(self.GromacsTopologyFileObject.molecules)):
            self.molecules.append( self.ConvertMoleculeDirectivesToTopology( self.GromacsTopologyFileObject.molecules[mol] ) )

        return

    def grabMissingInfo(self, molecule, header, func, particle1, particle2, particle3 = None, particle4 = None):
        """
        Locates the indicies of the missing data from parameters to be plugged into molecule definitions.

        """

        if "bonds" in header:

            # Find the location of bondtypes in the index matrix
            i = self.index[0].index('[ bondtypes ]')

            # Filter out the bondtypes with different func numbers
            bondtypes = list()
            for bondtype in self.parameters[i]:
                if bondtype.func == func:
                    bondtypes.append(bondtype)

            if len(bondtypes) == 0:
                raise "Did not find any bondtypes"

            # Convert each particle to an atomname in [ atoms ]
            for atom in molecule.atoms:
                if particle1 == int(atom.ID):
                    atomname1 = atom.metadata.atomtype
                if particle2 == int(atom.ID):
                    atomname2 = atom.metadata.atomtype

            # Convert each atomname to a bondtype
            for atomtype in self.parameters[1]:
                if atomname1 == atomtype.name:
                    bondtype1 = atomtype.bondtype
                if atomname2 == atomtype.name:
                    bondtype2 = atomtype.bondtype

            # Convert two bondtypes to a bond parameter
            for bondtype in bondtypes:
                if (bondtype1 == bondtype.atomtype1 and bondtype2 == bondtype.atomtype2) or (bondtype1 == bondtype.atomtype2 and bondtype2 == bondtype.atomtype1):
                    j = self.parameters[i].index(bondtype)

        if "pairs" in header:

            # Find the location of pairtypes in the index matrix
            i = self.index[0].index('[ pairtypes ]')

            # Filter out the pairtypes with different func numbers
            pairtypes = list()
            for pairtype in self.parameters[i]:
                if pairtype.func == func:
                    pairtypes.append(pairtype)

            if len(pairtypes) == 0:
                raise "Did not find any pairtypes"

            # Convert each particle to an atomname in [ atoms ]
            for atom in molecule.atoms:
                if particle1 == int(atom.ID):
                    atomname1 = atom.metadata.atomtype
                if particle2 == int(atom.ID):
                    atomname2 = atom.metadata.atomtype

            # Convert each atomname to a bondtype
            for atomtype in self.parameters[1]:
                if atomname1 == atomtype.name:
                    bondtype1 = atomtype.bondtype
                if atomname2 == atomtype.name:
                    bondtype2 = atomtype.bondtype

            # Convert two bondtypes to a pair parameter
            for angletype in angletypes:
                if (bondtype1 == angletype.atomtype1 and bondtype2 == angletype.atomtype2 and bondtype3 == angletype.atomtype3) or (bondtype1 == angletype.atomtype3 and bondtype2 == angletype.atomtype2 and bondtype3 == angletype.atomtype1):
                    j = self.parameters[i].index(pairtype)

        if "angles" in header:

            # Find the location of angletypes in the index matrix
            i = self.index[0].index('[ angletypes ]')

            # Filter out the angletypes with different func numbers
            angletypes = list()
            for angletype in self.parameters[i]:
                if angletype.func == func:
                    angletypes.append(angletype)

            if len(angletypes) == 0:
                raise "Did not find any angletypes"

            # Convert each particle to an atomname in [ atoms ]
            for atom in molecule.atoms:
                if particle1 == int(atom.ID):
                    atomname1 = atom.metadata.atomtype
                if particle2 == int(atom.ID):
                    atomname2 = atom.metadata.atomtype
                if particle3 == int(atom.ID):
                    atomname3 = atom.metadata.atomtype

            # Convert each atomname to a bondtype
            for atomtype in self.parameters[1]:
                if atomname1 == atomtype.name:
                    bondtype1 = atomtype.bondtype
                if atomname2 == atomtype.name:
                    bondtype2 = atomtype.bondtype
                if atomname3 == atomtype.name:
                    bondtype3 = atomtype.bondtype

            # Convert three bondtypes to an angle parameter
            for angletype in angletypes:
                if (bondtype1 == angletype.atomtype1 and bondtype2 == angletype.atomtype2 and bondtype3 == angletype.atomtype3) or (bondtype1 == angletype.atomtype3 and bondtype2 == angletype.atomtype2 and bondtype3 == angletype.atomtype1):
                    j = self.parameters[i].index(angletype)

        if "dihedrals" in header:

            # Find the location of dihedraltypes in the index matrix
            i = self.index[0].index('[ dihedraltypes ]')

            # Filter out the dihedraltypes with different func numbers
            dihedraltypes = list()
            for dihedraltype in self.parameters[i]:
                if dihedraltype.func == func:
                    dihedraltypes.append(dihedraltype)

            if len(dihedraltypes) == 0:
                raise "Did not find any dihedraltypes"

            # Convert each particle to an atomname in [ atoms ]
            for atom in molecule.atoms:
                if particle1 == int(atom.ID):
                    atomname1 = atom.metadata.atomtype
                if particle2 == int(atom.ID):
                    atomname2 = atom.metadata.atomtype
                if particle3 == int(atom.ID):
                    atomname3 = atom.metadata.atomtype
                if particle4 == int(atom.ID):
                    atomname4 = atom.metadata.atomtype

            # Convert each atomname to a bondtype
            for atomtype in self.parameters[1]:
                if atomname1 == atomtype.name:
                    bondtype1 = atomtype.bondtype
                if atomname2 == atomtype.name:
                    bondtype2 = atomtype.bondtype
                if atomname3 == atomtype.name:
                    bondtype3 = atomtype.bondtype
                if atomname4 == atomtype.name:
                    bondtype4 = atomtype.bondtype

            # Convert four bondtypes to a dihedral parameter
            for dihedraltype in dihedraltypes:
                if ((bondtype1 == dihedraltype.atomtype1) or ("X" == dihedraltype.atomtype1)) and ((bondtype2 == dihedraltype.atomtype2)  or ("X" == dihedraltype.atomtype2)) and ((bondtype3 == dihedraltype.atomtype3) or ("X" == dihedraltype.atomtype3)) and ((bondtype4 == dihedraltype.atomtype4) or ("X" == dihedraltype.atomtype4)):
                    j = self.parameters[i].index(dihedraltype)
                elif ((bondtype1 == dihedraltype.atomtype4) or ("X" == dihedraltype.atomtype4)) and ((bondtype2 == dihedraltype.atomtype3) or ("X" == dihedraltype.atomtype3)) and ((bondtype3 == dihedraltype.atomtype2) or ("X" == dihedraltype.atomtype2)) and ((bondtype4 == dihedraltype.atomtype1) or ("X" == dihedraltype.atomtype1)):
                    j = self.parameters[i].index(dihedraltype)

        if "constraints" in header:

            # Find the location of constrainttypes in the index matrix
            i = self.index[0].index('[ constrainttypes ]')

            # Filter out the constrainttypes with different func numbers
            constrainttypes = list()
            for constrainttype in self.parameters[i]:
                if constrainttype.func == func:
                    constrainttypes.append(constrainttype)

            if len(constrainttypes) == 0:
                raise "Did not find any constrainttypes"

            # Convert each particle to an atomname in [ atoms ]
            for atom in molecule.atoms:
                if particle1 == int(atom.ID):
                    atomname1 = atom.metadata.atomtype
		if particle2 == int(atom.ID):
                    atomname2 = atom.metadata.atomtype

            # Convert each atomname to a bondtype
            for atomtype in self.parameters[1]:
                if atomname1 == atomtype.name:
                    bondtype1 = atomtype.bondtype
                if atomname2 == atomtype.name:
                    bondtype2 = atomtype.bondtype

            # Convert two bondtypes to a constraint parameter
            for constrainttype in constrainttypes:
                if (bondtype1 == constrainttype.atomtype1 and bondtype2 == constrainttype.atomtype2 and bondtype3 == constrainttype.atomtype3) or (bondtype1 == constrainttype.atomtype3 and bondtype2 == constrainttype.atomtype2 and bondtype3 == constrainttype.atomtype1):
                    j = self.parameters[i].index(constrainttype)

        try:
            j
        except:
            raise "%(header)s Type Not Found"%vars()

        return i, j

    def splitline(self, line):
        """
        Function to split a line into a comment and a list of parameters

        """

        if ";" in line:
            fields = line.split(';')
            line = fields[0]
            comment = fields[1].strip()
        else:
            comment = ''
        fields = line.split()

        return fields, comment

    def findMoleculeDirective(self, currentMoleculeDirective, molecule):
        """
        Function to determine if a molecule directive has already been created
        if it has been created then it will just set the currentMoleculeDirective to the created one
        if not then it will use the new currentMoleculeDirective

        """

        if len(molecule.forces) == 0:
            molecule.forces.append(currentMoleculeDirective)
        else:
            for moleculeDirective in molecule.forces:
                if type(currentMoleculeDirective) == type(moleculeDirective):
                    index = molecule.forces.index(moleculeDirective)
                    currentMoleculeDirective = molecule.forces[index]
            if currentMoleculeDirective.getNumForces() == 0:
                molecule.forces.append(currentMoleculeDirective)
        return currentMoleculeDirective

    def ConvertMoleculeDirectivesToTopology(self, MoleculeDirectives):
        """Takes a list of MoleculeDirectives (Directive objects), and uses it
        to create a TopologySystem object with the correct forces.

        """

        # For each of the molecule directives, tell you how atoms are bonded.
        # Going through each of the directives
        # Start in a directive, go down the list of atoms, there are atom indices
        # each index corresponds to an atom of a particular type
        # each will have atom type.
        # How each of the atoms are bonded in a force
        # option in writing a gro file -- if you read in any include file atoms,
        # how do you know whether it was an include
        # create an extra variable in the parameter info object / print me

        # Each Parameter info object has a flag that indicates if it
        # comes from an include file or not.  The topology system
        # object also has such a flag.  When you create the force
        # objects from the parameter info objects, you check the
        # topology system flag.  If any parameter objects that create
        # force object are not from include files, then the entire
        # topology object must be printed.

        molecule = Topology()
        molIndex = list()


        # Fill in system default parameters
        try:
            molecule.NBFunc = self.parameters[0][0].func
            molecule.CombinationRule = self.parameters[0][0].cr
            #genpairs = self.parameters[0][0].genpairs
            molecule.LJCorrection = self.parameters[0][0].fudgeLJ
            molecule.CoulombCorrection = self.parameters[0][0].fudgeQQ
        except:
            pass
        # molecule has name, atoms (list), forces (list), constraints (list), atomgroups (list)

	for directive in MoleculeDirectives:
	    if len(MoleculeDirectives)>0:
                # Fill in moleculetype
                if "moleculetype" in directive.name:
                    for line in directive.lines:
                        fields = line.split()
                        molname = fields[0]
                        nrexcl = fields[1]
                        molecule.name = molname
                        molecule.nrexcl = nrexcl
                
                # Fill in atoms
                if "atoms" in directive.name:
                    for line in directive.lines:
                        fields, comment = self.splitline(line)
                        particle = int(fields[0])
                        atomtype = fields[1]
                        resnum = int(fields[2])
                        resname = fields[3]
                        atomname = fields[4]
                        cgnr = int(fields[5])
                        charge = float(fields[6]) * units.elementary_charge
                        try:
                            mass = float(fields[7]) * units.amu
                        except:
                            mass = 0.0000 * units.amu ####################### WHY DO WE HAVE THIS?
                        atomnum = particle
                        for atomtypenis in self.parameters[1]:
                            Z = atomtypenis.Z
                            if atomtype == atomtypenis.name:
                                if self.parameters[0][0].cr == 1:
                                    sigma = (atomtypenis.W / atomtypenis.V)**(1/6)
                                    epsilon = atomtypenis.V / (4*sigma**6)
                                elif self.parameters[0][0].cr == 2 or self.parameters[0][0].cr == 3:
                                    sigma = atomtypenis.V
                                    epsilon = atomtypenis.W
                        # B state parameters
                        try:
                                Bfields = comment.split()
                                Bname = Bfields[0]
                                Bcharge = float(Bfields[1]) * units.elementary_charge
                                if len(Bfields) == 3:
                                        Bmass = float(Bfields[2]) * units.amu
                                        Batom = TopAtom(particle, atomtype, resnum, resname, atomname, Z, cgnr, Bcharge, Bmass, sigma, epsilon)
                                else:
                                        Batom = TopAtom(particle, atomtype, resnum, resname, atomname, Z, cgnr, Bcharge, mass, sigma, epsilon)
    
                                molecule.BAtoms.append(Batom)
                        except:
                                pass	
                        atom = TopAtom(particle, atomtype, resnum, resname, atomname, Z, cgnr, charge, mass, sigma, epsilon)
                        molecule.atoms.append(atom)
    
                # Fill in bonds
                if "bonds" in directive.name:
                    for line in directive.lines:
                        fields, comment = self.splitline(line)
                        particle1 = int(fields[0])
                        particle2 = int(fields[1])
                        func = int(fields[2])
    
                        if func == 1:
                            newMoleculeDirective = BondForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 3:
                                length = float(fields[3]) * units.nanometers
                                k = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**(-2)
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                length = self.parameters[i][j].distance
                                k = self.parameters[i][j].kspring
                            parameters = [particle1, particle2, length, k]
    
                        if func == 2:
                            newMoleculeDirective = G96BondForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 3:
                                length = float(fields[3]) * units.nanometers
                                k = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**(-4)
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                length = self.parameters[i][j].distance
                                k = self.parameters[i][j].kspring
                            parameters = [particle1, particle2, length, k]
    
                        if func == 3:
                            newMoleculeDirective = MorseBondForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 3:
                                length = float(fields[3]) * units.nanometers
                                D = float(fields[4]) * units.kilojoules_per_mole
                                beta = float(fields[5]) * units.nanometers**(-1)
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                D = self.parameters[i][j].D
                                beta = self.parameters[i][j].beta
                            parameters = [particle1, particle2, length, D, beta]
    
                        if func == 4:
                            newMoleculeDirective = CubicBondForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 3:
                                length = float(fields[3]) * units.nanometers
                                C2 = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**(-2)
                                C3 = float(fields[5]) * units.kilojoules_per_mole * units.nanometers**(-3)
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                C2 = self.parameters[i][j].C2
                                C3 = self.parameters[i][j].C3
                            parameters = [particle1, particle2, C2, C3]
    
                        if func == 5:
                            newMoleculeDirective = ConnectionBondForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            parameters = [particle1, particle2]
    
                        if func == 6:
                            newMoleculeDirective = HarmonicBondForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 3:
                                length = float(fields[3]) * units.nanometers
                                k = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**(-2)
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                length = self.parameters[i][j].distance
                                k = self.parameters[i][j].kspring
                            parameters = [particle1, particle2, length, k]
    
                        if func == 7:
                            raise "Bond type not implemented"
    
                        if func == 8:
                            raise "Bond type not implemented"
    
                        if func == 9:
                            raise "Bond type not implemented"
    
                        currentMoleculeDirective.addForce(*parameters)
    
                # Fill in pairs
                if "pairs" in directive.name:
                    for line in directive.lines:
                        fields, comment = self.splitline(line)
                        atomtype1 = int(fields[0])
                        atomtype2 = int(fields[1])
                        func = int(fields[2])
    
                        if func == 1:
                            newMoleculeDirective = LJ1Force()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 3:
                                if self.parameters[0][0].cr == 1:
                                    V = float(fields[3]) * units.kilojoules_per_mole * units.nanometers**6
                                    W = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**12
                                elif self.parameters[0][0].cr == 2 or self.parameters[0][0].cr == 3:
                                    V = float(fields[3]) * units.nanometers
                                    W = float(fields[4]) * units.kilojoules_per_mole
                            elif "pairtypes" in self.index[0]:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                V = self.parameters[i][j].V
                                W = self.parameters[i][j].W
                            else:
                                V = None * units.dimensionless
                                W = None * units.dimensionless
                            parameters = [atomtype1, atomtype2, V, W]
    
                        if func == 2:
                            newMoleculeDirective = LJ2Force()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 3:
                                fudgeQQ = float(fields[3])
                                qi = float(fields[4]) * units.elementary_charge
                                qj = float(fields[5]) * units.elementary_charge
                                if self.parameters[0][0].cr == 1:
                                    V = float(fields[6]) * units.nanometers
                                    W = float(fields[7]) * units.kilojoules_per_mole
                                elif self.parameters[0][0].cr == 2 or self.parameters[0][0].cr == 3:
                                    V = float(fields[6]) * units.kilojoules_per_mole * units.nanometers**6
                                    W = float(fields[7]) * units.kilojoules_per_mole * units.nanometers**12
                            elif "pairtypes" in self.index[0]:
                                i, j = grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                V = self.parameters[i][j].V
                                W = self.parameters[i][j].W
                            else:
                                V = None * units.dimensionless
                                W = None * units.dimensionless
                            parameters = [atomtype1, atomtype2, qi, qj, V, W]
    
                        currentMoleculeDirective.addForce(*parameters)
    
                # Fill in nonbonded pairs
                if "pairs_nb" in directive.name:
                    for line in directive.lines:
                        fields, comment = self.splitline(line)
                        func = int(fields[2])
    
                        if func == 1:
                            newMoleculeDirective = LGNBForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 3:
                                qi = float(fields[3]) * units.elementary_charge
                                qj = float(fields[4]) * units.elementary_charge
                                if self.parameters[0][0].cr == 1:
                                    V = float(fields[5]) * units.nanometers**6
                                    W = float(fields[6]) * units.kilojoules_per_mole
                                elif self.parameters[0][0].cr == 2 or self.parameters[0][0].cr == 3:
                                    V = float(fields[5]) * units.kilojoules_per_mole * units.nanometers**6
                                    W = float(fields[6]) * units.kilojoules_per_mole * units.nanometers**12
                            elif "pairtypes" in self.index[0]:
                                i, j = grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                V = self.parameters[i][j].V
                                W = self.parameters[i][j].W
                            else:
                                V = None * units.dimensionless
                                W = None * units.dimensionless
                            parameters = [atomtype1, atomtype2, qi, qj, V, W]
    
                        currentMoleculeDirective.addForce(*parameters)
    
    
                # Fill in angles
                if "angles" in directive.name:
                    for line in directive.lines:
                        fields, comment = self.splitline(line)
                        particle1 = int(fields[0])
                        particle2 = int(fields[1])
                        particle3 = int(fields[2])
                        func = int(fields[3])
    
                        if func == 1:
                            newMoleculeDirective = AngleForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 4:
                                angle = float(fields[4]) * units.degrees
                                k = float(fields[5]) * units.kilojoules_per_mole * units.radians**(-2)
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3)
                                angle = self.parameters[i][j].theta
                                k = self.parameters[i][j].k
                            parameters = [particle1, particle2, particle3, angle, k]
    
                        if func == 2:
                            newMoleculeDirective = G96AngleForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 4:
                                angle = float(fields[4]) * units.degrees
                                k = float(fields[5]) * units.kilojoules_per_mole
                            else:
                                i, j = grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3)
                                angle = self.parameters[i][j].theta
                                k = self.parameters[i][j].k
                            parameters = [particle1, particle2, particle3, angle, k]
    
                        if func == 3:
                            newMoleculeDirective = CrossBondBondAngleForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 4:
                                r1e = float(fields[4]) * units.nanometers
                                r2e = float(fields[5]) * units.nanometers
                                krr = float(fields[6]) * units.kilojoules_per_mole * units.nanometers**(-2)
                            else:
                                i, j = grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3)
                                r1e = self.parameters[i][j].r1e
                                r2e = self.parameters[i][j].r2e
                                krr = self.parameters[i][j].krr
                            parameters = [particle1, particle2, particle3, r1e, r2e, krr]
    
                        if func == 4:
                            newMoleculeDirective = CrossBondAngleAngleForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 4:
                                r1e = float(fields[4]) * units.nanometers
                                r2e = float(fields[5]) * units.nanometers
                                r3e = float(fields[6]) * units.nanometers
                                kr = float(fields[7]) * units.kilojoules_per_mole * units.nanometers
                            else:
                                i, j = grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3)
                                r1e = self.parameters[i][j].r1e
                                r2e = self.parameters[i][j].r2e
                                r3e = self.parameters[i][j].r3e
                                kr = self.parameters[i][j].kr
                            parameters = [particle1, particle2, particle3, r1e, r2e, r3e, kr]
    
                        if func == 5:
                            newMoleculeDirective = UreyBradleyAngleForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 4:
                                angle = float(fields[4]) * units.degrees
                                k = float(fields[5]) * units.kilojoules_per_mole
                                r13 = float(fields[6]) * units.nanometers
                                kUB = float(fields[7]) * units.kilojoules_per_mole
                            else:
                                i, j = grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3)
                                angle = self.parameters[i][j].theta
                                k = self.parameters[i][j].k
                                r13 = self.parameters[i][j].r13
                                kUB = self.parameters[i][j].kUB
                            parameters = [particle1, particle2, particle3, angle, k, r13, kUB]
    
                        if func == 6:
                            newMoleculeDirective = QuarticAngleForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 4:
                                angle = fields[4]
                                C0 = float(fields[5]) * units.kilojoules_per_mole
                                C1 = float(fields[6]) * units.kilojoules_per_mole * units.radians**(-1)
                                C2 = float(fields[7]) * units.kilojoules_per_mole * units.radians**(-2)
                                C3 = float(fields[8]) * units.kilojoules_per_mole * units.radians**(-3)
                                C4 = float(fields[9]) * units.kilojoules_per_mole * units.radians**(-4)
                            else:
                                i, j = grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3)
                                C0 = self.parameters[i][j].C0
                                C1 = self.parameters[i][j].C1
                                C2 = self.parameters[i][j].C2
                                C3 = self.parameters[i][j].C3
                                C4 = self.parameters[i][j].C4
                            parameters = [particle1, particle2, particle3, angle, C0, C1, C2, C3, C4]
    
                        if func == 8:
                            raise "Angle type not implemented"
    
                        currentMoleculeDirective.addForce(*parameters)
    
                # Fill in dihedrals
                if "dihedrals" in directive.name:
                    for line in directive.lines:
                        fields, comment = self.splitline(line)
                        particle1 = int(fields[0])
                        particle2 = int(fields[1])
                        particle3 = int(fields[2])
                        particle4 = int(fields[3])
                        func = int(fields[4])
    
                        if func == 1:
                            newMoleculeDirective = PeriodicTorsionForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 5:
                                phi = float(fields[5]) * units.degrees
                                kphi = float(fields[6]) * units.kilojoules_per_mole
                                multiplicity = int(fields[7])
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3, particle4)
                                phi = self.parameters[i][j].phi
                                kphi = self.parameters[i][j].kphi
                                multiplicity = self.parameters[i][j].multiplicity
                            parameters = [particle1, particle2, particle3, particle4, multiplicity, phi, kphi]
    
                        if func == 2:
                            newMoleculeDirective = NonPeriodicTorsionForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 5:
                                xi = float(fields[5]) * units.degrees
                                kxi = float(fields[6]) * units.kilojoules_per_mole * units.radians**(-2)
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3, particle4)
                                xi = self.parameters[i][j].xi
                                kxi = self.parameters[i][j].kxi
                            parameters = [particle1, particle2, particle3, particle4, xi, kxi]
    
                        if func == 3:
                            newMoleculeDirective = RBTorsionForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 5:
                                C0 = float(fields[5]) * units.kilojoules_per_mole
                                C1 = float(fields[6]) * units.kilojoules_per_mole
                                C2 = float(fields[7]) * units.kilojoules_per_mole
                                C3 = float(fields[8]) * units.kilojoules_per_mole
                                C4 = float(fields[9]) * units.kilojoules_per_mole
                                C5 = float(fields[10]) * units.kilojoules_per_mole
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3, particle4)
                                C0 = self.parameters[i][j].C0
                                C1 = self.parameters[i][j].C1
                                C2 = self.parameters[i][j].C2
                                C3 = self.parameters[i][j].C3
                                C4 = self.parameters[i][j].C4
                                C5 = self.parameters[i][j].C5
                            parameters = [particle1, particle2, particle3, particle4, C0, C1, C2, C3, C4, C5]
    
                        if func == 5:
                            newMoleculeDirective = FourierTorsionForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                            if len(fields) > 5:
                                C0 = float(fields[5]) * units.kilojoules_per_mole
                                C1 = float(fields[6]) * units.kilojoules_per_mole
                                C2 = float(fields[7]) * units.kilojoules_per_mole
                                C3 = float(fields[8]) * units.kilojoules_per_mole
                                C4 = float(fields[9]) * units.kilojoules_per_mole
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3, particle4)
                                C0 = self.parameters[i][j].C0
                                C1 = self.parameters[i][j].C1
                                C2 = self.parameters[i][j].C2
                                C3 = self.parameters[i][j].C3
                                C4 = self.parameters[i][j].C4
                            parameters = [particle1, particle2, particle3, particle4, C0, C1, C2, C3, C4]
    
                        if func == 8:
                            raise "Dihedral type not implemented"
    
                        currentMoleculeDirective.addForce(*parameters)
    
                # Fill in exclusions
                if "exclusions" in directive.name:
                    for line in directive.lines:
                        fields, comment = self.splitline(line)
                        particle = int(fields.pop(0))
                        exclusionList = list()
                        for excludedParticle in fields:
                            exclusionList.append(int(excludedParticle))
                        for topAtom in molecule.atoms:
                            if topAtom.ID == particle:
                                topAtom.exclusionList = exclusionList
    
                # Fill in constraints
                if "constraints" in directive.name:
                    for line in directive.lines:
                        fields, comment = self.splitline(line)
                        particle1 = int(fields[0])
                        particle2 = int(fields[1])
                        func = int(fields[3])
    
                        if func == 1:
                            length = float(fields[4])
    
                        if func == 2:
                            lengthA = float(fields[4])
                            lengthB = float(fields[5])
    
                # Fill in settles
                if "settles" in directive.name:
                    for line in directive.lines:
                        fields, comment = self.splitline(line)
                        particle = int(fields[0])
                        func = int(fields[1])
                        doh = float(fields[2]) * units.nanometers
                        dhh = float(fields[3]) * units.nanometers
                        newMoleculeDirective = SettleForce()
                        currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule)
                        parameters = [particle, doh, dhh]
                        currentMoleculeDirective.addForce(*parameters)

            molIndex.append(directive.getName())

        self.index.append(molIndex)

        return molecule


    def ConvertParameterDirectiveToGromacsParameter(self, parmDirective):
        """Takes a parmDirective (Directive object), and uses it
        to create a GromacsParameter object to store the forcefield parameters.

        """

        parmList = list()
        index = 0
        while index < len(parmDirective.lines):
            line = parmDirective.lines[index]
            fields, comment = self.splitline(line)

            if "defaults" in parmDirective.getName():

                func = int(fields[0])
                cr = int(fields[1])
                genpairs = fields[2]
                fudgeLJ = float(fields[3])
                fudgeQQ = float(fields[4])
                newGroParm = GromacsDefaultParameterInfo(func, cr, genpairs, fudgeLJ, fudgeQQ, comment)

            if "atomtypes" in parmDirective.getName():
                i = 0
                name = fields[i]
                if len(fields)==8:
                    i+=1
                    bondtype = fields[i]
                    i+=1
                    Z = float(fields[i])
                elif len(fields)==7:
                    i+=1
                    bondtype = fields[i]
                    Z = ''
                else:
                    bondtype = ''
                    Z = ''
                mass = float(fields[i+1]) * units.amu
                charge = float(fields[i+2]) * units.elementary_charge
                ptype = fields[i+3]
                if self.parameters[0][0].cr == 1:
                    V = float(fields[i+4]) * units.kilojoules_per_mole*units.nanometers**6
                    W = float(fields[i+5]) * units.kilojoules_per_mole*units.nanometers**12
                    newGroParm = GromacsAtom1ParameterInfo(name, mass, charge, ptype, V, W, bondtype, Z,comment)
                elif self.parameters[0][0].cr == 2 or self.parameters[0][0].cr == 3:
                    V = float(fields[i+4]) * units.nanometers
                    W = float(fields[i+5]) * units.kilojoules_per_mole
                    newGroParm = GromacsAtom23ParameterInfo(name, mass, charge, ptype, V, W, bondtype, Z, comment)


            if "bondtypes" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                func = int(fields[2])

                if func == 1:
                    distance = float(fields[3]) * units.nanometers
                    kspring = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**(-2)
                    newGroParm = GromacsBondParameterInfo(atomtype1, atomtype2, func, distance, kspring, comment)

                if func == 2:
                    distance = float(fields[3]) * units.nanometers
                    kspring = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**(-4)
                    newGroParm = GromacsG96BondParameterInfo(atomtype1, atomtype2, func, distance, kspring, comment)

                if func == 3:
                    distance = float(fields[3]) * units.nanometers
                    D = float(fields[4]) * units.kilojoules_per_mole
                    beta = float(fields[5]) * units.nanometers**(-1)
                    newGroParm = GromacsMorseBondParameterInfo(atomtype1, atomtype2, func, distance, D, beta, comment)

                if func == 4:
                    distance = float(fields[3]) * units.nanometers
                    C2 = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**(-2)
                    C3 = float(fields[5]) * units.kilojoules_per_mole * units.nanometers**(-3)
                    newGroParm = GromacsCubicBondParameterInfo(atomtype1, atomtype2, func, distance, C2, C3, comment)

                if func == 5:
                    newGroParm = GromacsConnectionBondParameterInfo(atomtype1, atomtype2, func, comment)

                if func == 6:
                    distance = float(fields[3]) * units.nanometers
                    kspring = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**(-2)
                    newGroParm = GromacsHarmonicBondParameterInfo(atomtype1, atomtype2, func, distance, kspring, comment)

            if "pairtypes" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                func = int(fields[2])

                if func == 1:
                    if self.parameters[0][0] == 1:
                        V = float(fields[3]) * units.kilojoules_per_mole * units.nanometers**6
                        W = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**12
                        newGroParm = GromacsPairTypes11ParameterInfo(atomtype1, atomtype2, func, V, W, comment)
                    elif self.parameters[0][0] == 2 or self.parameters[0][0].cr == 3:
                        V = float(fields[3]) * units.nanometers**6
                        W = float(fields[4]) * units.kilojoules_per_mole
                        newGroParm = GromacsPairTypes123ParameterInfo(atomtype1, atomtype2, func, V, W, comment)

                if func == 2:
                    fudgeQQ = float(fields[3])
                    qi = float(fields[4]) * units.elementary_charge
                    qj = float(fields[5]) * units.elementary_charge
                    if self.parameters[0][0] == 1:
                        V = float(fields[6]) * units.nanometers**6
                        W = float(fields[7]) * units.kilojoules_per_mole
                        newGroParm = GromacsPairTypes2ParameterInfo(atomtype1, atomtype2, func, fudgeQQ, qi, qj, V, W, comment)
                    elif self.parameters[0][0] == 2 or self.parameters[0][0].cr == 3:
                        V = float(fields[6]) * units*kilojoules_per_mole * units.nanometers**6
                        W = float(fields[7]) * units.kilojoules_per_mole * units.nanometers**12
                        newGroParm = GromacsPairTypes2ParameterInfo(atomtype1, atomtype2, func, fudgeQQ, qi, qj, V, W, comment)

            # Possibly not needed?  It does not appear in grompp.h
            if "pairs_nb" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                func = int(fields[2])
                qi = float(fields[3]) * units.elementary_charge
                qj = float(fields[4]) * units.elementary_charge

                if func == 1:
                    if self.parameters[0][0] == 1 or self.parameters[0][0].cr == 3:
                        V = float(fields[5]) * units.kilojoules_per_mole * units.nanometers**6
                        W = float(fields[6]) * units.kilojoules_per_mole * units.nanometers**12
                        newGroParm = GromacsPairTypesNB13ParameterInfo(atomtype1, atomtype2, func, qi, qj, V, W, comment)
                    elif self.parameters [0] == 2:
                        V = float(fields[5]) * units.nanometers
                        W = float(fields[6]) * units.kilojoules_per_mole
                        newGroParm = GromacsPairTypesNB2ParameterInfo(atomtype1, atomtype2, func, qi, qj, V, W, comment)

            if "angletypes" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                atomtype3 = fields[2]
                func = int(fields[3])

                if func == 1:
                    theta = float(fields[4]) * units.degrees
                    k = float(fields[5]) * units.kilojoules_per_mole * units.radians**(-2)
                    newGroParm = GromacsAngleParameterInfo(atomtype1, atomtype2, atomtype3, func, theta, k, comment)

                if func == 2:
                    theta = float(fields[4]) * units.degrees
                    k = float(fields[5]) * units.kilojoules_per_mole
                    newGroParm = GromacsG96AngleParameterInfo(atomtype1, atomtype2, atomtype3, func, theta, k, comment)

                if func == 3:
                    r1e = float(fields[4]) * units.nanometers
                    r2e = float(fields[5]) * units.nanometers
                    krr = float(fields[6]) * units.kilojoules_per_mole*units.nanometers**(-2)
                    newGroParm = GromacsCrossBondBondAngleParameterInfo(atomtype1, atomtype2, atomtype3, func, r1e, r2e, krr, comment)

                if func == 4:
                    r1e = float(fields[4]) * units.nanometers
                    r2e = float(fields[5]) * units.nanometers
                    r3e = float(fields[6]) * units.nanometers
                    kr = float(fields[7]) * units.kilojoules_per_mole * units.nanometers**(-2)
                    newGroParm = GromacsCrossBondAngleAngleParameterInfo(atomtype1, atomtype2, atomtype3, func, r1e, r2e, r3e, kr, comment)

                if func == 5:
                    theta = float(fields[4]) * units.degrees
                    k = float(fields[5]) * units.kilojoules_per_mole
                    r13 = float(fields[6]) * units.nanometers
                    kUB = float(fields[7]) * units.kilojoules_per_mole
                    newGroParm = GromacsUreyBradleyAngleParameterInfo(atomtype1, atomtype2, atomtype3, func, theta, k, r13, kUB, comment)

                if func == 6:
                    theta = float(fields[4]) * units.degrees
                    C0 = float(fields[5]) * units.kilojoules_per_mole
                    C1 = float(fields[6]) * units.kilojoules_per_mole * units.radians**(-1)
                    C2 = float(fields[7]) * units.kilojoules_per_mole * units.radians**(-2)
                    C3 = float(fields[8]) * units.kilojoules_per_mole * units.radians**(-3)
                    C4 = float(fields[9]) * units.kilojoules_per_mole * units.radians**(-4)
                    newGroParm = GromacsQuarticAngleParameterInfo(atomtype1, atomtype2, atomtype3, func, theta, C0, C1, C2, C3, C4, comment)

            if "dihedraltypes" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                atomtype3 = fields[2]
                atomtype4 = fields[3]
                func = int(fields[4])

                if func == 1:
                    phi = float(fields[5]) * units.degrees
                    kphi = float(fields[6]) * units.kilojoules_per_mole
                    multiplicity = int(fields[7])
                    newGroParm = GromacsProperDihedralParameterInfo(atomtype1, atomtype2, atomtype3, atomtype4, func, phi, kphi, multiplicity, comment)

                if func == 2:
                    xi = float(fields[5]) * units.degrees
                    kxi = float(fields[6]) * units.kilojoules_per_mole * units.radians**(-2)
                    multiplicity = int(fields[7])
                    newGroParm = GromacsImproperDihedralParameterInfo(atomtype1, atomtype2, atomtype3, atomtype4, func, xi, kxi, comment)

                if func == 3:
                    C0 = float(fields[5]) * units.kilojoules_per_mole
                    C1 = float(fields[6]) * units.kilojoules_per_mole
                    C2 = float(fields[7]) * units.kilojoules_per_mole
                    C3 = float(fields[8]) * units.kilojoules_per_mole
                    C4 = float(fields[9]) * units.kilojoules_per_mole
                    C5 = float(fields[10]) * units.kilojoules_per_mole
                    newGroParm = GromacsRBDihedralParameterInfo(atomtype1, atomtype2, atomtype3, atomtype4, func, C0, C1, C2, C3, C4, C5, comment)

                if func == 5:
                    C1 = float(fields[5]) * units.kilojoules_per_mole
                    C2 = float(fields[6]) * units.kilojoules_per_mole
                    C3 = float(fields[7]) * units.kilojoules_per_mole
                    C4 = float(fields[8]) * units.kilojoules_per_mole
                    newGroParm = GromacsFourierDihedralParameterInfo(atomtype1, atomtype2, atomtype3, atomtype4, func, C1, C2, C3, C4, comment)

            if "exclusions" in parmDirective.getName():

                indices = list()
                atomtype1 = fields[0]
                for i in range(1, len(fields)-1):
                    indices.append(int(fields[i]))
                newGroParm = GromacsExclusionParameterInfo(atomtype1, indices, comment)

            if "constrainttypes" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                func = int(fields[2])

                if func == 1:
                    distance = float(fields[3]) * units.nanometers
                    newGroParm = GromacsConstraintParameterInfo(atomtype1, atomtype2, func, distance, comment)

                if func == 2:
                    distance = float(fields[3]) * units.nanometers
                    newGroParm = GromacsConstraintNCParameterInfo(atomtype1, atomtype2, func, distance, comment)

            if "settles" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                atomtype3 = fields[2]
                func = int(fields[3])

                if func == 1:
                    distanceOH = float(fields[4]) * units.nanometers
                    distanceHH = float(fields[5]) * units.nanometers
                    newGroParm = GromacsSettleParameterInfo(atomtype1, atomtype2, atomtype3, func, distanceOH, distanceHH, comment)

            if "vsites2" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                atomtype3 = fields[2]
                func = int(fields[3])

                if func == 1:
                    a = float(fields[4])
                    newGroParm = GromacsDummy2ParameterInfo(atomtype1, atomtype2, atomtype3, func, a, comment)

            if "vsites3" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                atomtype3 = fields[2]
                atomtype4 = fields[3]
                func = int(fields[4])

                if func == 1:
                    a = float(fields[5])
                    b = float(fields[6])
                    newGroParm = GromacsDummy3ParameterInfo(atomtype1, atomtype2, atomtype3, atomtype4, func, a, b, comment)

                if func == 2:
                    a = float(fields[5])
                    d = float(fields[6]) * units.nanometers
                    newGroParm = GromacsDummy3fdParameterInfo(atomtype1, atomtype2, atomtype3, atomtype4, func, a, d, comment)

                if func == 3:
                    theta = float(fields[5]) * units.degrees
                    d = float(fields[6]) * units.nanometers
                    newGroParm = GromacsDummy3fadParameterInfo(atomtype1, atomtype2, atomtype3, atomtype4, func, theta, d, comment)

                if func == 4:
                    a = float(fields[5])
                    b = float(fields[6])
                    c = float(fields[7]) * units.nanometers**(-1)
                    newGroParm = GromacsDummy3outParameterInfo(atomtype1, atomtype2, atomtype3, atomtype4, func, a, b, c, comment)

            if "vsites4" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                atomtype3 = fields[2]
                atomtype4 = fields[3]
                atomtype5 = fields[4]
                func = int(fields[5])

                if func == 1:
                    a = float(fields[6])
                    b = float(fields[7])
                    c = float(fields[8])
                    # Gromacs parameter needs to be implemented


            if "vsitesn" in parmDirective.getName():
                atomtype1 = fields[0]
                func = int(fields[1])

                if func == 1:
                    atomIndicies = int(fields[2]) # There can be more indicies
                    # Gromacs parameter needs to be implemented

                if func == 2:
                    atomIndicies = int(fields[2]) # There can be more indicies
                    # Gromacs parameter needs to be implemented

                if func == 3:
                    atomIndicies = int(fields[2]) # There can be more indicies
                    # Gromacs parameter needs to be implemented

            if "position_restraints" in parmDirective.getName():
                atomtype1 = fields[0]
                func = int(fields[1])

                if func == 1:
                    kx = float(fields[1]) * units.kilojoules_per_mole * units.nanometers**(-2)
                    ky = float(fields[2]) * units.kilojoules_per_mole * units.nanometers**(-2)
                    kz = float(fields[3]) * units.kilojoules_per_mole * units.nanometers**(-2)
                    newGroParm = GromacsPositionRestraintParameterInfo(atomtype1, func, kx, ky, kz, comment)

            if "distance_restraints" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                func = int(fields[2])

                if func == 1:
                    Type = fields[3]
                    label = fields[4]
                    low = float(fields[5]) * units.nanometers
                    up1 = float(fields[6]) * units.nanometers
                    up2 = float(fields[7]) * units.nanometers
                    weight = float(fields[8])
                    newGroParm = GromacsDistanceRestraintParameterInfo(atomtype1, atomtype2, func, Type, label, low, up1, up2, weight, comment)

            if "orientation_restraints" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                func = int(fields[2])

                if func == 1:
                    exp = fields[3]
                    label = fields[4]
                    alpha = float(fields[5])
                    c = float(fields[6])
                    obs = float(fields[7]) * units.amu
                    weight = float(fields[8]) * units.amu**(-1)
                    newGroParm = GromacsOrientationRestraintParameterInfo(atomtype1, atomtype2, func, exp, label, alpha, c, obs, weight, comment)

            if "angle_restraints" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                atomtype3 = fields[2]
                atomtype4 = fields[3]
                func = int(fields[4])

                if func == 1:
                    theta = float(fields[5]) * units.degrees
                    kc = float(fields[6]) * units.kilojoules_per_mole
                    multiplicity = int(fields[7])
                    newGroParm = GromacsAngleRestraintParameterInfo(atomtype1, atomtype2, atomtype3, atomtype4, func, theta, kc, multiplicity, comment)

            if "angle_restraints_z" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                func = int(fields[2])

                if func == 1:
                    theta = float(fields[3]) * units.degrees
                    kc = float(fields[4]) * units.kilojoules_per_mole
                    multiplicity = int(fields[5])
                    newGroParm = GromacsAngleRestraintParameterInfo(atomtype1, atomtype2, func, theta, kc, multiplicity, comment)

            if "nonbond_params" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                func = int(fields[2])

                if func == 1:
                    if self.parameters == 1 or self.parameters[0][0].cr == 3:
                        V = float(fields[3]) * units.kilojoules_per_mole * units.nanometers**6
                        W = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**12
                        newGroParm = GromacsLJNonbondedParameterInfo(atomtype1, atomtype2, func, V, W, comment)
                    elif self.parameters == 2:
                        V = float(fields[3]) * units.nanometers
                        W = float(fields[4]) * units.kilojoules_per_mole
                        newGroParm = GromacsLJNonbondedParameterInfo(atomtype1, atomtype2, func, V, W, comment)

                if func == 2:
                    a = float(fields[3]) * kilojoules_per_mole
                    b = float(fields[4]) * units.nanometers**(-1)
                    c6 = float(fields[5]) * units.kilojoules_per_mole * units.nanometers**6
                    newGroParm = GromacsBuckinghamNonbondedParameterInfo(atomtype1, atomtype2, func, a, b, c6, comment)

            parmList.append(newGroParm)
            index += 1
        return parmList


    def ConvertTopologyToDirectives(self, Topology):
        """
        Creating a list (by molecule) of lists (of Directives) for a single molecule.

        """

        myDirectivesList = list()

        # Create [ atoms ] directive from Topology.atoms
        atomsDirectiveName = '[ atoms ]'
        header = ';  nr  type    resnr   residue  atom    cgnr   charge     mass  typeB    chargeB      massB\n'
        atomsDirective = self.GromacsTopologyFileObject.Directive(atomsDirectiveName, header)

        exclusionDirectiveName = '[ exclusions ]'
        exclusionDirective = self.GromacsTopologyFileObject.Directive(exclusionDirectiveName, header='')

        for atom in Topology.atoms:
            convertedTopAtom = '%5d%6s%8d%8s%8s%8d%11.5f%10.3f\n'%(atom.ID, atom.metadata.atomtype, atom.metadata.resnum, atom.metadata.resname, atom.atomname, atom.metadata.cgnr, atom.charge._value, atom.mass._value)
            atomsDirective.lines.append(convertedTopAtom)

            # Test for exclusions
            if atom.exclusionList != list():
                exclusion = str(atom.ID)
                for excludedParticle in atom.exclusionList:
                    exclusion += '\t' + str(excludedParticle)
                exclusion += '\n'
                exclusionDirective.lines.append(exclusion)

        # Create [ moleculetype ] from the Topology
        moleculeTypeDirectiveName = '[ moleculetype ]'
        header = '; molname      nrexcl\n'
        moleculeTypeDirective = self.GromacsTopologyFileObject.Directive(moleculeTypeDirectiveName, header)
        molname = Topology.name
        nrexcl = Topology.nrexcl
        moleculeTypeDirective.lines.append('%s%9s\n'%(molname, nrexcl))
        myDirectivesList.append(moleculeTypeDirective)



        myDirectivesList.append(atomsDirective)

        # Create the remaining directives from Topology.forces
        for force in Topology.forces:
            forceInteractionType = type(force).__name__
            forceDirectiveType = ForceDictionary[forceInteractionType]
            forceDirectiveName = '[ %s ]'%(forceDirectiveType)
            forceDirective = self.GromacsTopologyFileObject.Directive(forceDirectiveName, header='')

            if forceDirectiveType == 'bonds':

                if forceInteractionType == 'BondForce':
                    forceDirective.header = '; i    j  func       b0          kb\n'
                    # Cycle through each bond
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 1
                        length = fields[2]._value
                        k = fields[3]._value
                        gromacsForce = '%(particle1)s   %(particle2)s       %(func)d    %(length)f   %(k)f\n'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'G96BondForce':
                    forceDirective.header = '; i    j  func       b0          kb\n'
                    # Cycle through each bond
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 2
                        length = fields[2]._value
                        k = fields[3]._value
                        gromacsForce = '%(particle1)s   %(particle2)s   %(func)d    %(length)f   %(k)f\n'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'MorseBondForce':
                    forceDirective.header = '; i    j  func       b0          D         beta\n'
                    # Cycle through each bond
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 3
                        length = fields[2]._value
                        D = fields[3]._value
                        beta = fields[4]._value
                        gromacsForce = '%(particle1)s   %(particle2)s   %(func)d    %(length)f   %(D)f   %(beta)f\n'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'CubicBondForce':
                    forceDirective.header = '; i    j  func       b0         C2         C3\n'
                    # Cycle through each bond
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 4
                        length = fields[2]._value
                        C2 = fields[3]._value
                        C3 = fields[4]._value
                        gromacsForce = '%(particle1)s   %(particle2)s   %(func)d    %(length)f   %(C2)f   %(C3)f\n'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'ConnectionBondForce':
                    forceDirective.header = '; i    j  func\n'
                    # Cycle through each bond
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 5
                        gromacsForce = '%(particle1)s   %(particle2)s   %(func)d\n'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'HarmonicBondForce':
                    forceDirective.header = '; i    j  func       b0          kb\n'
                    # Cycle through each bond
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 6
                        length = fields[2]._value
                        k = fields[3]._value
                        gromacsForce = '%(particle1)s   %(particle2)s   %(func)d    %(length)f   %(k)f\n'%vars()
                        forceDirective.lines.append(gromacsForce)

            if forceDirectiveType == 'pairs':

                if forceInteractionType == 'LJ1Force':
                    forceDirective.header = '; i    j  func       sigma         epsilon\n'
                    # Cycle through each pair
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 1
                        try:
                            V = str(fields[2]._value)
                            W = str(fields[3]._value)
                        except:
                            V = ''
                            W = ''
                        gromacsForce = '%(particle1)s   %(particle2)s   %(func)d    %(V)s   %(W)s\n'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'LJ2Force':
                    forceDirective.header = '; i    j  func       fudgeQQ     qi       qj       sigma         epsilon\n'
                    # Cycle through each pair
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 2
                        fudgeQQ = fields[2]._value
                        qi = fields[3]._value
                        qj = fields[4]._value
                        try:
                            V = str(fields[5]._value)
                            W = str(fields[6]._value)
                        except:
                            V = ''
                            W = ''
                        gromacsForce = '%(particle1)s   %(particle2)s   %(func)d    %(fudgeQQ)f   %(qi)f   %(qj)f   %(V)s   %(W)s\n'%vars()
                        forceDirective.lines.append(gromacsForce)

            if forceDirectiveType == 'pairs_nb':

                if forceInteractionType == 'LJNBForce':
                    forceDirective.header = '; i    j  func         qi       qj       sigma         epsilon\n'
                    # Cycle through each pair
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 1
                        qi = fields[2]._value
                        qj = fields[3]._value
                        try:
                            V = str(fields[4]._value)
                            W = str(fields[5]._value)
                        except:
                            V = ''
                            W = ''
                        gromacsForce = '%(particle1)s   %(particle2)s   %(func)d    %(qi)f   %(qj)f   %(V)s   %(W)s\n'%vars()
                        forceDirective.lines.append(gromacsForce)

            if forceDirectiveType == 'angles':

                if forceInteractionType == 'AngleForce':
                    forceDirective.header = ';  i    j    k  func       th0       cth\n'
                    # Cycle through each angles
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        func = 1
                        theta = fields[3]._value
                        k = fields[4]._value
                        gromacsForce = '%(particle1)s   %(particle2)s   %(particle3)s   %(func)d    %(theta)f   %(k)f\n'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'G96AngleForce':
                    forceDirective.header = ';  i    j    k  func       th0       cth\n'
                    # Cycle through each angles
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        func = 2
                        theta = fields[3]._value
                        k = fields[4]._value
                        gromacsForce = '%(particle1)s   %(particle2)s   %(particle3)s   %(func)d    %(theta)f   %(k)f\n'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'CrossBondBondAngleForce':
                    forceDirective.header = ';  i    j    k  func       r1e       r2e         cth\n'
                    # Cycle through each angles
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        func = 3
                        r1e = fields[3]._value
                        r2e = fields[4]._value
                        krr = fields[5]._value
                        gromacsForce = '%(particle1)s   %(particle2)s   %(particle3)s   %(func)d    %(r1e)f   %(r2e)f   %(krr)f\n'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'CrossBondAngleAngleForce':
                    forceDirective.header = ';  i    j    k  func       r1e       r2e         r3e        cth\n'
                    # Cycle through each angles
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        func = 4
                        r1e = fields[3]._value
                        r2e = fields[4]._value
                        r3e = fields[5]._value
                        kr = fields[6]._value
                        gromacsForce = '%(particle1)s   %(particle2)s   %(particle3)s   %(func)d    %(r1e)f   %(r2e)f   %(r3e)f   %(kr)f\n'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'UreyBradleyAngleForce':
                    forceDirective.header = ';  i    j    k  func       th0       cth         r13        kUB\n'
                    # Cycle through each angles
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        func = 5
                        theta = fields[3]._value
                        k = fields[4]._value
                        r13 = fields[5]._value
                        kUB = fields[6]._value
                        gromacsForce = '%(particle1)s   %(particle2)s   %(particle3)s   %(func)d    %(theta)f   %(k)f   %(r13)f   %(kUB)f\n'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'QuarticAngleForce':
                    forceDirective.header = ';  i    j    k  func       th0        c0         c1         c2         c3         c4\n'
                    # Cycle through each angles
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        func = 6
                        theta = fields[3]._value
                        C0 = fields[4]._value
                        C1 = fields[5]._value
                        C2 = fields[6]._value
                        C3 = fields[7]._value
                        C4 = fields[8]._value
                        gromacsForce = '%(particle1)s   %(particle2)s   %(particle3)s   %(func)d    %(theta)f   %(C0)f   %(C1)f   %(C2)f   %(C3)f   %(C4)f\n'%vars()
                        forceDirective.lines.append(gromacsForce)

            if forceDirectiveType == 'dihedrals':

                if forceInteractionType == 'PeriodicTorsionForce':
                    forceDirective.header = ';  ai  aj      ak      al   funct     th0         kth0       mult\n'
                    # Cycle through each torsion
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        particle4 = fields[3]
                        func = 1
                        phi = fields[5]._value
                        kphi = fields[6]._value
                        multiplicity = fields[4]
                        gromacsForce = '%(particle1)s   %(particle2)s   %(particle3)s   %(particle4)s   %(func)d    %(phi)f   %(kphi)f   %(multiplicity)d\n'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'NonPeriodicTorsionForce':
                    forceDirective.header = ';  ai  aj      ak      al   funct     xi         kxi\n'
                    # Cycle through each torsion
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        particle4 = fields[3]
                        func = 2
                        xi = fields[4]._value
                        kxi = fields[5]._value
                        gromacsForce = '%(particle1)s   %(particle2)s   %(particle3)s   %(particle4)s   %(func)d    %(xi)f   %(kxi)f\n'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'RBTorsionForce':
                    forceDirective.header = ';   ai    aj    ak    al funct       c0          c1          c2          c3          c4          c5\n'
                    # Cycle through each torsion
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        particle4 = fields[3]
                        func = 3
                        C0 = fields[4]._value
                        C1 = fields[5]._value
                        C2 = fields[6]._value
                        C3 = fields[7]._value
                        C4 = fields[8]._value
                        C5 = fields[9]._value
                        gromacsForce = '%(particle1)s   %(particle2)s   %(particle3)s   %(particle4)s   %(func)d    %(C0)f   %(C1)f   %(C2)f   %(C3)f   %(C4)f   %(C5)f\n'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'FourierTorsionForce':
                    forceDirective.header = ';   ai    aj    ak    al funct       c0          c1          c2          c3          c4\n'
                    # Cycle through each torsion
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        particle4 = fields[3]
                        func = 5
                        C1 = fields[4]._value
                        C2 = fields[5]._value
                        C3 = fields[6]._value
                        C4 = fields[7]._value
                        gromacsForce = '%(particle1)s   %(particle2)s   %(particle3)s   %(particle4)s   %(func)d    %(C1)f   %(C2)f   %(C3)f   %(C4)f\n'%vars()
                        forceDirective.lines.append(gromacsForce)

            if forceDirectiveType == 'settles':
                if forceInteractionType == 'SettleForce':
                    forceDirective.header = '; i    funct   length1   length2'
                    # Cycle through each settle bond
                    for i in range(force.getNumForces()):
                        fields = force.getForceParameters(i)
                        particle = fields[0]
                        func = 1
                        length1 = fields[1]._value
                        length2 = fields[2]._value
                        gromacsForce = '%(particle)s   %(func)d   %(length1)f   %(length2)f'%vars()

            myDirectivesList.append(forceDirective)

        if exclusionDirective.lines != []:
            myDirectivesList.append(exclusionDirective)

        return myDirectivesList

    def CreateDefaultsDirective(self):
        directiveName = '[ defaults ]'
        directiveHeader = '; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n'
        directive = self.GromacsTopologyFileObject.Directive(directiveName, directiveHeader)

        for mol in self.molecules:
            try:
                nbfunc = mol.NBFunc
                combrule = mol.CombinationRule
                genpairs = 'yes'
                fudgeLJ = mol.LJCorrection
                fudgeQQ = mol.CoulombCorrection
                directive.lines.append('%(nbfunc)d      %(combrule)d      %(genpairs)s      %(fudgeLJ)f      %(fudgeQQ)f\n'%vars())
		return directive
            except:
                pass

        return "Error: No default parameters found" 


    def CreateAtomtypesDirective(self):
        """

        """
        directiveName = '[ atomtypes ]'
        directiveHeader = ';name     mass    charge   ptype          sigma      epsilon\n'
        directive = self.GromacsTopologyFileObject.Directive(directiveName, directiveHeader)

        for mol in self.molecules:
            for atom in mol.atoms:
                name = atom.metadata.atomtype
                mass = 0.0000
                charge = 0.0000
                ptype = 'A'
                sigma = atom.sigma._value
                epsilon = atom.epsilon._value
                line = '%(name)s  %(mass)f   %(charge)f   %(ptype)s   %(sigma)f   %(epsilon)f\n'%vars()
                if line not in directive.lines:
                    directive.lines.append(line)
        directive.lines.sort()
        directiveList = list()
        directiveList.append(directive)
        return directiveList

    def CreateSystemDirective(self, molecules):
        """

        """

        myDirectivesList = []

        # Create [ system ]
        systemDirectiveName = '[ system ]'
        systemDirective = self.GromacsTopologyFileObject.Directive(systemDirectiveName, header='')
        ### Needs to be fixed ###
        systemName = 'system\n'
        #########################
        systemDirective.lines.append(systemName)
        myDirectivesList.append(systemDirective)

        # Create [ molecules ] section
        moleculeDirectiveName = '[ molecules ]'
        header = '; Compound        nmols\n'
        moleculeDirective = self.GromacsTopologyFileObject.Directive(moleculeDirectiveName, header)
        for molecule in molecules:
            moleculeName = molecule.name
            nmolecules = len(molecule.atoms)
            moleculeLine = '%s %16d\n'%(moleculeName, nmolecules)
            moleculeDirective.lines.append(moleculeLine)
        myDirectivesList.append(moleculeDirective)

        return myDirectivesList


    def writeTopologyFile(self, filename):
        """Write the Gromacs topology *.top file to file.

        """

        # Process filename
        if filename[-4:] != '.top':
            filename += '.top'

        # Rebuild the Directives from the information in the Topology objects.

        self.GromacsTopologyFileObject = GromacsTopologyParser()

        # Create [ default ] and [ atomtypes ]
        self.GromacsTopologyFileObject.directives.append( self.CreateDefaultsDirective() )
        self.GromacsTopologyFileObject.directives.extend( self.CreateAtomtypesDirective() )

        # Translate each set of Topology objects to molecule directives
        for mol in range(len(self.molecules)):
            self.GromacsTopologyFileObject.directives.extend( self.ConvertTopologyToDirectives(self.molecules[mol] ) )

        # Create System directives
        self.GromacsTopologyFileObject.directives.extend( self.CreateSystemDirective(self.molecules) )

        self.GromacsTopologyFileObject.write(filename)

        return

