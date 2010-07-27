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
"FourierTorsionForce":"dihedrals"
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

        self.parameters = list() # These are GromacsParameter[i] objects, akin to Force objects,
                                 # but holding Force parameters for various atomtype combinations 
        self.molecules = list()  # molecules[i] is the ith Topology object
        self.index = list()

        if (topfile):
            self.GromacsTopologyFileObject = GromacsTopologyFile(topfile=topfile)
        else:
            self.GromacsTopologyFileObject = None
        return

    def setName(self, name):
        self.name = name
        return

    def readTopologyFile(self, topfile):
        """Read in the contents of the topology file.
        NOTES

        """

        # Read in the *.top file, creating a new GromacsTopologyFile object
        self.GromacsTopologyFileObject = GromacsTopologyFile(topfile=topfile)

        parmIndex = list()
        # Translate all the parameter directives to GromacsParameter objects
        for parm in self.GromacsTopologyFileObject.ParameterDirectives:
            self.parameters.append(self.ConvertParameterDirectiveToGromacsParameter(parm))
            parmIndex.append(parm.getName())
        self.index.append(parmIndex)

        #incredibly sloppy function that combines the contents of duplicate parameter entries
        set = {}
        copy = [set.setdefault(e,e) for e in self.index[0] if e not in set]

        reallist = list()
        for i in range(len(self.index[0])):
            subtemp = list()
            for j in range(i, len(self.index[0])):
                if self.index[0][i] == self.index[0][j]:
                    parmInfo = self.parameters[j]
                    subtemp.extend(parmInfo)
            reallist.append(subtemp)

        self.index[0] = copy
        while len(reallist) > len(self.index[0]):
            reallist.pop()
        self.parameters = reallist

        # Translate each set of molecule directives to Topology objects
        for mol in range(len(self.GromacsTopologyFileObject.MoleculeDefinitionDirectives)):
            self.molecules.append( self.ConvertMoleculeDirectivesToTopology( self.GromacsTopologyFileObject.MoleculeDefinitionDirectives[mol] ) )

        return

    def grabMissingInfo(self, molecule, header, func, particle1, particle2, particle3 = None, particle4 = None):

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
            i = self.index[0].index('[ constrainttypes ]')

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

    def findMoleculeDirective(self, currentMoleculeDirective, molecule, forceType):
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
            if forceType == "bonds":
                if currentMoleculeDirective.getNumBonds() == 0:
                    molecule.forces.append(currentMoleculeDirective)
            if forceType == "pairs":
                if currentMoleculeDirective.getNumForces() == 0:
                    molecule.forces.append(currentMoleculeDirective)
            if forceType == "pairs_nb":
                if currentMoleculeDirective.getNumForces() == 0:
                    molecule.forces.append(currentMoleculeDirective)
            if forceType == "angles":
                if currentMoleculeDirective.getNumAngles() == 0:
                    molecule.forces.append(currentMoleculeDirective)
            if forceType == "dihedrals":
                if currentMoleculeDirective.getNumTorsions() == 0:
                    molecule.forces.append(currentMoleculeDirective)
            if forceType == "exclusions":
                raise "Not Implemented Yet"
                """
                if currentMoleculeDirective.getNumForces() == 0:
                    molecule.forces.append(currentMoleculeDirective)
                """
            if forceType == "constraints":
                raise "Not Implemented Yet"
                """
                if currentMoleculeDirective.getNumForces() == 0:
                    molecule.forces.append(currentMoleculeDirective)
                """
            if forceType == "settles":
                raise "Not Implemented Yet"
                """
                if currentMoleculeDirective.getNumForces() == 0:
                    molecule.forces.append(currentMoleculeDirective)
                """
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

        # molecule has name, atoms (list), forces (list), constraints (list), atomgroups (list)
        for directive in MoleculeDirectives:
            if len(directive.lines) > 0:
                # Fill in MoleculeType
                if "moleculetype" in directive.name:
                    for line in directive.lines:
                        fields = line.split()
                        molname = fields[0]
                        nrexcl = fields[1]
                        molecule.name = molname

                # Fill in atoms
                if "atoms" in directive.name:
                    for line in directive.lines:
                        fields, comment = self.splitline(line)
                        try: particle = int(fields[0])
                        except: pdb.set_trace()
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
                            mass = 0.0000 * units.amu
                        atomnum = particle
                        freeEnergyAtom = False
                        for atomtypenis in self.parameters[1]:
                            if atomtype == atomtypenis.name:
                                V = atomtypenis.sigma
                                W = atomtypenis.epsilon
                        atom = TopAtom(particle, atomtype, resnum, resname, atomname, cgnr, charge, mass, V, W, freeEnergyAtom)
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
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "bonds")
                            if len(fields) > 3:
                                length = float(fields[3]) * units.nanometers
                                k = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**(-2)
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                length = self.parameters[i][j].distance
                                k = self.parameters[i][j].kspring
                            currentMoleculeDirective.addBond(particle1, particle2, length, k)

                        if func == 2:
                            raise "Bond type not implemented"
                            """
                            newMoleculeDirective = G96BondForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "bonds")
                            if len(fields) > 3:
                                length = float(fields[3]) * units.nanometers
                                k = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**(-4)
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                length = self.parameters[i][j].distance
                                k = self.parameters[i][j].kspring
                            currentMoleculeDirective.addBond(particle1, particle2, length, k)
                            """

                        if func == 3:
                            raise "Bond type not implemented"
                            """
                            newMoleculeDirective = MorseBondForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "bonds")
                            if len(fields) > 3:
                                length = float(fields[3]) * units.nanometers
                                D = float(fields[4]) * units.kilojoules_per_mole
                                beta = float(fields[5]) * units.nanometers**(-1)
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                D = self.parameters[i][j].D
                                beta = self.parameters[i][j].beta
                            currentMoleculeDirective.addBond(particle1, particle2, length, D, beta)
                            """

                        if func == 4:
                            raise "Bond type not implemented"
                            """
                            newMoleculeDirective = CubicBondForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "bonds")
                            if len(fields) > 3:
                                length = float(fields[3]) * units.nanometers
                                C2 = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**(-2)
                                C3 = float(fields[5]) * units.kilojoules_per_mole * units.nanometers**(-3)
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                C2 = self.parameters[i][j].C2
                                C3 = self.parameters[i][j].C3
                            currentMoleculeDirective.addBond(particle1, particle2, length, C2, C3)
                            """

                        if func == 5:
                            raise "Bond type not implemented"
                            """
                            newMoleculeDirective = ConnectionBondForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "bonds")
                            currentMoleculeDirective.addBond(particle1, particle2)
                            """

                        if func == 6:
                            raise "Bond type not implemented"
                            """
                            newMoleculeDirective = HarmonicBondForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "bonds")
                            if len(fields) > 3:
                                length = float(fields[3]) * units.nanometers
                                k = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**(-2)
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                length = self.parameters[i][j].distance
                                k = self.parameters[i][j].kspring
                            currentMoleculeDirective.addBond(particle1, particle2, length, k)
                            """

                        if func == 7:
                            raise "Bond type not implemented"
                            # Not implemented in GromacsParameter.py

                        if func == 8:
                            raise "Bond type not implemented"
                            # Not implemented in GromacsParameter.py

                        if func == 9:
                            raise "Bond type not implemented"
                            # Not implemented in GromacsParameter.py

                # Fill in pairs
                if "pairs" in directive.name:
                        fields, comment = self.splitline(line)
                        atomtype1 = int(fields[0])
                        atomtype2 = int(fields[1])
                        func = int(fields[2])

                        if func == 1:
                            newMoleculeDirective = LJ1Force() # depends on combination rule
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "pairs")
                            if len(fields) > 3:
                                if self.parameters[0][0].cr == 1:
                                    V = float(fields[3]) * units.kilojoules_per_mole * units.nanometers**6
                                    W = float(fields[4]) * units.kilojoules_per_mole * units.nanometers**12
                                elif self.parameters[0][0].cr == 2 or self.parameters[0][0].cr == 3:
                                    V = float(fields[3]) * units.nanometers**6
                                    W = float(fields[4]) * units.kilojoules_per_mole
                            elif "pairtypes" in self.index[0]:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                V = self.parameters[i][j].V
                                W = self.parameters[i][j].W
                            else:
                                V = None
                                W = None
                            currentMoleculeDirective.addForce(atomtype1, atomtype2, V, W) # Depends on combination rule

                        if func == 2:
                            newMoleculeDirective = LJ2Force()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "pairs")
                            if len(fields) > 3:
                                fudgeQQ = fields[3]
                                qi = float(fields[4]) * units.elementary_charge
                                qj = float(fields[5]) * units.elementary_charge
                                if self.parameters[0][0].cr == 1:
                                    V = float(fields[6]) * units.nanometers**6
                                    W = float(fields[7]) * units.kilojoules_per_mole
                                elif self.parameters[0][0].cr == 2 or self.parameters[0][0].cr == 3:
                                    V = float(fields[6]) * units.kilojoules_per_mole * units.nanometers**6
                                    W = float(fields[7]) * units.kilojoules_per_mole * units-nanometers**12
                            elif "pairtypes" in self.index[0]:
                                i, j = grabMissingInfo(molecule, directive.name, func, particle1, particle2)
                                V = self.parameters[i][j].V
                                W = self.parameters[i][j].W
                            else:
                                V = None
                                W = None
                            currentMoleculeDirective.addForce(atomtype1, atomtype2, qi, qj, V, W)

                # Fill in nonbonded pairs
                if "pairs_nb" in directive.name:
                        fields, comment = self.splitline(line)
                        func = int(fields[2])

                        if func == 1:
                            newMoleculeDirective = LGNBForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "pairs_nb")
                            # Need to finish

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
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "angles")
                            if len(fields) > 4:
                                angle = float(fields[4]) * units.degrees
                                k = float(fields[5]) * units.kilojoules_per_mole * units.radians**(-2)
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3)
                                angle = self.parameters[i][j].theta
                                k = self.parameters[i][j].k
                            currentMoleculeDirective.addAngle(particle1, particle2, particle3, angle, k)

                        if func == 2:
                            raise "Angle type not implemented"
                            """
                            newMoleculeDirective = G96AngleForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "angles")
                            if len(fields) > 4:
                                angle = float(fields[4]) * units.degrees
                                k = float(fields[5]) * units.kilojoules_per_mole
                            else:
                                i, j = grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3)
                                angle = self.parameters[i][j].theta
                                k = self.parameters[i][j].k
                            currentMoleculeDirective.addAngle(particle1, particle2, particle3, angle, k)
                            """

                        if func == 3:
                            raise "Angle type not implemented"
                            """
                            newMoleculeDirective = CrossBondBondAngleForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "angles")
                            if len(fields) > 4:
                                r1e = float(fields[4]) * units.nanometers
                                r2e = float(fields[5]) * units.nanometers
                                krr = float(fields[6]) * units.kilojoules_per_mole * units.nanometers**(-2)
                            else:
                                i, j = grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3)
                                r1e = self.parameters[i][j].r1e
                                r2e = self.parameters[i][j].r2e
                                krr = self.parameters[i][j].krr
                            currentMoleculeDirective.addAngle(particle1, particle2, particle3, r1e, r2e, krr)
                            """

                        if func == 4:
                            raise "Angle type not implemented"
                            """
                            newMoleculeDirective = CrossBondAngleAngleForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "angles")
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
                            currentMoleculeDirective.addAngle(particle1, particle2, particle3, r1e, r2e, r3e, kr)
                            """

                        if func == 5:
                            raise "Angle type not implemented"
                            """
                            newMoleculeDirective = UreyBradleyAngleForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "angles")
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
                            currentMoleculeDirective.addAngle(particle1, particle2, particle3, angle, k, r13, kUB)
                            """

                        if func == 6:
                            raise "Angle type not implemented"
                            """
                            newMoleculeDirective = QuarticAngleForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "angles")
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
                            currentMoleculeDirective.addAngle(particle1, particle2, particle3, angle, C0, C1, C2, C3, C4)
                            """

                        if func == 8:
                            raise "Angle type not implemented"
                            # Not implemented in GromacsParameter.py

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
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "dihedrals")
                            if len(fields) > 5:
                                phi = float(fields[5]) * units.degrees
                                kphi = float(fields[6]) * units.kilojoules_per_mole
                                multiplicity = int(fields[7])
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3, particle4)
                                phi = self.parameters[i][j].phi
                                kphi = self.parameters[i][j].kphi
                                multiplicity = self.parameters[i][j].multiplicity
                            currentMoleculeDirective.addTorsion(particle1, particle2, particle3, particle4, multiplicity, phi, kphi)

                        if func == 2:
                            newMoleculeDirective = NonPeriodicTorsionForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "dihedrals")
                            if len(fields) > 5:
                                xi = float(fields[5]) * units.degrees
                                kxi = float(fields[6]) * units.kilojoules_per_mole * units.radians**(-2)
                            else:
                                i, j = self.grabMissingInfo(molecule, directive.name, func, particle1, particle2, particle3, particle4)
                                xi = self.parameters[i][j].xi
                                kxi = self.parameters[i][j].kxi
                            currentMoleculeDirective.addTorsion(particle1, particle2, particle3, particle4, xi, kxi)

                        if func == 3:
                            newMoleculeDirective = RBTorsionForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "dihedrals")
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
                            currentMoleculeDirective.addTorsion(particle1, particle2, particle3, particle4, C0, C1, C2, C3, C4, C5)

                        if func == 5:
                            newMoleculeDirective = FourierTorsionForce()
                            currentMoleculeDirective = self.findMoleculeDirective(newMoleculeDirective, molecule, "dihedrals")
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
                            currentMoleculeDirective.addTorsion(particle1, particle2, particle3, particle4, C0, C1, C2, C3, C4)

                        if func == 8:
                            raise "Dihedral type not implemented"

                # Fill in exclusions
                if "exclusions" in directive.name:
                    for line in directive.lines:
                        fields, comment = self.splitline(line)
                        # Need to finish

                # Fill in constraints
                if "constraints" in directive.name:
                    for line in directive.lines:
                        fields, comment = self.splitline(line)
                        func = int(fields[4])

                        if func == 1:
                            raise "Constraint type not implemented"

                        if func == 2:
                            raise "Constraint type not implemented"

                # Fill in settles
                if "settles" in directive.name:
                    for line in directive.lines:
                        fields, comment = self.splitline(line)
                        # Need to finish

            molIndex.append(directive.getName())
        self.index.append(molIndex)

        return molecule


    def ConvertParameterDirectiveToGromacsParameter(self, parmDirective):
        """Takes a parmDirective (Directive object), and uses it
        to create a GromacsParameter object to store the forcefield parameters."""

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
                if len(fields) == 8:
                    i += 2
                    bondtype = fields[i-1]
                    atomicmass = fields[i]
                elif len(fields) == 7:
                    try:
                        int(fields[i+1])
                        i += 1
                        atomicmass = fields[i+1]
                    except:
                        i += 1
                        bondtype = fields[i]
                else:
                    bondtype = ''
                    atomicmass = 0.0000
                mass = float(fields[i+1]) * units.amu
                charge = float(fields[i+2]) * units.elementary_charge
                ptype = fields[i+3]
                if self.parameters[0][0].cr == 1:
                    sigma = float(fields[i+4]) * units.kilojoules_per_mole*units.nanometers**6
                    epsilon = float(fields[i+5]) * units.kilojoules_per_mole*units.nanometers**12
                    newGroParm = GromacsAtom1ParameterInfo(name, bondtype, mass, charge, ptype, sigma, epsilon, comment)
                elif self.parameters[0][0].cr == 2 or self.parameters[0][0].cr == 3:
                    sigma = float(fields[i+4]) * units.nanometers
                    epsilon = float(fields[i+5]) * units.kilojoules_per_mole
                    newGroParm = GromacsAtom23ParameterInfo(name, bondtype, mass, charge, ptype, sigma, epsilon, comment)

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
                        V = float(fields[3]) * units.kilojoules_per_mole*units.nanometers**6
                        W = float(fields[4]) * units.kilojoules_per_mole*units.nanometers**12
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
                        V = float(fields[6]) * units*kilojoules_per_mole*units.nanometers**6
                        W = float(fields[7]) * units.kilojoules_per_mole*units.nanometers**12
                        newGroParm = GromacsPairTypes2ParameterInfo(atomtype1, atomtype2, func, fudgeQQ, qi, qj, V, W, comment)

            if "pairs_nb" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                func = int(fields[2])
                qi = float(fields[3]) * units.elementary_charge
                qj = float(fields[4]) * units.elementary_charge

                if func == 1:
                    if self.parameters[0][0] == 1 or self.parameters[0][0].cr == 3:
                        V = float(fields[5]) * units.kilojoules_per_mole*units.nanometers**6
                        W = float(fields[6]) * units.kilojoules_per_mole*units.nanometers**12
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
                    kr = float(fields[7]) * units.kilojoules_per_mole*units.nanometers**(-2)
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
                    kxi = float(fields[6]) * units.kilojoules_per_mole*units.radians**(-2)
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

            if "virtual_sites2" in parmDirective.getName():
                atomtype1 = fields[0]
                atomtype2 = fields[1]
                atomtype3 = fields[2]
                func = int(fields[3])

                if func == 1:
                    a = float(fields[4])
                    newGroParm = GromacsDummy2ParameterInfo(atomtype1, atomtype2, atomtype3, func, a, comment)

            if "virtual_sites3" in parmDirective.getName():
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

            if "position_restraints" in parmDirective.getName():
                atomtype1 = fields[0]
                func = int(fields[1])

                if func == 1:
                    kx = float(fields[1]) * units.kilojoules_per_mole*units.nanometers**(-2)
                    ky = float(fields[2]) * units.kilojoules_per_mole*units.nanometers**(-2)
                    kz = float(fields[3]) * units.kilojoules_per_mole*units.nanometers**(-2)
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
                        V = float(fields[3]) * units.kilojoules_per_mole*units.nanometers**6
                        W = float(fields[4]) * units.kilojoules_per_mole*units.nanometers**12
                        newGroParm = GromacsLJNonbondedParameterInfo(atomtype1, atomtype2, func, V, W, comment)
                    elif self.parameters == 2:
                        V = float(fields[3]) * units.nanometers
                        W = float(fields[4]) * units.kilojoules_per_mole
                        newGroParm = GromacsLJNonbondedParameterInfo(atomtype1, atomtype2, func, V, W, comment)

                if func == 2:
                    a = float(fields[3]) * kilojoules_per_mole
                    b = float(fields[4]) * units.nanometers**(-1)
                    c6 = float(fields[5]) * units.kilojoules_per_mole*units.nanometers**6
                    newGroParm = GromacsBuckinghamNonbondedParameterInfo(atomtype1, atomtype2, func, a, b, c6, comment)

            parmList.append(newGroParm)
            index += 1
        return parmList


    def ConvertTopologyToMoleculeDirectives(self, Topology):
        """
        Creating a list (by molecule) of lists (of Directives) for a single molecule

        """

        myDirectiveList = list()

        # Create [ moleculetype ] from the Topology
        moleculeTypeDirectiveName = '[ moleculetype ]\n'
        header = '; molname      nrexcl\n'
        moleculeTypeDirective = self.GromacsTopologyFileObject.Directive(moleculeTypeDirectiveName, header)
        molname = Topology.name
        nrexcl = '' # Needs to come from the number of exclusions
        moleculeTypeDirective.lines.append('%s%9s\n'%(molname, nrexcl))
        myDirectiveList.append(moleculeTypeDirective)

        # Create [ atoms ] directive from Topology.atoms
        atomsDirectiveName = '[ atoms ]\n'
        header = ';  nr  type    resnr   residue  atom    cgnr   charge     mass  typeB    chargeB      massB\n'
        atomsDirective = self.GromacsTopologyFileObject.Directive(atomsDirectiveName, header)
        for atom in Topology.atoms:
            convertedTopAtom = '%5d%6s%8d%8s%8s%8d%11.5f%10.3f\n'%(atom.ID, atom.metadata.atomtype, atom.metadata.resnum, atom.metadata.resname, atom.atomname, atom.metadata.cgnr, atom.charge._value, atom.mass._value)
            atomsDirective.lines.append(convertedTopAtom)
        myDirectiveList.append(atomsDirective)

        # Create the remaining directives from Topology.forces
        for force in Topology.forces:
            forceInteractionType = type(force).__name__
            forceDirectiveType = ForceDirectionary[forceInteractionType]
            forceDirectiveName = '[ %s ]\n'%(forceInteractionType)
            forceDirective = self.GromacsTopologyFileObject.Directive(forceDirectiveName, header='')

            if forceType == 'bonds':

                if forceInteractionType == 'BondForce':
                    forceDirective.header = '; i    j  func       b0          kb'
                    # Cycle through each bond
                    for i in range(len(force.getNumBonds())):
                        fields = force.getBondParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 1
                        length = fields[2]._value
                        k = fields[3]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(func)d%(length)f%(k)f'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'G96BondForce':
                    forceDirective.header = '; i    j  func       b0          kb'
                    # Cycle through each bond
                    for i in range(len(force.getNumBonds())):
                        fields = force.getBondParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 2
                        length = fields[2]._value
                        k = fields[3]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(func)d%(length)f%(k)f'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'MorseBondForce':
                    forceDirective.header = '; i    j  func       b0          D         beta'
                    # Cycle through each bond
                    for i in range(len(force.getNumBonds())):
                        fields = force.getBondParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 3
                        length = fields[2]._value
                        D = fields[3]._value
                        beta = fields[4]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(func)d%(length)f%(D)f%(beta)f'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'CubicBondForce':
                    forceDirective.header = '; i    j  func       b0         C2         C3'
                    # Cycle through each bond
                    for i in range(len(force.getNumBonds())):
                        fields = force.getBondParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 4
                        length = fields[2]._value
                        C2 = fields[3]._value
                        C3 = fields[4]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(func)d%(length)f%(C2)f%(C3)f'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'ConnectionBondForce':
                    forceDirective.header = '; i    j  func'
                    # Cycle through each bond
                    for i in range(len(force.getNumBonds())):
                        fields = force.getBondParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 5
                        gromacsForce = '%(particle1)s%(particle2)s%(func)d'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'HarmonicBondForce':
                    forceDirective.header = '; i    j  func       b0          kb'
                    # Cycle through each bond
                    for i in range(len(force.getNumBonds())):
                        fields = force.getBondParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 6
                        length = fields[2]._value
                        k = fields[3]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(func)d%(length)f%(k)f'%vars()
                        forceDirective.lines.append(gromacsForce)

            if forceType == 'pairs':

                if forceInteractionType == 'LJ1Force':
                    forceDirective.header = '; i    j  func       sigma         epsilon'
                    # Cycle through each pair
                    for i in range(len(force.getNumForces())):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 1
                        V = fields[2]._value
                        W = fields[3]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(func)d%(V)f%(W)f'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'LJ2Force':
                    forceDirective.header = '; i    j  func       fudgeQQ     qi       qj       sigma         epsilon'
                    # Cycle through each pair
                    for i in range(len(force.getNumForcess())):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 2
                        fudgeQQ = fields[2]._value
                        qi = fields[3]._value
                        qj = fields[4]._value
                        V = fields[5]._value
                        W = fields[6]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(func)d%(fudgeQQ)f%(qi)f%(qj)f%(V)f%(W)f'%vars()
                        forceDirective.lines.append(gromacsForce)

            if forceType == 'pairs_nb':

                if forceInteractionType == 'LJNBForce':
                    forceDirective.header = '; i    j  func         qi       qj       sigma         epsilon'
                    # Cycle through each pair
                    for i in range(len(force.getNumForcess())):
                        fields = force.getForceParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        func = 1
                        qi = fields[2]._value
                        qj = fields[3]._value
                        V = fields[4]._value
                        W = fields[5]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(func)d%(qi)f%(qj)f%(V)f%(W)f'%vars()
                        forceDirective.lines.append(gromacsForce)

            if forceType == 'angles':

                if forceInteractionType == 'AngleForce':
                    forceDirective.header = ';  i    j    k  func       th0       cth'
                    # Cycle through each angles
                    for i in range(len(force.getNumAngless())):
                        fields = force.getAngleParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        func = 1
                        theta = fields[3]._value
                        k = fields[4]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(particle3)s%(func)d%(theta)f%(k)f'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'G96AngleForce':
                    forceDirective.header = ';  i    j    k  func       th0       cth'
                    # Cycle through each angles
                    for i in range(len(force.getNumAngless())):
                        fields = force.getAngleParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        func = 2
                        theta = fields[3]._value
                        k = fields[4]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(particle3)s%(func)d%(theta)f%(k)f'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'CrossBondBondAngleForce':
                    forceDirective.header = ';  i    j    k  func       r1e       r2e         cth'
                    # Cycle through each angles
                    for i in range(len(force.getNumAngless())):
                        fields = force.getAngleParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        func = 3
                        r1e = fields[3]._value
                        r2e = fields[4]._value
                        krr = fields[5]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(particle3)s%(func)d%(r1e)f%(r2e)f%(krr)f'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'CrossBondAngleAngleForce':
                    forceDirective.header = ';  i    j    k  func       r1e       r2e         r3e        cth'
                    # Cycle through each angles
                    for i in range(len(force.getNumAngless())):
                        fields = force.getAngleParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        func = 4
                        r1e = fields[3]._value
                        r2e = fields[4]._value
                        r3e = fields[5]._value
                        kr = fields[6]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(particle3)s%(func)d%(r1e)f%(r2e)f%(r3e)f%(kr)f'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'UreyBradleyAngleForce':
                    forceDirective.header = ';  i    j    k  func       th0       cth         r13        kUB'
                    # Cycle through each angles
                    for i in range(len(force.getNumAngless())):
                        fields = force.getAngleParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        func = 5
                        theta = fields[3]._value
                        k = fields[4]._value
                        r13 = fields[5]._value
                        kUB = fields[6]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(particle3)s%(func)d%(theta)f%(k)f%(r13)f%(kUB)f'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'QuarticAngleForce':
                    forceDirective.header = ';  i    j    k  func       th0        c0         c1         c2         c3         c4'
                    # Cycle through each angles
                    for i in range(len(force.getNumAngless())):
                        fields = force.getAngleParameters(i)
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
                        gromacsForce = '%(particle1)s%(particle2)s%(particle3)s%(func)d%(theta)f%(C0)f%(C1)f%(C2)f%(C3)f%(C4)f'%vars()
                        forceDirective.lines.append(gromacsForce)

            if forceType == 'dihedrals':

                if forceInteractionType == 'PeriodicTorsionForce':
                    forceDirective.header = ';  ai  aj      ak      al   funct     th0         kth0       mult'
                    # Cycle through each torsion
                    for i in range(len(force.getNumTorsions())):
                        fields = force.getTorsionParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        particle4 = fields[3]
                        func = 1
                        phi = fields[4]._value
                        kphi = fields[5]._value
                        multiplicity = fields[6]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(particle3)s%(particle4)s%(func)d%(phi)f%(kphi)f%(multiplicity)d'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'NonPeriodicTorsionForce':
                    forceDirective.header = ';  ai  aj      ak      al   funct     xi         kxi'
                    # Cycle through each torsion
                    for i in range(len(force.getNumTorsions())):
                        fields = force.getTorsionParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        particle4 = fields[3]
                        func = 2
                        xi = fields[4]._value
                        kxi = fields[5]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(particle3)s%(particle4)s%(func)d%(xi)f%(kxi)f'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'RBTorsionsForce':
                    forceDirective.header = ';   ai    aj    ak    al funct       c0          c1          c2          c3          c4          c5'
                    # Cycle through each torsion
                    for i in range(len(force.getNumTorsions())):
                        fields = force.getTorsionParameters(i)
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
                        gromacsForce = '%(particle1)s%(particle2)s%(particle3)s%(particle4)s%(func)d%(C0)f%(C1)f%(C2)f%(C3)f%(C4)f%(C5)f'%vars()
                        forceDirective.lines.append(gromacsForce)

                if forceInteractionType == 'FourierTorsionForce':
                    forceDirective.header = ';   ai    aj    ak    al funct       c0          c1          c2          c3          c4'
                    # Cycle through each torsion
                    for i in range(len(force.getNumTorsions())):
                        fields = force.getTorsionParameters(i)
                        particle1 = fields[0]
                        particle2 = fields[1]
                        particle3 = fields[2]
                        particle4 = fields[3]
                        func = 4
                        C1 = fields[4]._value
                        C2 = fields[5]._value
                        C3 = fields[6]._value
                        C4 = fields[7]._value
                        gromacsForce = '%(particle1)s%(particle2)s%(particle3)s%(particle4)s%(func)d%(C1)f%(C2)f%(C3)f%(C4)f'%vars()
                        forceDirective.lines.append(gromacsForce)

            myDirectiveList.append(forceDirective)

        return myDirectivesList


    def ConvertParametersToDirectives(self, GromacsParm):
        """
        Creates a Directive from a list of GromacsParameters

        """

        myDirectivename = '' # How will we do this part? Perhaps we can pull it from index (which must be filled out whenever a GromacsTopology is created?)
        myDirective = self.GromacsTopologyFileObject.Directive(myDirectiveName, header='') # Does this need to be self.Directive?
        for parm in GromacsParm:
            myDirective.lines.append(parm.directiveString())
        return myDirective


    def writeTopologyFile(self, filename, ExpandIncludes=False, RebuildDirectives=False):
        """Write the Gromacs topology *.top file to file.

        """

        # Rebuild the Directives from the information in the Topology objects.
        if (RebuildDirectives):

            self.GromacsTopologyFileObject = GromacsTopologyFile()

            # Translate all the GromacsParameter objects to parameter directives
            #for parm in self.parameters:
            #    self.GromacsTopologyFileObject.ParameterDirectives.append( self.ConvertGromacsParameterToParameterDirective(parm) )

            # Translate each set of Topology objects to molecule directives
            for mol in range(len(self.molecules)):
                self.GromacsTopologyFileObject.MoleculeDefinitionDirectives.append( self.ConvertTopologyToMoleculeDirectives(self.molecules[mol] ) )

        self.GromacsTopologyFileObject.write(filename, ExpandIncludes=ExpandIncludes)

        return


class GromacsTopologyFile(object):

    def __init__(self, topfile=None, defines=None):
        """A class to store and modify information in Gromacs *.top topology files.

        A Gromacs *.top topology file stores not only the connectivity for each molecule in the system, but (usually as #includes),
        all of the forcefield parameters for all possible (or at least all common) combinations of bonded atom types.

        Some information about default #define parameters, from gmx 3.1.4 manual:

        "You can use any defines to control options in your customized topology files. Options that are
        already available by default are:

        define = -DFLEX_SPC
            Will tell grompp to include FLEX SPC in stead of SPC into your topology, this is necessary to make conjugate
            gradient work and will allow steepest descent to minimize further.
        define = -DPOSRE
            Will tell grompp to include posre.itp into your topology, used for position restraints.
        """

        if (defines):
            self.defines = defines
        else:
            self.defines = list()  # A list of strings that have been #define'd

        # Add default #defines
        self.addDefine('FLEX_SPC')
        self.addDefine('POSRE')

        # *.top file contents will be stored as lines
        self.lines = []

        # lines after preprocessing
        self.processed_lines = []

        # Store an ordered list of Directive objects, to be created after parsing the lines 
        self.directives = []
        self.defineDirectives = list()


        self.ParameterDirectives = list()
        self.MoleculeDefinitionDirectives = list()  # this is a list of lists (for multiple molecules)
        self.SystemDirectives = list()

        # Read in the contents of the topology file, if supplied
        if not topfile == None:
            self.read(topfile)

        """
        self.includes = []      # raw lines of text, i.e. '#include "ffamber03.itp"'
        self.ifdefs = []        # list of list of raw lines of text for each #ifdef (#ifdef'ed  #includes are in here too)
        self.molecules = []     # a list of TopologyFileMolecule object classes
        self.systemTitle = ''   # a one-line title for the system
        self.nmols = {}         # a dictionary of {molecule name: #mols}

        self.nmolsOrder = []    # It's very important to preserve the *ordering* of the defined [ moleculetype ]
                                # at the end of the *.top file.  This list specifies the ordering of the self.nmols.keys()

        if not topfile == None:
            [self.lines, self.includes, self.ifdefs, \
             self.molecules, self.systemTitle, self.nmols, self.nmolsOrder] = self.read(topfile)
        """


    def addDefine(self, defineString):
        """Add a #define string to the list of self.defines."""

        if self.defines.count(defineString) == 0:
            self.defines.append(defineString)
        return 


    def read(self, filename):
        """
        Read in the contents of a Gromacs topology file.

        WHAT ARE WE PARSING?

        The contents of a *.top file is formatted in blocks of Directives.  A Directive looks like this: [ my_directive ]
        From the gmx 3.1.4 manual, here's the description of what to expect:

            * semicolon (;) and newline surround comments
            * on a line ending with \ the newline character is ignored.
            * directives are surrounded by [ and ]
            * the topology consists of three levels:
              - the parameter level (see Table 5.3)
              - the molecule level, which should contain one or more molecule definitions (see Table 5.4)
              - the system level: [ system ], [ molecules ]


        HOW DO WE PARSE THIS?

        The basic strategy is straighforward:  For each Directive's block of text, we'll create an ordered list of
        Directive() objects, which are private container classess.  The only wrinkles here are that Gromacs *.top topology
        file have preprocesser statements:

            #include "myfile.itp"

                Inserts a file's contents at a specific position in the *.top file

            #ifdef USE_THIS
            <statements>
            #endif

                Only reads the statements if the variable USE_THIS is the defined (in the *.mdp, usually)

        To deal with this, we do our *own* preprocessing when we read in the file, in two passes:

            1) We parse all #ifdef statements according to our internal list in self.defines (these can be set by the user before reading the *.top file using methods addDefine() and delDefine()  )
            2) We consider each #include as a special kind of Directive class called an IncludeDirective, which in turn stores it's own GromacsTopology object (!)  Arbitrary levels of include-nesting thus be read, and their parameters unpacked into the container classes ParameterInfo or MoleculeDefinitionInfo classes.  However, since the original structure of nested #includes (wel, at least the top level) is still intact, we can write a *.top to file using only the highest-level (i.e. original) top level.

        There also is a convention for comments that we assume is in place:

            1) Any free-floating "; <comment>" line gets its own Directive (CommentDirective)
            2) Any "; <comment>"s betwewn a directive statement and directive data is considered to be the "header" for that directive

        """

        # Read in the *.top topology file lines
        self.lines = self.readLines(filename)

        # Perform pre-processing of the #ifdef and #endif lines
        self.processedLines = self.preprocess()

        # Parse the contents into a list of Directives
        self.parseDirectives()

        # Organize the Directives into Parameter, MoleculeDefinition, and System groups
        self.organizeDirectives()

    def readLines(self, filename):

        fin = open(filename, 'r')
        lines = fin.readlines()
        fin.close()
        return lines


    def preprocess(self):

        processedLines = []

        # Does this need to be modified?###########################################
        for line in self.lines:
            processedLines.append(line)
        # Concatenate any line continuations
        i = 0
        while i < len(processedLines)-1:
            if len(processedLines[i].strip()) > 0:
                if processedLines[i].strip()[-1] == '\\':
                    processedLines[i] = processedLines[i].strip() + processedLines.pop(i+1)
            i += 1
        return processedLines


    def parseDirectives(self, debug=False):
        """Parse the lines of the topology file contents into an ordered list of Directive containers."""

        self.directives = list()

        index = 0

        if debug:
            print '*** Parsing the following lines: ***'
            for line in self.processedLines:
                print line.strip()

        while index < len(self.processedLines):

            if debug: print '### line',index, ':', self.processedLines[index].strip()

            if len(self.processedLines[index]) > 0:

                # Is the line an include?
                if self.processedLines[index].count('#include') > 0 and self.processedLines[index][0] == '#':
                    # Then create an IncludeDirective
                    if debug: print '### line', index,': creating an IncludeDirective' 
                    self.directives.append( self.IncludeDirective(self.processedLines[index]) )

                    # Go to the next line
                    index += 1

                # Is the line a comment (or a set of comments)?
                elif self.processedLines[index][0] == ';':
                    if debug: print '### line', index,': found comment:', self.processedLines[index].strip()
                    linetxt = ''
                    while self.processedLines[index][0] == ';':
                        if debug: print '### line', index,': found continued comments:', self.processedLines[index].strip()
                        linetxt += self.processedLines[index]
                        index += 1
                        if not (index < len(self.processedLines)):
                            break

                    # Then create an CommentDirective
                    if debug: print '### Adding comment:', linetxt
                    self.directives.append( self.CommentDirective(linetxt) )

                    # Go to the next line?   No - index is already incremented in the while loop.

                # Is the line a define?
                elif self.processedLines[index].count('#define') > 0:
                    # Then create an DefineDirective
                    if debug: print '### line', index,': creating DefineDirective'
                    fields = self.processedLines[index][8:].split(' ', 1)
                    name = fields[0].strip()
                    if len(fields) > 1:
                        line = fields[1]
                    else:
                        line = ''
                    self.directives.append( self.DefineDirective(name, line) )
                    # Go to the next line
                    index += 1

                # Is the line a blank line?
                elif self.processedLines[index].strip() == '':
                    # ignore it
                    if debug: print '### line', index,': ignoring blank line'

                    # Go to the next line
                    index += 1

                # Is the line the start of a directive?
                elif self.processedLines[index][0] == '[' and self.processedLines[index].count(']') > 0:

                    if debug: print '### line', index, ': Found a directive:', self.processedLines[index].strip()
                    myDirectiveName = self.processedLines[index].strip()
                    myDirective = self.Directive(myDirectiveName, header='')
                    index += 1
                    # until we hit a blank line (or the end of the file) ...
                    while len(self.processedLines[index].strip()) > 0:
                        # Test for an ifdef
                        if self.processedLines[index].count('#ifdef') > 0:
                            fields = self.processedLines[index].split()
                            if len(fields) > 1:
                                if self.defines.count(fields[1]) > 0:
                                    self.processedLines.pop(index)
                                while self.processedLines[index].strip() != '#else':
                                    self.processedLines.pop(index)
                                if self.processedLines[index].strip() == '#else':
                                    self.processedLines.pop(index)
                                    if self.processedLines[index][0] == '[':
                                        break
                                    if self.defines.count(fields[1]) > 0: # need to be able to test for all defines, not just the last one
                                        while self.processedLines[index].strip() != '#endif':
                                            self.processedLines.pop(index)
                        if self.processedLines[index].strip() == '#endif':
                            self.processedLines.pop(index)
                            index -= 1

                        if self.processedLines[index][0] == ';':
                            if debug: print '### line', index, ': Found a directive header:', self.processedLines[index].strip()
                            myDirective.header += self.processedLines[index]
                        # ... and data lines get appended to Directive.lines.
                        elif self.processedLines[index].count('#define') > 0:
                            # Then create an DefineDirective
                            if debug: print '### line', index,': creating DefineDirective'
                            fields = self.processedLines[index][8:].split(' ', 1)
                            name = fields[0].strip()
                            if len(fields) > 1:
                                line = fields[1]
                            else:
                                line = ''
                            self.directives.append( self.DefineDirective(name, line) )
                            # Also, add the #define string to our list of defines
                            self.addDefine(name)
                        else:
                            if debug: print '### line', index, ': Found a directive line:', self.processedLines[index].strip()
                            myDirective.lines.append( self.processedLines[index].replace('\t', '    ') )
                        index += 1

                        if not (index < len(self.processedLines)):
                            break

                    self.directives.append( myDirective )

                    # Go to the next line -- no already incremented in the while loop above.


                else:
                    # Go to the next line
                    index += 1


    def organizeDirectives(self):
        """Organize the Directives into Parameter, MoleculeDefinition, and System groups"""

        # Make an expanded list of Directives from the nested IncludeDirectives, skipping CommentDirectives and DefineDirectives
        defineList = self.getAllDefines([])
        expandedDirectives = self.getAllDirectives([])

        # group the System Directives first
        index = 0
        while index < len(expandedDirectives):
            if expandedDirectives[index].name.count( '[ system ]' ) > 0:
                self.SystemDirectives.append( expandedDirectives.pop(index) )
            elif expandedDirectives[index].name.count( '[ molecules ]' ) > 0:
                self.SystemDirectives.append( expandedDirectives.pop(index) )
            else:
                index += 1

        # group the Parameter Directives next
        index = 0
        while index < len(expandedDirectives) and ( expandedDirectives[index].name.count( '[ moleculetype ]' ) == 0 ):
            self.ParameterDirectives.append( expandedDirectives.pop(index) )

        # finally group lists of lists of molecule directives
        index = 0
        while index < len(expandedDirectives):
            if expandedDirectives[index].name.count( '[ moleculetype ]' ) > 0:
                self.MoleculeDefinitionDirectives.append( [] )
            for line in expandedDirectives[index].lines:
                for define in defineList:
                    if define.name in line:
                        i = expandedDirectives[index].lines.index(line)
                        expandedDirectives[index].lines.remove(line)
                        line = line.replace(define.name, define.line)
                        expandedDirectives[index].lines.insert(i, line)
            self.MoleculeDefinitionDirectives[-1].append(expandedDirectives.pop(index))
        return

    def getAllDefines(self, defineList):
        """Get the complete list of defines by traversing all nexted IncludeDirectives.

        """

        for d in self.directives:
            if d.classType == 'DefineDirective':
                defineList.append(d)
            elif d.classType == 'IncludeDirective':
                defineList = d.GromacsTopologyFileObject.getAllDefines(defineList)
        return defineList

    def getAllDirectives(self, appendToList, IgnoreComments=True, IgnoreDefines=True):
        """Get the complete list of directives by traversing all nested IncludeDirectives.

        """

        for d in self.directives:
            if d.classType == 'Directive':
                appendToList.append(d)
            elif d.classType == 'IncludeDirective':
                appendToList = d.GromacsTopologyFileObject.getAllDirectives(appendToList)
        return appendToList

    def testDirectives(self):
        """Print out all the Directives, in order.

        """

        for d in self.directives:
            print '%%% Directive:', d.classType
            print d

    def write(self, filename, ExpandIncludes=False):
        """Write the GromacsTopologyFile to file.

        """

        fout = open(filename, 'w')
        for d in self.directives:
            if d.classType == 'IncludeDirective':
                if ExpandIncludes:
                    fout.write( repr(d) )
                else:
                    fout.write( d.name + '\n' ) # directives need trailing '\n' for padding
            else:
                fout.write( repr(d) )
        fout.close()

        return

    #=============================================================================================
    # Containers for *.top file DIRECTIVES, e.g.: [ defaults ], [ atomtypes ]
    #=============================================================================================

    class Directive(object):
        """A base class for for the Directive container"""

        def __init__(self, name, header=None):

            self.classType = 'Directive'
            self.name = name        # Can be 'name', or be '[ name ]\n' (It will be fixed to be the latter).
            self.header = header    # any "; <comment>\n" text between the name and data (can have multiple newlines).
            self.lines = []
            return

        def fixName(self):
            """if name is 'name', convert it to '[ name ]' """
            self.name = '[ ' + self.name.replace('[','').replace(']','').strip() + ' ]\n'
            return

        def getName(self):
            return self.name

        def setName(self, name):
            self.name = name
            return

        def getHeader(self):
            return self.header

        def setHeader(self, header):
            self.header = header
            return

        def __repr__(self):
            """Returns a string representation that can be printed or written to file."""

            outtxt = self.name + '\n' + self.header
            for line in self.lines:
                outtxt += line
            outtxt += '\n'  # each directive needs a trailing newline for coding
            return outtxt



    #=============================================================================================
    # special DIRECTIVES

    class CommentDirective(Directive):
        """ any line that starts with a semicolon """

        def __init__(self, comment):

            self.classType = 'CommentDirective'
            self.name = comment
            return

        def __repr__(self):
            """Returns a string representation that can be printed or written to file."""
            return self.name+'\n'  # each directive needs a trailing newline for coding

    class DefineDirective(Directive):
        """ any line that starts with '#define'  """

        def __init__(self, name, line):

            self.classType = 'DefineDirective'
            self.name = name
            self.line = line
            return

        def getName(self):
            return self.name

        def getLine(self):
            return self.line

        def __repr__(self):
            """Returns a string representation that can be printed or written to file."""
            return '#define ' + self.name+'\t'+self.line+'\n'  # each directive needs a trailing newline for coding


    class IncludeDirective(Directive):
        """A special Directive container that points to the filename of the #include file,
        and has as an attribute a self.GromacsTopologyFileObject for that #include file. """

        def __init__(self, includeString):
            self.classType = 'IncludeDirective'
            self.name  = includeString   # should be entire '#include\n'
            self.includeFileName = self.name.replace('#include ','').replace('"','').strip()
            if not os.path.exists(self.includeFileName):
                self.includeFileName = os.path.join(os.environ['GMXLIB'], self.includeFileName)
            self.GromacsTopologyFileObject = GromacsTopologyFile(topfile=self.includeFileName)
            return

        def __repr__(self):
            """Returns a string representation that can be printed or written to file."""
            outtxt = ''
            for d in self.GromacsTopologyFileObject.directives:
                outtxt += repr(d)
            outtxt += '\n'  # each directive needs a trailing newline for coding
            return outtxt
