import pdb
import sys
import os
import string
import re
import copy
from collections import deque
from Topology.Atom import *
from Topology.Molecule import *
from Topology.Types import *
from Topology.Force import *
from Topology.HashMap import *

class GromacsTopologyParser(object):
    """
    A class containing methods required to read in a Gromacs(4.5.3) Topology File   
    """
    def __init__(self, sys = None, defines=None):
        """
        Initializes a GromacsTopologyParse object which serves to read in a Gromacs
        topology into the abstract representation.
    
        Args:
            sys: Global system passed by reference so values can be stored globally
            defines: Sets of default defines to use while parsing.
        """
        self.sys = sys 
        self.includes = set()       # set storing includes
        self.defines = dict()       # list of defines
        self.comments = list()      # list of comments

        self.atomtypes = HashMap()
        self.bondtypes = HashMap()
        self.pairtypes = HashMap()
        self.angletypes = HashMap()
        self.dihedraltypes = HashMap()
        self.constrainttypes = HashMap()

        if defines:
            self.defines.union(defines)
        self.defines["FLEX_SPC"] = None
        self.defines["POSRE"] = None
   
    def parseTopology(self, topfile, verbose = False):
        """
        Parses a Gromacs topology reading into the abstract

        Args:
            topfile: filename of the file to be parsed
            verbose: verbose output
        """
        lines = self.preprocess(topfile, verbose)
        if lines:
             self.readTopology(lines, verbose)
 
    def readTopology(self, expanded, verbose = False):
        """
        Read in a previously preprocessed Gromacs topology file     

        Args:
            expanded: lines from a preprocessed topology file
            verbose: verbose output
        """
        sysDirective = re.compile(r"""
          \[[ ]{1}
          ((?P<defaults>defaults)
          |
          (?P<atomtypes>atomtypes)
          |
          (?P<bondtypes>bondtypes)
          |
          (?P<pairtypes>pairtypes)
          |
          (?P<angletypes>angletypes)
          |
          (?P<dihedraltypes>dihedraltypes)
          |
          (?P<constrainttypes>constrainttypes)
          |
          (?P<nonbond_params>nonbond_params)
          |
          (?P<system>system))
          [ ]{1} \]
        """, re.VERBOSE)
        i = 0       
        while i < len(expanded):
            match = sysDirective.match(expanded[i])
            if match:
                if verbose:
                    print match.groups()
                if match.group('defaults'):
                    if verbose:
                        print "Parsing [ defaults ]..."
                    expanded.pop(i)
                    
                    fields = expanded[i].split()
                    self.sys.nbFunc = int(fields[0])
                    self.sys.combinationRule = int(fields[1])
                    self.sys.genpairs = fields[2]
                    self.sys.ljCorrection = float(fields[3])
                    self.sys.coulombCorrection = float(fields[4])
                    
                    expanded.pop(i)

                elif match.group('atomtypes'):
                    if verbose:
                        print "Parsing [ atomtypes ]..."
                    expanded.pop(i)                    
                
                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newAtomType = None

                        if len(split) == 7:
                            if self.sys.combinationRule == 1:
                                
                                # TODO: double check the following equations
                                sigma = (float(split[6])/float(split[5]))**(1/6)
                                epsilon = float(split[5])/(4*sigma**6)

                                newAtomType = AtomCR1Type(split[0].strip(),         # atomtype or name
                                        split[1].strip(),                           # bondtype
                                        -1,                                         # Z
                                        float(split[2]) * units.amu,                # mass
                                        float(split[3]) * units.elementary_charge,  # charge
                                        split[4],                                   # ptype
                                        sigma * units.kilojoules_per_mole * units.nanometers**(6),      # sigma
                                        epsilon * units.kilojoules_per_mole * units.nanometers**(12))   # epsilon

                            elif self.sys.combinationRule == (2 or 3):
                                newAtomType = AtomCR23Type(split[0].strip(),        # atomtype or name
                                        split[1].strip(),                           # bondtype
                                        -1,                                         # Z    
                                        float(split[2]) * units.amu,                # mass
                                        float(split[3]) * units.elementary_charge,  # charge
                                        split[4],                                   # ptype
                                        float(split[5]) * units.nanometers ,        # sigma
                                        float(split[6]) * units.kilojoules_per_mole)# epsilon

                        self.sys.atomtypes.add(newAtomType)

                elif match.group('bondtypes'):
                    if verbose:
                        print "Parsing [ bondtypes ]..."                    
                    expanded.pop(i)
                    
                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newBondType = None
            
                        # Bond 
                        if int(split[2]) == 1:
                            newBondType = BondType(split[0], 
                                    split[1],
                                    split[2], 
                                    float(split[3]) * units.nanometers, 
                                    float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-2))

                        # G96Bond
                        elif int(split[2]) == 2:
                            newBondType = G96BondType(split[0],
                                        split[1],
                                        split[2],
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-4))

                        # Morse
                        elif int(split[2]) == 3:
                            newBondType = MorseBondType(split[0],
                                        split[1],
                                        split[2],
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole,
                                        float(split[5]) * units.nanometers**(-1))


                        # Cubic 
                        elif int(split[2]) == 4:
                            newBondType = CubicBondType(split[0],
                                        split[1],
                                        split[2],
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-2),
                                        float(split[5]) * units.kilojoules_per_mole * units.nanometers**(-3))

                        # Harmonic 
                        elif int(split[2]) == 6:
                            newBondType = HarmonicBondType(split[0],
                                        split[1],
                                        split[2],
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-2))
                        else:
                            # TODO make a more descriptive error
                            print "could not find bond type"

                        if newBondType:
                            self.bondtypes.add(newBondType)

                elif match.group('pairtypes'):
                    if verbose:
                        print "Parsing [ pairtypes ]..."                   
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        
                        split = expanded.pop(i).split()
                        newPairType = None
                        if int(split[2]) == 1:
                            # LJ/Coul. 1-4 (Type 1)
                            if len(split) == 5:
                                if self.sys.combinationRule == 1:
                                    newPairType = LJ1PairCR1Type(split[0],
                                            split[1],
                                            split[2],
                                            float(split[3]) * units.kilojoules_per_mole * units.nanometers**(6),
                                            float(split[4]) * units.kilojoules_per_mole * units.nanometers**(12))

                                elif self.sys.combinationRule == (2 or 3):
                                    newPairType = LJ1PairCR1Type(split[0],
                                            split[1],
                                            split[2],
                                            float(split[3]) * units.nanometers,
                                            float(split[4]) * units.kilojoules_per_mole)

                            # LJ/C. pair NB
                            elif len(split) == 7:
                                if self.sys.combinationRule == 1:
                                    newPairType = LJ1PairCR1Type(split[0],
                                            split[1],
                                            split[2],
                                            float(split[3]) * units.elementary_charge,
                                            float(split[4]) * units.elementary_charge,
                                            float(split[5]) * units.kilojoules_per_mole * units.nanometers**(6),
                                            float(split[6]) * units.kilojoules_per_mole * units.nanometers**(12))

                                elif self.sys.combinationRule == (2 or 3):
                                    newPairType = LJ1PairCR1Type(split[0],
                                            split[1],
                                            split[2],
                                            float(split[3]) * units.elementary_charge,
                                            float(split[4]) * units.elementary_charge,
                                            float(split[5]) * units.nanometers,
                                            float(split[6]) * units.kilojoules_per_mole)


                        # LJ/Coul. 1-4 (Type 2) 
                        elif int(split[2]) == 2:
                            if self.sys.combinationRule == 1:
                                newPairType = LJ1PairCR1Type(split[0],
                                        split[1],
                                        split[2],
                                        split[3],
                                        float(split[4]) * units.elementary_charge,
                                        float(split[5]) * units.elementary_charge,
                                        float(split[6]) * units.kilojoules_per_mole * units.nanometers**(6),
                                        float(split[7]) * units.kilojoules_per_mole * units.nanometers**(12))

                            elif self.sys.combinationRule == (2 or 3):
                                 newPairType = LJ1PairCR1Type(split[0],
                                        split[1],
                                        split[2],
                                        split[3],
                                        float(split[4]) * units.elementary_charge,
                                        float(split[5]) * units.elementary_charge,
                                        float(split[6]) * units.nanometers,
                                        float(split[7]) * units.kilojoules_per_mole)

                        else:
                            # TODO: need better error handling 
                            print "could not find pair type"

                        if newPairType:
                            self.pairtypes.add(newPairType)    
                        



                elif match.group('angletypes'):
                    if verbose:
                        print "Parsing [ angletypes ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newAngleType = None

                        # Angle
                        if int(split[3]) == 1:
                            newAngleType = AngleType(split[0],
                                    split[1],
                                    split[2],
                                    split[3],
                                    float(split[4]) * units.degrees,
                                    float(split[5]) * units.kilojoules_per_mole * units.radians**(-2))

                        # G96Angle
                        elif int(split[3]) == 2:
                            newAngleType = G96AngleType(split[0],
                                    split[1],
                                    split[2],
                                    split[3],
                                    float(split[4]) * units.degrees,
                                    float(split[5]) * units.kilojoules_per_mole)

                        # Cross bond-bond 
                        elif int(split[3]) == 3:
                            newAngleType = CrossBondBondAngleType(split[0],
                                    split[1],
                                    split[2],
                                    split[3],
                                    float(split[4]) * units.nanometers,
                                    float(split[5]) * units.nanometers,
                                    float(split[6]) * units.kilojoules_per_mole * units.nanometers**(-2))


                        # Cross bond-angle
                        elif int(split[3]) == 4:
                            newAngleType = CrossBondAngleAngleType(split[0],
                                    split[1],
                                    split[2],
                                    split[3],
                                    float(split[4]) * units.nanometers,
                                    float(split[5]) * units.nanometers,
                                    float(split[6]) * units.nanometers,
                                    float(split[7]) * units.kilojoules_per_mole * units.nanometers**(-2))

                        # Urey-Bradley
                        elif int(split[3]) == 5:
                            newAngleType = UreyBradleyAngleType(split[0],
                                    split[1],
                                    split[2],
                                    split[3],
                                    float(split[4]) * units.degrees,
                                    float(split[5]) * units.kilojoules_per_mole,
                                    float(split[6]) * units.nanometers,
                                    float(split[7]) * units.kilojoules_per_mole)

                    
                        # Quartic
                        elif int(split[3]) == 6:
                            newAngleType = QuarticAngleType(split[0],
                                    split[1],
                                    split[2],
                                    split[3],
                                    float(split[4]) * degrees,
                                    float(split[5]) * units.kilojoules_per_mole,
                                    float(split[6]) * units.kilojoules_per_mole * units.radians*(-1),
                                    float(split[7]) * units.kilojoules_per_mole * units.radians*(-2),
                                    float(split[8]) * units.kilojoules_per_mole * units.radians*(-3),
                                    float(split[9]) * units.kilojoules_per_mole * units.radians*(-4))


                        else:
                            print "could not find angle type"

                        if newAngleType:
                            self.angletypes.add(newAngleType)
                                    

                elif match.group('dihedraltypes'):
                    if verbose:
                        print "Parsing [ dihedraltypes ]..."                    
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        #newPairType = PairType(expanded.pop(i).split())
                        #self.pairtypes.add(newPairType)
                        expanded.pop(i)

                elif match.group('constrainttypes'):
                    if verbose:
                        print "Parsing [ constrainttypes ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        #newPairType = PairType(expanded.pop(i).split())
                        #self.pairtypes.add(newPairType)
                        expanded.pop(i)
 
                elif match.group('nonbond_params'):
                    if verbose:
                        print "Parsing [ nonbond_params ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                       expanded.pop(i) 

                elif match.group('system'):
                    if verbose:
                        print "Parsing [ system ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        expanded.pop(i)
                else:
                    i += 1
            else:
                i += 1

        molDirective = re.compile(r"""
          \[[ ]{1}
          ((?P<moleculetype>moleculetype)
          |
          (?P<atoms>atoms)
          |
          (?P<bonds>bonds)
          |
          (?P<pairs>pairs)
          |
          (?P<angles>angles)
          |
          (?P<dihedrals>dihedrals)
          |
          (?P<constraints>constraints)
          |
          (?P<settles>settles)
          |
          (?P<exclusions>exclusions)
          |
          (?P<molecules>molecules))
          [ ]{1} \]
        """, re.VERBOSE)
        i = 0
        moleculeName = None
        currentMolecule = None
        while i < len(expanded):
            match = molDirective.match(expanded[i])
            if match:
                if match.group('moleculetype'):
                    if verbose:
                        print "Parsing [ moleculetype ]..."
                    expanded.pop(i)

                    split = expanded[i].split()

                    moleculeName = split[0]
                    currentMolecule = Molecule(moleculeName)
                    self.sys.addMolecule(currentMolecule)
                    currentMoleculeType = self.sys.molecules[moleculeName]
                    currentMoleculeType.nrexcl = int(split[1])
                    expanded.pop(i)

                elif match.group('atoms'):
                    if verbose:
                        print "Parsing [ atoms ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        atom = Atom(int(split[0]),          # AtomNum
                                int(split[2]),              # resNum
                                split[3].strip(),           # resName
                                split[4].strip())           # atomName
                        atom.atomtype.insert(0,split[1].strip())
                        atom.cgnr = int(split[5])
                        atom.charge.insert(0, float(split[6]) * units.elementary_charge)
                        atom.mass.insert(0, float(split[7]) * units.amu)

                        if len(split) == 11:
                            atom.atomtype.insert(1, split[8].strip())
                            atom.charge.insert(1, float(split[9]) * units.elementary_charge)
                            atom.mass.insert(1, float(split[10]) * units.amu)


                        index = 0
                        for atomType in atom.atomtype:
                            # Searching for a matching atomType to pull values from
                            tempType = AbstractAtomType(atom.atomtype[index])
                            atomType = self.sys.atomtypes.get(tempType)
                            if atomType:
                                atom.Z = atomType.Z
                                if not atom.bondtype:
                                    if atomType.bondtype:
                                        atom.bondtype = atomType.bondtype
                                    else:
                                        sys.stderr.write("Warning: A suspicious parameter was found in atom/atomtypes. Visually inspect before using.\n")
                                if atom.mass[index]._value < 0:
                                    if atomType.mass._value >= 0:
                                        atom.mass.insert(index, atomType.mass)
                                    else:
                                        sys.stderr.write("Warning: A suspicious parameter was found in atom/atomtypes. Visually inspect before using.\n")
                                
                                
                                # Assuming ptype = A
                                #atom.ptype = atomType.ptype

                                atom.sigma.insert(index, atomType.sigma)
                                atom.epsilon.insert(index, atomType.epsilon)
                                
                            else:
                                sys.stderr.write("Warning: A corresponding AtomType was not found. Insert missing values yourself.\n")
                            index +=1
                        
                        currentMolecule.addAtom(atom)

                elif match.group('bonds'):
                    if verbose:
                        print "Parsing [ bonds ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newBondForce = None
                            
                        if len(split) == 3:
                            tempType = AbstractBondType(split[0], split[1], split[2])
                            bondType = self.sys.bondtypes.get(tempType)
                            
                            if isinstance(bondType, BondType):
                                split.append(length)
                                split.append(k)

                            if isinstance(bondType, G96BondType):
                                split.append(length)
                                split.append(k)

                            if isinstance(bondType, CubicBondType):
                                split.append(length)
                                split.append(C2)
                                split.append(C3)

                            if isinstance(bondType, MorseBondType):
                                split.append(length)
                                split.append(D)
                                split.append(beta)

                            if isinstance(bondType, HarmonicBondType):
                                split.append(length)
                                split.append(k)

                        if int(split[2]) == 1: 
                            try:
                                newBondForce = Bond(int(split[0]),
                                        int(split[1]),
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-2))
                            except:
                                newBondForce = Bond(int(split[0]),
                                        int(split[1]),
                                        split[3],
                                        split[4])

                        if int(split[2]) == 2:
                            try:
                                newBondForce = G96Bond(int(split[0]),
                                        int(split[1]),
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-4))
                            except:
                                newBondForce = G96Bond(int(split[0]),
                                        int(split[1]),
                                        split[3],
                                        split[4])

                        if int(split[2]) == 3:
                            try:
                                newBondForce = MorseBond(int(split[0]),
                                        int(split[1]),
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole,
                                        float(split[5]) * units.nanometers**(-1))
                            except:
                                newBondForce = MorseBond(int(split[0]),
                                        int(split[1]),
                                        split[3],
                                        split[4],
                                        split[5])

                        if int(split[2]) == 4:
                            try:
                                newBondForce = CubicBond(int(split[0]),
                                        int(split[1]),
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-2),
                                        float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-3))
                            except:
                                newBondForce = CubicBond(int(split[0]),
                                        int(split[1]),
                                        split[3],
                                        split[4], 
                                        split[5])

                        if int(split[2]) == 5:
                            try:
                                newBondForce = HarmonicBond(int(split[0]),
                                        int(split[1]),
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-2))
                            except:
                                newBondForce = HarmonicBond(int(split[0]),
                                        int(split[1]),
                                        split[3],
                                        split[4])
                            

                        currentMoleculeType.bondForceSet.add(newBondForce)
                        self.sys.forces.add(newBondForce)

                elif match.group('pairs'):
                    if verbose:
                        print "Parsing [ pairs ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newPairForce = None

                        if int(split[2]) == 1:
                            if len(split) == 3:
                                # this probably won't work due to units
                                newPairForce = AbstractPair(int(split[0]), int(split[1]))

                        currentMoleculeType.pairForceSet.add(newPairForce)
                        self.sys.forces.add(newPairForce)
                                    

                elif match.group('angles'):
                    if verbose:
                        print "Parsing [ angles ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newAngleForce = None

                        if len(split) == 4:
                            tempType = AbstractAngleType(split[0], split[1], split[2], split[3])
                            angleType = self.sys.angletypes.get(tempType)
                     
                            if isinstance(angleType, AngleType):
                                split.append(theta)
                                split.append(k)

                            if isinstance(angleType, G96BondType):
                                split.append(theta)
                                split.append(k)

                            if isinstance(angleType, CrossBondBondAngleType):
                                split.append(r1)
                                split.append(r2)
                                split.append(k)

                            if isinstance(angleType, CrossBondAngleAngleType):
                                split.append(r1)
                                split.append(r2)
                                split.append(r3)
                                split.append(k)

                            if isinstance(angleType, UreyBradleyAngleType):
                                split.append(theta)
                                split.append(k)
                                split.append(r)
                                split.append(kUB)

                            if isinstance(angleType, QuarticAngleType):
                                split.append(theta)
                                split.append(C0)
                                split.append(C1)
                                split.append(C2)
                                split.append(C3)
                                split.append(C4)


                        # Angle 
                        if int(split[3]) == 1:
                            try:
                                newAngleForce = Angle(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        float(split[4]) * units.degrees,
                                        float(split[5]) * units.kilojoules_per_mole * units.radians**(-2))
                            except:
                                newAngleForce = Angle(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        split[4],
                                        split[5])
                        
                        # G96Angle
                        elif int(split[3]) == 2:
                            try:
                                newAngleForce = G96Angle(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        float(split[4]) * units.degrees,
                                        float(split[5]) * units.kilojoules_per_mole)
                            except:
                                newAngleForce = G96Angle(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        split[4],
                                        split[5])

                        # Cross Bond-Bond
                        elif int(split[3]) == 3:
                            try:
                                newAngleForce = CrossBondBondAngle(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        float(split[4]) * units.nanometers,
                                        float(split[5]) * units.nanometers,
                                        float(split[6]) * units.kilojoules_per_mole * units.nanometers**(-2))
                            except:
                                newAngleForce = CrossBondBondAngle(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        split[4],
                                        split[5],
                                        split[6])

                        # Cross Bond-Angle
                        elif int(split[3]) == 4:
                            try:
                                newAngleForce = CrossBondAngleAngle(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        float(split[4]) * units.nanometers,
                                        float(split[5]) * units.nanometers,
                                        float(split[6]) * units.nanometers,
                                        float(split[7]) * units.kilojoules_per_mole * units.nanometers**(-2))
                            except:
                                newAngleForce = CrossBondAngleAngle(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        split[4],
                                        split[5],
                                        split[6],
                                        split[7])


                        # Urey-Bradley
                        elif int(split[3]) == 5:
                            try:
                                newAngleForce = UreyBradleyAngle(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        float(split[4]) * units.degrees,
                                        float(split[5]) * units.kilojoules_per_mole,
                                        float(split[6]) * units.nanometers,
                                        float(split[7]) * units.kilojoules_per_mole)
                            except:
                                newAngleForce = UreyBradleyAngle(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        split[4],
                                        split[5],
                                        split[6],
                                        split[7])
                        
                        # Quartic
                        elif int(split[3]) == 6:
                            try:
                                newAngleForce = Angle(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        float(split[4]) * units.degrees,
                                        float(split[5]) * units.kilojoules_per_mole,
                                        float(split[6]) * units.kilojoules_per_mole * units.radians**(-1),
                                        float(split[7]) * units.kilojoules_per_mole * units.radians**(-2),
                                        float(split[8]) * units.kilojoules_per_mole * units.radians**(-3),
                                        float(split[9]) * units.kilojoules_per_mole * units.radians**(-4))
                            except:
                                newAngleForce = Angle(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        split[4],
                                        split[5],
                                        split[6],
                                        split[7],
                                        split[8],
                                        split[9])

                        currentMoleculeType.angleForceSet.add(newAngleForce)
                        self.sys.forces.add(newAngleForce)
                                
                            

                elif match.group('dihedrals'):
                    if verbose:
                        print "Parsing [ dihedrals ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newDihedralForce = None

                        if len(split) == 5:
                            tempType = AbstractDihedralType(split[0], split[1], split[2], split[3], split[4])
                            dihedralType = self.sys.dihedraltypes.get(tempType)

                            if isinstance(dihedralType, ProperDihedral1Type):
                                split.append(phi)
                                split.append(k)
                                split.append(multiplicity)

                            if isinstance(dihedralype, ProperDihedral9Type):
                                split.append(phi)
                                split.append(k)
                                split.append(multicplicity)

                            if isinstance(dihedralType, ImproperDihedral2Type):
                                split.append(xi)
                                split.append(k)

                            if isinstance(dihedralType, RBDihedralType):
                                split.append(C0)
                                split.append(C1)
                                split.append(C2)
                                split.append(C3)
                                split.append(C4)
                                split.append(C5)

                            if isinstance(dihedralType, ImproperDihedral4Type):
                                split.append(phi)
                                split.append(k)
                                split.append(multiplicity)


                        # Proper Dihedral 1
                        if int(split[4]) == 1:
                            try:
                                newDihedralForce = ProperDihedral1(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        int(split[3]),
                                        float(split[5]) * units.degrees,
                                        float(split[6]) * units.kilojoules_per_mole,
                                        int(split[7]))
                            except:
                                newDihedralFroce = ProperDihedral1(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        int(split[3]),
                                        split[5],
                                        split[6],
                                        split[7])

                        # Improper Dihedral 2 
                        elif int(split[4]) == 2:
                            try:
                                newDihedralForce = ImproperDihedral2(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        int(split[3]),
                                        float(split[5]) * units.degrees,
                                        float(split[6]) * units.kilojoules_per_mole * units.radians**(-2))
                            except:
                                newDihedralFroce = ImproperDihedral2(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        int(split[3]),
                                        split[5],
                                        split[6])

                        # RBDihedral
                        elif int(split[4]) == 3:
                            try:
                                newDihedralForce = RBDihedral(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        int(split[3]),
                                        float(split[5]) * units.kilojoules_per_mole,
                                        float(split[6]) * units.kilojoules_per_mole,
                                        float(split[7]) * units.kilojoules_per_mole,
                                        float(split[8]) * units.kilojoules_per_mole,
                                        float(split[9]) * units.kilojoules_per_mole,
                                        float(split[10]) * units.kilojoules_per_mole)
                            except:
                                newDihedralFroce = RBDihedral(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        int(split[3]),
                                        split[5],
                                        split[6],
                                        split[7],
                                        split[8],
                                        split[9],
                                        split[10])

                        # Improper Dihedral 4
                        elif int(split[4]) == 4:
                            try:
                                newDihedralForce = ImproperDihedral4(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        int(split[3]),
                                        float(split[5]) * units.degrees,
                                        float(split[6]) * units.kilojoules_per_mole,
                                        int(split[7]))
                            except:
                                newDihedralFroce = ImproperDihedral4(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        int(split[3]),
                                        split[5],
                                        split[6],
                                        split[7])

                        # Proper Dihedral 9
                        elif int(split[4]) == 9:
                            try:
                                newDihedralForce = ProperDihedral9(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        int(split[3]),
                                        float(split[5]) * units.degrees,
                                        float(split[6]) * units.kilojoules_per_mole,
                                        int(split[7]))
                            except:
                                newDihedralFroce = ProperDihedral9(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        int(split[3]),
                                        split[5],
                                        split[6],
                                        split[7])

                        currentMoleculeType.dihedralForceSet.add(newDihedralForce)
                        self.sys.forces.add(newDihedralForce)

                elif match.group('constraints'):
                    if verbose:
                        print "Parsing [ constraints ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        expanded.pop(i)

                elif match.group('settles'):
                    if verbose:
                        print "Parsing [ settles ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newSettlesForce = None
            
                        if len(split) == 4:
                            newSettlesForce = Settles(int(split[0]),
                                    float(split[2]) * units.nanometers,
                                    float(split[3]) * units.nanometers)

                        currentMoleculeType.settles = newSettlesForce
                        self.sys.forces.add(newSettlesForce)


                elif match.group('exclusions'): 
                    if verbose:
                        print "Parsing [ exclusions ]..."
                    expanded.pop(i)
               
                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        newExclusion = Exclusions(expanded.pop(i).split())
                        
                        currentMoleculeType.exclusions.add(newExclusion)
                        self.sys.forces.add(newExclusion)   

                elif match.group('molecules'):
                    if verbose:
                        print "Parsing [ molecules ]..."
                    expanded.pop(i)
                    while i < len(expanded) and not (expanded[i].count('[')):
                        split = expanded.pop(i).split()
                        tempMolecule = self.sys.molecules[split[0]].moleculeSet[0]
                        max = int(split[1])
                        n = 1
                        while n < max:
                            mol = copy.deepcopy(tempMolecule)
                            self.sys.addMolecule(mol)
                            n += 1
                else:
                    i += 1

            else:
                i += 1

    def preprocess(self, filename, verbose = False):
        """
        Preprocess a topology file
        
        Preprocesses for #include, #define, #ifdef directives and handles them
        
        Args:
            filename: the filename to preprocess
            verbose: verbose output
        """
        expanded = list()
        lineReg = re.compile(r"""
         [ ]*                   # omit spaces
         (?P<directive>[^\\;\n]*)   # search for directive section
         [ ]*                   # omit spaces
         (?P<ignore_newline>\\)?        # look for the \ for newline
         (?P<comment>;.*)?          # look for a comment
        """, re.VERBOSE)
        preReg = re.compile(r"""
         (?:\#include[ ]+)
          (?:"|<)
          (?P<include>.+\.itp|.+\.top)
          (?:"|>)
         |  
         (?:\#define[ ]+)
          (?P<defineLiteral>[\S]+)
          (?P<defineData>[ ]+.+)?
         |
         (?:\#undef[ ]+)
          (?P<undef>[\S]+)
         |
         (?:\#ifdef[ ]+)
          (?P<ifdef>[\S]+)
         |
         (?:\#ifndef[ ]+)
          (?P<ifndef>[\S]+)
         |
         (?P<else>\#else)
         |
         (?P<endif>\#endif)
        """, re.VERBOSE)
        condDepth = 0
        write = [True]  
        fd = self.open(filename)
        if fd != None:
            lines = deque(fd)   # initial read of root file
            fd.close()
            while(lines):       # while the queue isn't empty continue preprocessing
                line = lines.popleft()
                if verbose:
                    print "=========================" 
                    print "Original line:", line
                match = lineReg.match(line) # ensure that the line matches what we expect!
                if match:
                    line = match.group('directive')
                    comment = match.group('comment')
                    if verbose:
                        print "Directive:", line
                        print "Comment:", comment
                    while True:     # expand the lines
                        if match.group('ignore_newline'):
                            if verbose: print "Ignore newline found"
                            if lines:
                                temp = lines.popleft()
                                match = lineReg.match(temp)
                                line += match.group('directive')
                            else:
                                print "WARNING: Previous line continues yet EOF encountered"
                                break   # EOF so we can't loop again        
                        else:
                            break   # line terminates   
                    line = line.strip() # just in case theres whitespace we missed
                        
                    match = preReg.match(line)
                    if match != None:
                        if match.group('ifdef') and write[condDepth]:
                            if verbose:
                                print 'Found an ifdef:', match.group('ifdef')
                            condDepth += 1
                            ifdef = match.group('ifdef')
                            if ifdef in self.defines:
                                write.append(True)
                            else:
                                write.append(False) 
                        elif match.group('ifndef') and write[condDepth]:
                            if verbose:
                                print 'Found an ifndef:', match.group('ifndef')
                            condDepth += 1
                            ifndef = match.group('ifndef')
                            if ifndef not in self.defines:
                                write.append(True)
                            else:
                                write.append(False) 
                        elif match.group('else'):
                            if verbose:
                                print 'Found an else'
                            if condDepth > 0:
                                write[condDepth] = not write[condDepth]
                            else:
                                print "WARNING: Found an else not associated with a conditional"
                        elif match.group('endif'):
                            if verbose:
                                print 'Found an endif'
                            if condDepth > 0:
                                condDepth -= 1
                                write.pop()
                            else:
                                print "ERROR: Found an endif not associated with a conditional"
                        elif write[condDepth]:
                            if match.group('include'):
                                if verbose:
                                    print "Found a include:", line
                                include = match.group('include')
                                fd = self.open(include)
                                if fd != None:
                                    tempList = list(fd)
                                    tempList.reverse()
                                    lines.extendleft(tempList)
                            elif match.group('defineLiteral'):
                                if verbose:
                                    print "Found a define:", line
                                define = match.group('defineLiteral')
                                if define not in self.defines:
                                    self.defines[define] = match.group('defineData')
                                else:
                                    print "WARNING: Overriding define:", define
                                    self.define[define] = match.group('defineData')
                            elif match.group('undef'):
                                if verbose:
                                    print "Found a undefine:", line
                                undef = match.group('undef')
                                if undef in self.defines:
                                    self.defines.pop(undef)     
                    elif write[condDepth]:
                        if line != '':
                            for define in self.defines:
                                if define in line:
                                    line = line.replace(define, self.defines[define])
                            if verbose:
                                print "Writing:", line
                            expanded.append(line)
                        if comment:
                            self.comments.append(comment)
                else:
                    print "ERROR: Unreadable line!" 
            return expanded
        else:
            fd.close() 
    
    def open(self, filename):
        """
        Open a file and add to includes list
    
        Checks if a file exists in the current directory or in the GMXLIB environment then either opens or fails
        
        Args:
            filename: the name of the file to open
    
        Returns:
            file descriptor of the file to be openend
        
        Raises:
            IOError when the file cannot be opened for any reason 
        """
        if filename in self.includes:
            print "WARNING: Omitting file ", filename, ". It has already been included!"
            return None
        temp = filename
        if os.path.exists(filename):
            try:
                fd = open(filename)
                self.includes.add(temp)
                return fd
            except IOError, (errno, strerror):
                sys.stderr.write("I/O error(%d): %s for local instance of '%s'\n" % (errno,strerror,filename))
        filename = os.path.join(os.environ['GMXLIB'], filename)
        try:
            fd = open(filename)
            self.includes.add(temp)
            return fd
        except IOError, (errno, strerror):
            sys.stderr.write("Unrecoverable I/O error(%d): %s: '%s'\n" %(errno, strerror,filename))
            sys.exit()
    
    def writeTopology(self, filename):
        """
        Write this topology to file
    
        Write out this topology in Gromacs format
        
        Args:
            filename: the name of the file to write out to
        """
        lines = list()
       
        # [ defaults ]
        lines.append('[ defaults ]\n')
        lines.append('; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n')
        lines.append('%d%16d%18s%20.4f%8.4f\n'%(self.sys.nbFunc, self.sys.combinationRule,
              self.sys.genpairs, self.sys.ljCorrection, self.sys.coulombCorrection))
        lines.append('\n')


        # [ atomtypes ]
        lines.append('[ atomtypes ]\n')
        lines.append(';type, bondtype, Z, mass, charge, ptype, sigma, epsilon\n')
        for atomtype in self.sys.atomtypes.itervalues():
            lines.append('%-11s%5s%6d%18.8f%18.8f%5s%18.8e%18.8e\n'%(atomtype.atomtype, atomtype.bondtype, atomtype.Z, atomtype.mass._value, atomtype.charge._value, atomtype.ptype, atomtype.sigma._value, atomtype.epsilon._value)) 
        lines.append('\n')
 
        # [ moleculetype]
        for moleculeType in self.sys.molecules.itervalues(): 
            lines.append('[ moleculetype ]\n')
            lines.append('%s%10d\n\n'%(moleculeType.name, moleculeType.nrexcl))
                
            # [ atoms ]
            lines.append('[ atoms ]\n')
            lines.append(';num, type, resnum, resname, atomname, cgnr, q, m\n')
            molecule = moleculeType.moleculeSet[0]
            count = 1
            for atom in molecule.atoms:
                try:
                    lines.append('%6d%18s%6d%8s%8s%6d%18.8f%18.8f%18s%18.8f%18.8f\n'%(count, atom.atomtype[0], atom.resNum, atom.resName, atom.atomName, atom.cgnr, atom.charge[0]._value, atom.mass[0]._value, atom.atomtype[1], atom.charge[1]._value, atom.mass[1]._value))
                except:
                    lines.append('%6d%18s%6d%8s%8s%6d%18.8f%18.8f\n'%(count, atom.atomtype[0], atom.resNum, atom.resName, atom.atomName, atom.cgnr, atom.charge[0]._value, atom.mass[0]._value))

                count +=1
            lines.append('\n')

            if moleculeType.bondForceSet:
                # [ bonds ]
                lines.append('[ bonds ]\n')
                lines.append(';   ai     aj funct  r               k\n')
                for bond in moleculeType.bondForceSet.itervalues():
    
                    if isinstance(bond, Bond):
                        type = 1
                        lines.append('%6d%7d%4d%18.8e%18.8e\n'%(bond.atom1, bond.atom2, type, bond.length._value, bond.k._value))
                    elif isinstance(bond, G96Bond):
                        type = 2
                        lines.append('%6d%7d%4d%5.8f%5.8f\n'%(bond.atom1, bond.atom2, type, bond.length._value, bond.k._value))
                    else:
                        print "ERROR (writeTopology): found unsupported bond type"
                lines.append('\n')       
       
            if moleculeType.pairForceSet:
                # [ pair ]
                lines.append('[ pairs ]\n')
                lines.append(';  ai    aj   funct\n')
                for pair in moleculeType.pairForceSet.itervalues():
  
                    if isinstance(pair, AbstractPair):
                        type = 1
                        lines.append('%6d%7d%4d\n'%(pair.atom1, pair.atom2, type))

                    else:
                        print "ERROR (writeTopology): found unsupported pair type"
                lines.append('\n')

 
            if moleculeType.angleForceSet:
                # [ angles ]
                lines.append('[ angles ]\n')
                lines.append(';   ai     aj     ak    funct   theta         cth\n')
                for angle in moleculeType.angleForceSet.itervalues():
                    if isinstance(angle, Angle):
                        type = 1
                        lines.append('%6d%7d%7d%7d%18.8e%18.8e\n'%(angle.atom1, angle.atom2, angle.atom3, type, angle.theta._value, angle.k._value))
                    else:
                        print "ERROR (writeTopology): found unsupported angle type"                    
                lines.append('\n')
   
            """ 
            # [ pairs]
            lines.append('[ pairs ]\n')
            lines.append(';   ai     aj    funct')
            
            for pair  in moleculeType.pairForceSet.itervalues():
                if isinstance(pair, LJ1PairCR1) or isinstance(pair, LJ1PairCR23)
                    type = 1
                    lines.append('%6d%7d%4d%18.8e%18.8e\n'%(pair.atom1, pair.atom2, type, pair.V._value, pair.W._value)

                elif isinstance(pair, LJ2PairCR1) or isinstance(pair, LJ2PairCR23):
                    type = 2
                    lines.append('%6d%7d%4d%18.8e%18.8e\n'%(pair.atom1, pair.atom2, type, pair.V._value, pair.W._value)) 
                
                elif isinstance(pair, LJNBCR1) or isinstance( pair, LJNBCR23):
                    type = 1
                    lines.append('%6d%7d%4d%18.8f%18.8f%18.8f%18.8f\n'%(pair.atom1, pair.atom2, type, pair.qi._value, pair.qj._value, pair.V._value, pair.W._value))    
                else:
                    print "Could not identify pair!"
            """
  
            if moleculeType.dihedralForceSet:          
                # [ dihedrals ]
                lines.append('[ dihedrals ]\n')
                lines.append(';    i      j      k      l   func\n')
                    
                for dihedral in moleculeType.dihedralForceSet.itervalues():
                    if isinstance(dihedral, ProperDihedral1):
                        type = 1
                        lines.append('%6d%7d%7d%7d%4d%18.8f%18.8f%4d\n'
                                %(dihedral.atom1, 
                                dihedral.atom2, 
                                dihedral.atom3, 
                                dihedral.atom4, 
                                type, 
                                dihedral.phi._value, 
                                dihedral.k._value, 
                                dihedral.multiplicity))
    
                    elif isinstance(dihedral, ImproperDihedral2):
                        type = 2
                        lines.append('%6d%7d%7d%7d%4d%18.8f%18.8f\n'
                                %(dihedral.atom1, 
                                dihedral.atom2, 
                                dihedral.atom3, 
                                dihedral.atom4, 
                                type, 
                                dihedral.xi._value, 
                                dihedral.k._value))
        
                    elif isinstance(dihedral, RBDihedral):
                        type = 3
                        lines.append('%6d%7d%7d%7d%4d%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f\n'
                                %(dihedral.atom1, 
                                dihedral.atom2,
                                dihedral.atom3,     
                                dihedral.atom4, 
                                type, 
                                dihedral.C0._value, 
                                dihedral.C1._value, 
                                dihedral.C2._value, 
                                dihedral.C3._value, 
                                dihedral.C4._value, 
                                dihedral.C5._value))
    
                    elif isinstance(dihedral, ImproperDihedral4):
                        type = 4
                        lines.append('%6d%7d%7d%7d%4d%18.8f%18.8f%4d\n'
                                %(dihedral.atom1, 
                                dihedral.atom2, 
                                dihedral.atom3, 
                                dihedral.atom4, 
                                type, 
                                dihedral.phi._value, 
                                dihedral.k._value, 
                                dihedral.multiplicity))
    
                    elif isinstance(dihedral, ProperDihedral9):
                        type = 9
                        lines.append('%6d%7d%7d%7d%4d%18.8f%18.8f%4d\n'
                                %(dihedral.atom1, 
                                dihedral.atom2, 
                                dihedral.atom3, 
                                dihedral.atom4, 
                                type, 
                                dihedral.phi._value, 
                                dihedral.k._value, 
                                dihedral.multiplicity))
       
                    else:
                        print "ERROR (writeTopology): found unsupported  dihedral type"
                lines.append('\n')
        
        if moleculeType.settles:
            settles = moleculeType.settles
            # [ settles ]
            lines.append('[ settles ]\n')
            lines.append('; i  funct   dOH  dHH\n')
            type = 1 
            lines.append('%6d%6d%18.8f%18.8f\n'%(settles.atom1, type, settles.dOH._value, settles.dHH._value))
            lines.append('\n')

        if moleculeType.exclusions:
            # [ exclusions ]
            lines.append('[ exclusions ]\n')
            for exclusion in moleculeType.exclusions.itervalues():
                lines.append('%6s%6s%6s\n'%(exclusion.exclusions[0], exclusion.exclusions[1], exclusion.exclusions[2]))
            lines.append('\n') 

        # [ system ]
        lines.append('[ system ]\n')
        lines.append('%s\n'%(self.sys.name))
        lines.append('\n')


        # [ molecules ]
        lines.append('[ molecules ]\n')
        lines.append('; Compound        nmols\n')
        for molType in self.sys.molecules:
            lines.append('%-15s%8d\n'%(molType, len(self.sys.molecules[molType].moleculeSet)))

        fout = open(filename, 'w')
        for line in lines:
            fout.write(line)
        fout.close()

    def getComments(self):
        """
        Get comments
        """
        return self.comments
    
    def getExpanded(self):
        """
        Get expanded list of directives
        """
        return self.expanded

    def getIncludes(self):
        """
        Get list of includes
        """
        return self.includes
    
    def getDefines(self):
        """
        Get list of defines
        """
        return self.defines
