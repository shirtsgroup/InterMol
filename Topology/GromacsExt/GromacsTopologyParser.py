import pdb
import sys
import os
import string
import re
from collections import deque
from cctools.Atom import *
from cctools.Molecule import *
from cctools.Types import *
from cctools.Force import *
from cctools.CaptureObj import *
from cctools.HashMap import *

class GromacsTopologyParser(object):
    """
    A class containing methods required to read in a Gromacs(4.5.1) Topology File   
    """
    def __init__(self, sys = None, defines=None):
        """
        A simple class to handle the parsing of Gromacs topology files
        
        Generates an instance of this object preprocessing topfile and adding
        the specified defines
        
        Args:
            topfile: The name of a topology file to parse
            defines: a set of default defines to add
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

        self.systemMembers = dict()


        if defines:
            self.defines.union(defines)
        self.defines["FLEX_SPC"] = None
        self.defines["POSRE"] = None
   
    def parseTopology(self, topfile, verbose = False):
        lines = self.preprocess(topfile, verbose)
        if lines:
             self.readTopology(lines, verbose)
 
    def readTopology(self, expanded, verbose = False):
        """
        Read in a Gromacs topology file     

        Args:
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
          (?P<system>system)
          |
          (?P<molecules>molecules))
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

                                newAtomType = AtomCR1Type(split[0],                 # atomtype or name
                                        round(split[1]),                            # bondtype
                                        None,                                       # Z
                                        float(split[2]) * units.amu,                # mass
                                        float(split[3]) * units.elementary_charge,  # charge
                                        split[4],                                   # ptype
                                        sigma * units.kilojoules_per_mole * units.nanometers**(6),    # sigma
                                        epsilon * units.kilojoules_per_mole * units.nanometers**(12))   # epsilon

                            elif self.sys.combinationRule == (2 or 3):
                                newAtomType = AtomCR23Type(split[0],                # atomtype or name
                                        round(split[1]),                            # bondtype
                                        None,                                       # Z    
                                        float(split[2]) * units.amu,                # mass
                                        float(split[3]) * units.elementary_charge,  # charge
                                        split[4],                                   # ptype
                                        float(split[5]) * units.nanometers ,        # sigma
                                        float(split[6]) * units.kilojoules_per_mole)# epsilon

                        self.atomtypes.add(newAtomType)

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

                elif match.group('molecules'):
                    if verbose:
                        print "Parsing [ molecules ]..."
                    expanded.pop(i)
                    counter = -1 # zero indexed offset
                    temp = self.sys.molecules.keys()
                    while i < len(expanded) and not (expanded[i].count('[')):
                        split = expanded.pop(i).split()
                        counter += 1                        
                        self.systemMembers[split[0]] = temp[counter]

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
          (?P<exclusions>exclusions))
          [ ]{1} \]
        """, re.VERBOSE)
        #verbose = True
        i = 0
        moleculeName = None 
        while i < len(expanded):
            match = molDirective.match(expanded[i])
            if match:
                #if verbose:
                    #print match.groups()
                if match.group('moleculetype'):
                    if verbose:
                        print "Parsing [ moleculetype ]..."
                    expanded.pop(i)

                    split = expanded[i].split()

                    moleculeName = split[0]
                    self.sys.molecules[self.systemMembers[moleculeName]].setNrexcl(int(split[1]))
 
                    expanded.pop(i)

                elif match.group('atoms'):
                    if verbose:
                        print "Parsing [ atoms ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        
                        for molecule in self.sys.molecules[self.systemMembers[moleculeName]].moleculeSet:
                            for atom in molecule.atoms:
                                if atom.atomtype == split[1]:
                                    tempType = AbstractAtomType(split[1])
                                
                                    atomType = self.atomtypes.get(tempType)
                                    if not atom.Z:
                                        atom.Z = atomType.Z
                                    if atomType.m == 0:
                                        atom.Z = round(float(split[2])/2)
                                    atom.mass = atomType.m
                                    atom.charge = atomType.q
                                    atom.ptype = atomType.ptype
                                    atom.sigma = atomType.V
                                    atom.epsilon = atomType.W
                                    
                                
                        """
                        TODO
                        if not atomtype:
                            print 'cannot find atomtype'
                        else:
                            pass   
                        """
                elif match.group('bonds'):
                    if verbose:
                        print "Parsing [ bonds ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newBondForce = None

                        if int(split[2]) == 1:
                            if len(split) == 5:
                                newBondForce = Bond(int(split[0]),
                                        int(split[1]),
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-2))

                        self.sys.molecules[self.systemMembers[moleculeName]].bondForceSet.add(newBondForce)
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
                                newPairForce = AbstractPair(split[0], split[1])

                        self.sys.forces.add(newPairForce)
                                

                elif match.group('angles'):
                    if verbose:
                        print "Parsing [ angles ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newAngleForce = None

                        if int(split[3]) == 1:
                            if len(split) == 6:
                                newAngleForce = Angle(split[0],
                                        split[1],
                                        split[2],
                                        split[4] * units.degrees,
                                        split[5] * units.kilojoules_per_mole * units.radians**(-2))
                        self.sys.molecules[self.systemMembers[moleculeName]].angleForceSet.add(newAngleForce)
                        self.sys.forces.add(newAngleForce)
                            
                            

                elif match.group('dihedrals'):
                    if verbose:
                        print "Parsing [ dihedrals ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newDihedralForce = None

                        if int(split[4]) == 1:
                            newDihedralForce = ProperDihedral1(split[0],
                                    split[1],
                                    split[2],
                                    split[3],
                                    split[5] * units.degrees,
                                    split[6] * units.kilojoules_per_mole,
                                    split[7])

                        elif int(split[4]) == 3:
                            newDihedralForce = RBDihedral(split[0],
                                    split[1],
                                    split[2],
                                    split[3],
                                    split[5] * units.kilojoules_per_mole,
                                    split[6] * units.kilojoules_per_mole,
                                    split[7] * units.kilojoules_per_mole,
                                    split[8] * units.kilojoules_per_mole,
                                    split[9] * units.kilojoules_per_mole,
                                    split[10] * units.kilojoules_per_mole)
                    
                        elif int(split[4]) == 4:
                            newDihedralForce = ImproperDihedral4(split[0],
                                    split[1],
                                    split[2],
                                    split[3],
                                    split[5] * units.degrees,
                                    split[6] * units.kilojoules_per_mole,
                                    split[7])

                        elif int(split[4]) == 9:
                            newDihedralForce = ProperDihedral9(split[0],
                                    split[1],
                                    split[2],
                                    split[3],
                                    split[5] * units.degrees,
                                    split[6] * units.kilojoules_per_mole,
                                    split[7])
                        
                        self.sys.molecules[self.systemMembers[moleculeName]].dihedralForceSet.add(newDihedralForce)
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
                            newSettlesForce = Settles(split[0],
                                    split[2] * units.nanometers,
                                    split[3] * units.nanometers)

                        self.sys.forces.add(newSettlesForce)


                elif match.group('exclusions'): 
                    if verbose:
                        print "Parsing [ exclusions ]..."
                    expanded.pop(i)
               
                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        newExclusion = Exclusions(expanded.pop(i).split())
                        
                        self.sys.forces.add(newExclusion)   
 
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
                    while True:     # expande the lines
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

        # [ system ]
        lines.append('[ system ]\n')
        lines.append('%s\n'%(self.sys.name))
        lines.append('\n')

        # [ moleculetype]
        for moleculeType in self.sys.molecules.itervalues(): 
            lines.append('[ moleculetype ]\n')
            lines.append('%s%10d\n'%(moleculeType.name, moleculeType.nrexcl))

            # [ atoms ]
            lines.append('[ atoms ]\n')

            for molecule in moleculeType.moleculeSet:
                for atom in molecule.atoms:
                    lines.append('%6d%18s%6d%8s%8s%6d%8.3f%12.5f\n'%(atom.atomNum, atom.atomtype, atom.resNum, atom.resName, atom.atomName, atom.Z, atom.charge._value, atom.mass._value))
                        

            # [ bonds ]
            lines.append('[ bonds ]\n')
            lines.append(';   ai     aj funct  r               k\n')
            for bond in moleculeType.bondForceSet:

                if isinstance(bond, Bond):
                    type = 1
                    lines.append('%6d%7d%4d%18.8e%18.8e\n'%(bond.atom1, bond.atom2, type, bond.length._value, bond.k._value))
                elif isinstance(bond, G96Bond):
                    type = 2
                    lines.append('%6d%7d%4d%5.8f%5.8f\n'%(bond.atom1, bond.atom2, type, bond.length._value, bond.k._value))
                else:
                    print "ERROR (writeTopology): found unsupported bond type"
       
        """ 
        # [ pairs ]
        lines.append('[ pairs ]\n')
        lines.append(';   ai     aj funct'\n)
        for pair in forceList.something.something:

            if type(pair) == (LJ1PairCR1 or LJ1PairCR23):
                type = 1
                atom1, atom2, V, W = pair.getForceParameters()
                lines.append('%6d%7d%4d%5.8f%5.8f\n'%(atom1, atom2, type, V, W))

            elif type(pair) == (LJ2PairCR1 or LJ2PairCR23):
                type = 2
                atom1, atom2, fudgeQQ, qi, qj, V, W = pair.getForceParameters()
                lines.append('%6d%7d%4d%5.8f%5.8f%5.8f%5.8f%5.8f\n'%(atom1, atom2, type, fudgeQQ, qi, qj, V, W))

            elif type(pair) == (LJNBPairCR1 or LJNBPairCR23):
                type = 1
                atom1, atom2, qi, qj, V, W = pair.getForceParameters()
                lines.append('%6d%7d%4d%5.8f%5.8f%5.8f%5.8f\n'%(atom1, atom2, type, qi, qj, V, W))
            else:
                print "ERROR (writeTopology): found unsupported pair type"

        # [ angles ]
        lines.append('[ angles ]\n')
        lines.append('; i j   k   funct'\n)
        for angle in forceList.something.something:

            if type(angle) == Angle:
                type = 1
                atom1, atom2, atom3, theta, k = pair.getForceParameters()
                lines.append('%6d%7d%7d%4d%5.8f%5.8f\n'%(atom1, atom2, atom3, type, theta, k))
            else:
                print "ERROR (writeTopology): found unsupported  angle type"

        # [ dihedrals ]
        lines.append('[ dihedrals ]\n')
        lines.append(';    i      j      k      l   func '\n)
        for dihedral in forceList.something.something:

            if type(dihedral) == ProperDihedral1:
                type = 1
                atom1, atom2, atom3, atom4, phi, k, multiplicity = pair.getForceParameters()
                lines.append('%6d%7d%7d%7d%4d%5.8f%5.8f%4d\n'%(atom1, atom2, atom3,atom4, type, phi, k, multiplicity))

            elif type(dihedral) == ImproperDihedral2:
                type = 2
                atom1, atom2, atom3, atom4, xi, k = pair.getForceParameters()
                lines.append('%6d%7d%7d%7d%4d%5.8f%5.8f\n'%(atom1, atom2, atom3,atom4, type, xi, k))

            elif type(dihedral) == RBDihedral:
                type = 3
                atom1, atom2, atom3, atom4, C0, C1, C2, C3, C4, C5 = pair.getForceParameters()
                lines.append('%6d%7d%7d%7d%4d%5.8f%5.8f%5.8f%5.8f%5.8f%5.8f\n'%(atom1, atom2, atom3,atom4, type, C0, C1, C2, C3, C4, C5))

            if type(dihedral) == ImproperDihedral4:
                type = 4
                atom1, atom2, atom3, atom4, phi, k, multiplicity = pair.getForceParameters()
                lines.append('%6d%7d%7d%7d%4d%5.8f%5.8f%4d\n'%(atom1, atom2, atom3,atom4, type, phi, k, multiplicity))

            if type(dihedral) == ProperDihedral9:
                type = 9
                atom1, atom2, atom3, atom4, phi, k, multiplicity = pair.getForceParameters()
                lines.append('%6d%7d%7d%7d%4d%5.8f%5.8f%4d\n'%(atom1, atom2, atom3,atom4, type, phi, k, multiplicity))

            else:
                print "ERROR (writeTopology): found unsupported  dihedral type"

        try:
            # [ settles ]
            lines.append('[ settles ]\n')
            lines.append('; i  funct   dOH  dHH\n')
            type = 1
            atom1, dOH, dHH = forceList.something.something.getForceParameters()
            lines.append('%6d%4d%5.8f%5.8f\n'%(atom1, type, dOH, dHH))
        except:
            pass

        try:
            # [ exclusions ]
            lines.append('[ exclusions ]\n')
            atom1, atom2, atom3 = forceList.something.something.getForceParameters()
            lines.append('%6d%6d%6d\n'%(atom1, atom2, atom3))
        except:
            pass
        """



        
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
