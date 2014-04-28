import pdb #FOR DEBUGGING PURPOSES
import shlex
import os
import copy
import string
import re
import numpy as np
from collections import deque

import intermol.unit as units
from intermol.Atom import *
from intermol.Molecule import *
from intermol.System import System
from intermol.Types import *
from intermol.Force import *
from intermol.HashMap import *


#   helper function
def endheadersection(blanksection,header,hlines):
    if blanksection:
        hlines = list()
        hlines.append(header)
        hlines.append('      :::\n')
    else:
        hlines[0] = header
    return hlines

class DesmondParser():
    """
    A class containing methods required to read in Desmond CMS File
    """

    def __init__(self,defines=None):
        """
        Initializes a DesmondParse object which serves to read in a CMS file
        into the abstract representation.

        Args:
        defines: Sets of default defines to use while parsing.
        """
        self.includes = set()       # set storing includes
        self.defines = dict()        # list of defines
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
        self.viparr = 1

        self.fblockpos = []
        self.a_blockpos = []
        self.b_blockpos = []
        self.ffio_blockpos = []

        self.stored_ffio_data = {}  # dictionary of stored ffio_entries

#LOAD FFIO BLOCKS IN FIRST (CONTAINS TOPOLOGY)

    def parse_ffio_block(self,lines,start,end):
            # read in a ffio_block that isn't ffio_ff and split it into the
            # commands and the values.
            # lots of room for additional error checking here, such as whether
            # each entry has the correct number of data values, whether they are the correct type, etc.

            # scroll to the next ffio entry
        while not 'ffio_' in lines[start]:
            # this is not an ffio block! or, we have reached the end of the file
            if ('ffio_' not in lines[start]): 
                start+=1
            if start >= end:
                return 'Done with ffio', 0, 0, 0, start

        lines[start].split()[1]
        components = re.split('\W', lines[start].split()[0]) # get rid of whitespace, split on nonword
        ff_type = components[0]
        ff_number = int(components[1])
        i = start+1
        entry_data = []
        while not ':::' in lines[i]:
            entry_data.append(lines[i].split()[0])
            i+=1
        i+=1 # skip the separator we just found
        entry_values = []
        while not ':::' in lines[i]:
            if lines[i].strip():  # skip the blank spaces.
                entry_values.append(lines[i])
            i+=1
        while '}' not in lines[i]:  # wait until we hit an end to the block
            i+=1
        i+=1 # step past the end of the block

        return ff_type,ff_number,entry_data,entry_values,i
    
    def store_ffio_data(self, ff_type, ff_number, entry_data, entry_values):

        self.stored_ffio_data[ff_type] = {}
        self.stored_ffio_data[ff_type]['ff_type'] = ff_type
        self.stored_ffio_data[ff_type]['ff_number'] = ff_number
        self.stored_ffio_data[ff_type]['entry_data'] = entry_data
        self.stored_ffio_data[ff_type]['entry_values'] = entry_values

    def retrive_ffio_data(self, ff_type):

        return self.stored_ffio_data[ff_type]['ff_type'],self.stored_ffio_data[ff_type]['ff_number'],self.stored_ffio_data[ff_type]['entry_data'],self.stored_ffio_data[ff_type]['entry_values']

    def load_ffio_block(self, lines, moleculeName, start, end, sysDirective, sysDirectiveAtm, verbose = False):

#        Loading in ffio blocks from Desmond format
#        Args:
#            lines: list of all data in CMS format
#            moleculeName: name of current molecule
#            start: beginning of where ffio_ff starts for each molecule
#            end: ending of where ffio_ff ends for each molecule
#           sysDirective: help locate positions of specific data in ffio blocks
#           sysDirectiveAtm: help locate positions of specific data in m_atoms
        i = start
        j = start
        vdwtypeskeys = []
        vdwtypes = []
        sitetypes = []

        split = []
        constraints = []
        temp = []

        currentMolecule = None
        currentMoleculeType = None
        newAtomType = None
        newBondForce = None
        newBondType = None
        newPairType = None
        newPairForce = None
        newAngleType = None
        newAngleForce = None
        newDihedralType = None
        newDihedralForce = None
        stemp = None
        etemp = None
        sigma = None
        epsilon = None

        #There are several sections which require sites information to
        #process.  We keep a flag for this so that we are aware when
        #we have seen the sites
        bPreambleRead = False
        stored_ffio_types = []  # a list of stored ffio_type to keep track 
                              # of the ordering later
        atomlist = OrderedSet()
        namecol = 0
        combrcol = 0
        vdwrcol = 0

        #DEFAULT VALUES WHEN CONVERTING TO GROMACS
        System._sys._nbFunc = 1
        System._sys._genpairs = 'yes'

        if verbose:
            print 'Parsing [ molecule %s]'%(moleculeName)
            print "Parsing [ ffio]"

        while i < end:
            if not bPreambleRead:
                # read the first section for the force field info
                while not (':::' in lines[i]):
                    if 's_ffio_name' in lines[i]:
                        namecol = i-start-1
                    elif 's_ffio_comb_rule' in lines[i]:
                        combrcol = i-start-1
                    elif 's_ffio_vdw_func' in lines[i]:
                        vdwtypercol = i-start-1
                    i+=1
                i+=1 # skip the ':::'    
                # figure out combination rule
                combrule = lines[i+combrcol]
                if re.search("GEOMETRIC", combrule, re.IGNORECASE):
                    if re.search("ARITHMETIC", combrule, re.IGNORECASE):
                        System._sys._combinationRule = 3
                    else:
                        System._sys._combinationRule = 2
                elif re.search("Arithmetic", combrule, re.IGNORECASE):
                    System._sys._combinationRule = 1
                vdwrule = lines[i+vdwtypercol] 
                # MISSING: need to identify vdw rule here -- currently assuming LJ12_6_sig_epsilon!

                # skip to the next ffio entry
                while not re.search('ffio',lines[i]):
                    i+=1
                bPreambleRead = True

            currentMolecule = Molecule(moleculeName)
            ff_type, ff_number, entry_data, entry_values, i = self.parse_ffio_block(lines,i,end)
            stored_ffio_types.append(ff_type)
            self.store_ffio_data(ff_type,ff_number, entry_data, entry_values)

        # Reorder so 'vdwtypes' is first, then 'sites'.  Could eventually get some simplification
        # by putting sites first, but too much rewriting for now.

        stored_ffio_types.insert(0, stored_ffio_types.pop(stored_ffio_types.index('ffio_sites')))
        stored_ffio_types.insert(0, stored_ffio_types.pop(stored_ffio_types.index('ffio_vdwtypes')))
    
        # now process all the data
        for type in stored_ffio_types:
            match = sysDirective.match(type)
            if not match:
                continue
                
            ff_type, ff_number, entry_data, entry_values = self.retrive_ffio_data(type)    

            if match.group('vdwtypes'): 
                # molecule name is at sites, but vdwtypes come
                # before sites. So we store info in vdwtypes and
                # edit it later at sites. Eventually, we should
                # probably move to a model where we store sections
                # we can't use yet, and then process them in the
                # order we want.

                if verbose:
                    print "Parsing [ vdwtypes]..."
                for j in range(ff_number):
                    vdwtypes.append(entry_values[j].split()[3:]) #THIS IS ASSUMING ALL VDWTYPES ARE STORED AS LJ12_6_SIG_EPSILON
                    vdwtypeskeys.append(entry_values[j].split()[1])

            elif match.group('sites'):   #correlate with atomtypes and atoms in GROMACS
                if verbose:                #also edit vdwtypes
                    print "Parsing [ sites]..."
                    
                #set indices to avoid continually calling list functions.    
                ivdwtype = entry_data.index('s_ffio_vdwtype')+1
                icharge = entry_data.index('r_ffio_charge')+1
                imass = entry_data.index('r_ffio_mass')+1
                
                if 'i_ffio_resnr' in entry_data:
                    iresnum = entry_data.index('i_ffio_resnr')+1
                    iresidue = entry_data.index('s_ffio_residue')+1
                    
                cgnr = 0
                for j in range(ff_number):
                    split = entry_values[j].split()
                    if split[1] == "atom":
                        if ('i_ffio_resnr' in entry_data): 
                            atom = Atom(int(split[0]), split[ivdwtype],
                                        int(split[iresnum]), 
                                        split[iresidue])
                        else:
                            # No residuenr, means we will have identical atoms sharing this.
                            atom = Atom(int(split[0]), split[ivdwtype])
                        atom.setAtomType(0, split[ivdwtype])
                        atom.setCharge(0, float(split[icharge])*units.elementary_charge) 
                        atom.setMass(0, float(split[imass]) * units.amu)
                        stemp = float(vdwtypes[vdwtypeskeys.index(split[ivdwtype])][0]) * units.angstroms #was in angstroms
                        etemp = float(vdwtypes[vdwtypeskeys.index(split[ivdwtype])][1]) * units.kilocalorie_per_mole #was in kilocal per mol
                        atom.setSigma(0, stemp)
                        atom.setEpsilon(0, etemp)
                        atom.cgnr = cgnr
                        cgnr+=1

                        currentMolecule.addAtom(atom)
                        if not System._sys._atomtypes.get(AbstractAtomType(atom.getAtomType().get(0))): #if atomtype not in System, add it
                            if System._sys._combinationRule == 1:
                                sigma = (etemp/stemp)**(1/6)
                                epsilon = (stemp)/(4*sigma**6)
                                newAtomType = AtomCR1Type(split[ivdwtypes],             #atomtype/name
                                              split[ivdwtype],                             #bondtype
                                              -1,                               #Z
                                              float(split[imass]) * units.amu,      #mass
                                              float(split[icharge]) * units.elementary_charge,  #charge--NEED TO CONVERT TO ACTUAL UNIT
                                              'A',                             #pcharge...saw this in top--NEED TO CONVERT TO ACTUAL UNITS
                                              sigma * units.kilocalorie_per_mole * angstroms**(6),
                                              epsilon * units.kilocalorie_per_mole * unit.angstro,s**(12))
                            elif (System._sys._combinationRule == 2) or (System._sys._combinationRule == 3):
                                newAtomType = AtomCR23Type(split[ivdwtype], #atomtype/name
                                              split[ivdwtype],                 #bondtype
                                              -1,                   #Z
                                              float(split[imass]) * units.amu,  #mass--NEED TO CONVERT TO ACTUAL UNITS
                                              float(split[icharge]) * units.elementary_charge,  #charge--NEED TO CONVERT TO ACTUAL UNIT
                                              'A',                  #pcharge...saw this in top--NEED TO CONVERT TO ACTUAL UNITS
                                              stemp,
                                              etemp)
                            System._sys._atomtypes.add(newAtomType)

                if len(self.a_blockpos) > 1:  #LOADING M_ATOMS
                    if self.a_blockpos[0] < start:
                        # generate the new molecules for this block; the number of molecules depends on
                        # The number of molecules depends on the number of entries in ffio_sites (ff_number)
                        NewMolecules = self.loadMAtoms(lines, self.a_blockpos[0], i, currentMolecule, ff_number, sysDirectiveAtm, verbose)
                        self.a_blockpos.pop(0)

                # now construct an atomlist with all the atoms
                index = 0
                for molecule in NewMolecules:
                    System._sys.addMolecule(molecule)
                    for atom in molecule._atoms:
                        # does this need to be a deep copy?
                        tmpatom = copy.deepcopy(atom)
                        tmpatom.atomIndex = index
                        atomlist.add(tmpatom)
                        index +=1

                currentMoleculeType = System._sys._molecules[moleculeName]
                currentMoleculeType.nrexcl = 3 #PLACEHOLDER FOR NREXCL...WE NEED TO FIND OUT WHERE IT IS

            elif match.group('bonds'): #may not have all bonds types here yet?
                forces = []
                if len(self.b_blockpos) > 1:  #LOADING M_BONDS
                    #print 'LENGTH OF B_BLOCKPOS: %d'%len(self.b_blockpos)
                    if self.b_blockpos[0] < start:
                        forces = self.loadMBonds(lines,self.b_blockpos[0], i, verbose)
                        currentMoleculeType.bondForceSet = forces[0]
                        System._sys._forces = forces[1]
                        self.b_blockpos.pop(0)
                if verbose:
                    print "Parsing [ bonds]..."

                for j in range(ff_number):
                    split = entry_values[j].split()
                    newBondForce = None
                    # integrate validation of harm_constrained and constraints better.
                    if re.match("HARM_CONSTRAINED", split[3], re.IGNORECASE):
                        try:
                            newBondType = BondType(atomlist[int(split[1])-1].atomName,
                                          atomlist[int(split[2])-1].atomName,
                                          1,
                                          float(split[4]) * units.angstroms, #UNITS IN ANGSTROMS--CHECK
                                          float(split[5]) * units.kilocalorie_per_mole * units.angstroms**(-2),
                                          1)
                        except:
                            newBondType = BondType(atomlist[int(split[1])-1].atomName,
                                          atomlist[int(split[2])-1].atomName,
                                          1,
                                          float(split[4]),
                                          float(split[5]),
                                          1)
                        try:
                            newBondForce = Bond(int(split[1]),
                                           int(split[2]),
                                           float(split[4]) * units.angstroms,
                                           float(split[5]) * units.kilocalorie_per_mole * units.angstroms**(-2),
                                           None,
                                           1)
                        except:
                            newBondForce = Bond(int(split[1]),
                                           int(split[2]),
                                           float(split[4]),
                                           float(split[5]),
                                           None,
                                           1)

                    elif re.match("HARM",split[3], re.IGNORECASE):
                        try:
                            newBondType = BondType(atomlist[int(split[1])-1].atomName,
                                          atomlist[int(split[2])-1].atomName,
                                          1,
                                          float(split[4]) * units.angstroms, #UNITS IN ANGSTROMS
                                          float(split[5]) * units.kilocalorie_per_mole * units.angstroms**(-2),
                                          0)
                        except:
                            newBondType = BondType(atomlist[int(split[1])-1].atomName,
                                          atomlist[int(split[2])-1].atomName,
                                          1,
                                          float(split[4]),
                                          float(split[5]),
                                          0)
                        try:
                            newBondForce = Bond(int(split[1]),
                                           int(split[2]),
                                           float(split[4]) * units.angstroms,
                                           float(split[5]) * units.kilocalorie_per_mole * units.angstroms**(-2),
                                           None,
                                           0)
                        except:
                            newBondForce = Bond(int(split[1]),
                                           int(split[2]),
                                           float(split[4]),
                                           float(split[5]),
                                           None,
                                           0)
                    else:
                        print "ERROR (readFile): found unsupported bond"

                    if newBondForce:
                        if newBondForce in currentMoleculeType.bondForceSet: #bondForceSet already contains i,j, bond order
                            oldBondForce = currentMoleculeType.bondForceSet.get(newBondForce)
                            currentMoleculeType.bondForceSet.remove(newBondForce)
                            newBondForce.order = oldBondForce.order
                        currentMoleculeType.bondForceSet.add(newBondForce)
                        System._sys._forces.add(newBondForce)
                    if newBondType and newBondType not in self.bondtypes:
                        self.bondtypes.add(newBondType)

            elif match.group('pairs'): #GromacsTopology didn't fix this, so fix this when changes are made
                if verbose:
                    print "Parsing [ pairs]..."
                funct1_lj = []
                funct1_cl = []
                funct1_both = []
                ljch=None
                clomb=None
                lj = False
                cl = False
                tempcnt = 0
                for j in range(ff_number):
                    split = entry_values[j].split()
                    if int(split[0]) == 1:
                        if re.match(split[3], "LJ"):
                            ljch = float(split[4])
                            lj = True
                            funct1_lj.append([int(split[1]),int(split[2])])
                        else:
                            clomb = float(split[4])
                            cl = True
                            funct1_cl.append([int(split[1]),int(split[2])])
                    elif lj:
                        if re.match(split[3], "Coulomb"):
                            if not clomb:
                                clomb = float(split[4])
                            cl = True
                            lj = False
                    elif cl:
                        if re.match(split[3], "LJ"):
                            if not ljch:
                                ljch = float(split[4])
                            cl = False
                            lj = True
                    if not int(split[0]) == 1:
                        if lj and [int(split[1]), int(split[2])] in funct1_cl:
                            funct1_cl.remove([int(split[1]), int(split[2])])
                            funct1_both.append([int(split[1]), int(split[2])])
                        elif cl and [int(split[1]), int(split[2])] in funct1_lj:
                            funct1_lj.remove([int(split[1]), int(split[2])])
                            funct1_both.append([int(split[1]), int(split[2])])
                        else:
                            if lj:
                                funct1_lj.append([int(split[1]),int(split[2])])
                            else:
                                funct1_cl.append([int(split[1]),int(split[2])])

                # Logic below can't be right, since this value will be written over.
                if ljch:
                    System._sys._ljCorrection = float(ljch)
                if clomb:
                    System._sys._coulombCorrection = float(clomb)

                for a in funct1_lj:     #PUT ALL PAIRS WITH ONLY LJ INTO LJ1
                    newPairForce = AbstractPair(a[0], a[1], "LJ")
                    currentMoleculeType.pairForceSet.add(newPairForce)
                    System._sys._forces.add(newPairForce)

                for a in funct1_cl: #PUT ALL PAIRS WITH ONLY CL INT LJ1
                    newPairForce = AbstractPair(a[0], a[1], "Coulomb")
                    currentMoleculeType.pairForceSet.add(newPairForce)
                    System._sys._forces.add(newPairForce)

                for a in funct1_both: #PUT ALL PAIRS WITH BOTH LJ AND COULOMB INTO NB
                    newPairForce = AbstractPair(a[0], a[1], "Both")
                    currentMoleculeType.pairForceSet.add(newPairForce)
                    System._sys._forces.add(newPairForce)

            elif match.group('angles'): #add more stuff later once you have more samples to work with
                if verbose:
                    print "Parsing [ angles]..."
                for j in range(ff_number):
                    split = entry_values[j].split()
                    newAngleForce = None
                    # todo: integrate constraints and angle constraint description together better.
                    if re.match("HARM_CONSTRAINED",split[4],re.IGNORECASE):  # this needs to go first because HARM is a substring
                        try:
                            newAngleForce = Angle(int(split[1]),
                                            int(split[2]),
                                            int(split[3]),
                                            float(split[5]) * units.degrees,
                                            float(split[6]) * units.kilocalorie_per_mole * units.radians**(-2),
                                            1)
                        except:
                            newAngleForce = Angle(int(split[1]),
                                            int(split[2]),
                                            int(split[3]),
                                            float(split[5]),
                                            float(split[6]),
                                            1)
                    elif re.match("HARM", split[4],re.IGNORECASE):
                        try:
                            newAngleForce = Angle(int(split[1]),
                                            int(split[2]),
                                            int(split[3]),
                                            float(split[5]) * units.degrees,
                                            float(split[6]) * units.kilocalorie_per_mole * units.radians **(-2))
                        except:
                            newAngleForce = Angle(int(split[1]),
                                            int(split[2]),
                                            int(split[3]),
                                            float(split[5]),
                                            float(split[6]))
                    elif re.match("UB", split[4],re.IGNORECASE):
                        try:
                            newAngleForce = Angle(int(split[1]),
                                            int(split[2]),
                                            int(split[3]),
                                            float(split[5]) * units.degrees,
                                            float(split[6]) * units.kilocalorie_per_mole * units.radians **(-2))
                        except:
                            newAngleForce = Angle(int(split[1]),
                                            int(split[2]),
                                            int(split[3]),
                                            float(split[5]),
                                            float(split[6]))

                    else:
                        print "ERROR (readFile): found unsupported angle in: ",
                        print lines[i]

                    if newAngleForce:
                        currentMoleculeType.angleForceSet.add(newAngleForce)
                        System._sys._forces.add(newAngleForce)

            elif match.group('dihedrals'):
                if verbose:
                    print "Parsing [ dihedrals]..."
                for j in range(ff_number):
                    split = entry_values[j].split()
                    newDihedralForce = None
                    #Proper Diehdral 1 ---NOT SURE ABOUT MULTIPLICITY
                    if re.match(split[5], "PROPER_HARM", re.IGNORECASE):
                        try:
                            newDihedralForce = ProperDihedral1(int(split[1]),
                                               int(split[2]),
                                               int(split[3]),
                                               int(split[4]),
                                               float(split[6]) * units.degrees,
                                               float(split[7]) * units.kilocalorie_per_mole * units.degrees**(-2))
                        except:
                            newDihedralFroce = ProperDihedral1(int(split[1]),
                                               int(split[2]),
                                               int(split[3]),
                                               int(split[4]),
                                               split[6],
                                               split[7],
                                               2)

                    #Improper Diehdral 2 ---NOT SURE ABOUT MULTIPLICITY
                    # These two should be the same function.  Check differences (polymer or protein defn, etc).
                    elif re.match(split[5], "IMPROPER_HARM", re.IGNORECASE):
                        try:
                            newDihedralForce = ImproperDihedral2(int(split[1]),
                                               int(split[2]),
                                               int(split[3]),
                                               int(split[4]),
                                               float(split[6]) * units.radians,
                                               float(split[7]) * units.kilocalorie_per_mole * units.radians**(-2))
                        except:
                            newDihedralForce = ImproperDihedral2(int(split[1]),
                                               int(split[2]),
                                               int(split[3]),
                                               int(split[4]),
                                               split[6],
                                               split[7])

                    #RB Dihedral (Assume for Improper_trig and Proper_trig for now)
                    elif re.match(split[5], "PROPER_TRIG", re.IGNORECASE) or re.match(split[5],"IMPROPER_TRIG", re.IGNORECASE) or re.match(split[5], "DIHEDRAL_TRIG", re.IGNORECASE):
                        # currently, I can't see a difference in Proper_trig and improper_trig angles.
                        try:
                            if float(split[6]) == 0:
                                sign = 1
                            elif float(split[6]) == 180:
                                sign = -1
                            else:
                                print("Currently, can't handle a phase angle that is not 0 or 180 in %s" % (split[5]))
                            c0,c1,c2,c3,c4,c5,c6 = ConvertDihedralFromProperTrigToRB(sign,float(split[7]),float(split[8]),
                                                                                     float(split[9]),float(split[10]),
                                                                                     float(split[11]),float(split[12]),
                                                                                     float(split[13]))
                        except:
                            c0 = 0
                            c1 = 0
                            c2 = 0
                            c3 = 0
                            c4 = 0
                            c5 = 0
                            c6 = 0
                            # do some additional error handling here?
                            print "ERROR (readFile): PROPER_TRIG terms not found"
                        try:    
                            newDihedralForce = RBDihedral(int(split[1]),
                                               int(split[2]),
                                               int(split[3]),
                                               int(split[4]),
                                               c0 * units.kilocalorie_per_mole,
                                               c1 * units.kilocalorie_per_mole,
                                               c2 * units.kilocalorie_per_mole,
                                               c3 * units.kilocalorie_per_mole,
                                               c4 * units.kilocalorie_per_mole,
                                               c5 * units.kilocalorie_per_mole,
                                               c6 * units.kilocalorie_per_mole)
                        except:
                            newDihedralForce = RBDihedral(int(split[1]),
                                               int(split[2]),
                                               int(split[3]),
                                               int(split[4]),
                                               c0,
                                               c1,
                                               c2,
                                               c3,
                                               c4,
                                               c5,
                                               c6)
                    elif (re.match(split[5], "OPLS_PROPER", re.IGNORECASE) or re.match(split[5], "OPLS_IMPROPER", re.IGNORECASE)):
                        try:
                            # Internal desmond formatting for OPLS_PROPER
                            #
                            # r_ffio_c0 = f0 (always 0)  split(6)
                            # r_ffio_c1 = f1             split(7)
                            # r_ffio_c2 = f2             split(8)
                            # r_ffio_c3 = f3             split(9)
                            # r_ffio_c4 = f4             split(10)
                            # r_ffio_c5 = f5 (always 0)  split(11)
                            # r_ffio_c6 = f6 (always 0)  split(12)
                            c0,c1,c2,c3,c4,c5,c6 = ConvertDihedralFromOPLSToRB(float(split[7]),float(split[8]),
                                                                               float(split[9]),float(split[10]))
                        except:
                            f1=0
                            f2=0
                            f3=0
                            f4=0
                            # do some additional error handling here.
                            print "ERROR (readFile): OPLS_PROPER terms not found"
                        try:
                            newDihedralForce = RBDihedral(int(split[1]),
                                               int(split[2]),
                                               int(split[3]),
                                               int(split[4]),
                                               c0 * units.kilocalorie_per_mole,
                                               c1 * units.kilocalorie_per_mole,
                                               c2 * units.kilocalorie_per_mole,
                                               c3 * units.kilocalorie_per_mole,
                                               c4 * units.kilocalorie_per_mole,
                                               c5 * units.kilocalorie_per_mole,
                                               0 * units.kilocalorie_per_mole)
                        except:
                            newDihedralForce = RBDihedral(int(split[1]),
                                               int(split[2]),
                                               int(split[3]),
                                               int(split[4]),
                                               c0 * units.kilocalorie_per_mole,
                                               c1 * units.kilocalorie_per_mole,
                                               c2 * units.kilocalorie_per_mole,
                                               c3 * units.kilocalorie_per_mole,
                                               c4 * units.kilocalorie_per_mole,
                                               c5 * units.kilocalorie_per_mole,
                                               c6 * units.kilocalorie_per_mole)
                    else:
                        print "ERROR (readFile): found unsupported dihedral in:",
                        print line[i]
                    if newDihedralForce:
                        currentMoleculeType.dihedralForceSet.add(newDihedralForce)
                        System._sys._forces.add(newDihedralForce)

                #9 proper dihedrals, funct = 1
                #3 improper dihedrals, funct = 2
                #Ryckaert-Bellemans type dihedrals, funct = 3 and pairs are removed

            elif match.group('constraints'):
                if verbose:
                    print "Parsing [ constraints]..."
                ctype = 1
                funct_pos = 0
                atompos = [] #position of atoms in constraints; spread all over the place
                lenpos = [] #position of atom length; spread all over the place
                tempatom = []
                templength = []
                templen = 0
                for j in range(len(entry_data)):
                    if entry_data[j] == 's_ffio_funct':
                        funct_pos = ctype
                    elif 'i_ffio' in entry_data[j]:
                        atompos.append(ctype)
                    elif 'r_ffio' in entry_data[j]:
                        lenpos.append(ctype)
                    ctype+=1

                for j in range(ff_number):
                    if 'HOH' in entry_values[j] or 'AH' in entry_values[j]:
                        split = entry_values[j].split()
                        tempatom = []
                        templength = []
                        for a in atompos:
                            if not '<>' in split[a]:
                                tempatom.append(int(split[a]))
                            else:
                                tempatom.append(None)
                        for l in lenpos:
                            if not '<>' in split[l]:
                                if 'HOH' in entry_values[j] and l == lenpos[0]:
                                    templength.append(float(split[l])*units.degrees) # Check units?
                                else:
                                    templength.append(float(split[l])*units.angstroms) # Check units?
                            else:
                                templength.append(None*units.angstroms)
                        if 'AH' in split[funct_pos]:
                            templen = int(list(split[funct_pos])[-1])
                        elif 'HOH' in split[funct_pos]:
                            templen = 2    # Different desmond files have different options here.
                        if templen == 1:
                            newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos])
                        elif templen == 2:
                            newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos],tempatom[2],templength[1],None,templength[2])
                        elif templen == 3:
                            newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos],tempatom[2],templength[1],tempatom[3],templength[2])
                        elif templen == 4:
                            newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos],tempatom[2],templength[1],tempatom[3],templength[2],tempatom[4],templength[3])
                        elif templen == 5:
                            newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos],tempatom[2],templength[1],tempatom[3],templength[2],tempatom[4],templength[3],tempatom[5],templength[4])
                        elif templen == 6:
                            newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos],tempatom[2],templength[1],tempatom[3],templength[2],tempatom[4],templength[3],tempatom[5],templength[4],tempatom[6],templength[5])
                        elif templen == 7:
                            newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos],tempatom[2],templength[1],tempatom[3],templength[2],tempatom[4],templength[3],tempatom[5],templength[4],tempatom[6],templength[5],tempatom[7],templength[6])
                        elif templen == 8:
                            newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos],tempatom[2],templength[1],tempatom[3],templength[2],tempatom[4],templength[3],tempatom[5],templength[4],tempatom[6],templength[5],tempatom[7],templength[6],tempatom[8],templength[7])
                    else:
                        print "ERROR (readFile): found unsupported constraint"
                    if newConstraint:
                        currentMoleculeType.constraints.add(newConstraint)
                        System._sys._forces.add(newConstraint)

            elif match.group('exclusions'):
                if verbose:
                    print "Parsing [ exclusions]..."
                for j in range(ff_number):
                    temp = entry_values[j].split()
                    temp.remove(temp[0])
                    newExclusion = Exclusions(temp)
                    currentMoleculeType.exclusions.add(newExclusion)
                    System._sys._forces.add(newExclusion)

            elif match.group('restraints'):
                if verbose:
                    print "Parsing [ restraints]..."

        else:  # no matches
            while '}' not in lines[i]:
                i+=1   # not the most robust if there is nesting in a particular pattern

    def loadMBonds(self, lines, start, end, verbose = False): #adds new bonds for each molecule in System

#        Loading in m_bonds in Desmond format
#        Args:
#            lines: list of all data in CMS format
#           start: beginning of where m_bonds starts for each molecule
#           end: ending of where m_bondsends for each molecule

        if verbose:
            print "Parsing [ m_bonds]..."
        bg = False
        newBondForce = None
        split = []
        i = start
        bondForceSet = HashMap()
        forces = OrderedSet()
        while i < end:
            if ':::' in lines[i]:
                if bg:
                    break
                else:
                    bg = True
                    i+=1
            if bg:
                split = lines[i].split()
                try:
                    newBondForce = Bond(int(split[1]),
                                   int(split[2]),
                                   float(0) * units.angstroms, #BOND ORDERS ARE DIFFERENT, SPLIT3 IS BOND ORDER. NOT ACCURATE CALCULATION
                                   float(0) * units.kilocalorie_per_mole * units.angstroms**(-2),
                                   int(split[3]),
                                   0)
                except:
                    newBondForce = Bond(int(split[1]),
                                   int(split[2]),
                                   float(0),
                                   float(0),
                                   int(split[3]),
                                   0)
                bondForceSet.add(newBondForce)
                forces.add(newBondForce)
            i+=1

        return [bondForceSet, forces]

    def loadMAtoms(self, lines, start, end, currentMolecule, slength, sysDirective, verbose = False): #adds positions and such to atoms in each molecule in System

#        Loading in m_atoms from Desmond format
#        Args:
#            lines: list of all data in CMS format
#           start: beginning of where m_atoms starts for each molecule
#           end: ending of where m_atoms ends for each molecule
#           currentMolecule
#           slength: number of unique atoms in m_atoms, used to calculate repetitions
#           sysDirective: help locate positions of specific data in m_atoms

        if verbose:
            print "Parsing [ m_atom ]..."
        i = start
        xcol = None
        ycol = None
        zcol = None
        rincol = None
        rncol = None
        aicol = None
        an1col = None
        an2col = None
        vxcol = None
        vycol = None
        vzcol = None
        bg = False
        aline1 = ""
        aline2 = ""

        mult = int(re.split('\W',lines[start].split()[0])[1])/slength

        while i < end:
            if ':::' in lines[i]:
                i+=1
                break
            else:
                match = sysDirective.match(lines[i])
                if match:
                    if match.group('First'):
                        start+=1
                    if match.group('xcoord'):
                        if verbose:
                            print "   Parsing [ xcoord]..."
                        xcol = i - start
                    elif match.group('ycoord'):
                        if verbose:
                            print "   Parsing [ ycoord]..."
                        ycol = i - start
                    elif match.group('zcoord'):
                        if verbose:
                            print "   Parsing [ zcoord]..."
                        zcol = i - start
                    elif match.group('rindex'):
                        if verbose:
                            print "   Parsing [ rindex]..."
                        rincol = i - start
                    elif match.group('rname'):
                        if verbose:
                            print "   Parsing [ rname]..."
                        rncol = i - start
                    elif match.group('aindex'):
                        if verbose:
                            print "   Parsing [ aindex]..."
                        aicol = i - start
                    elif match.group('aname1'):
                        if verbose:
                            print "   Parsing [ aname1]..."
                        an1col = i - start
                    elif match.group('aname2'):
                        if verbose:
                            print "   Parsing [ aname2]..."
                        an2col = i - start
                    elif match.group('xvelocity'):
                        if verbose:
                            print "   Parsing [ xvelocity]..."
                        vxcol = i - start
                    elif match.group('yvelocity'):
                        if verbose:
                            print "   Parsing [ yvelocity]..."
                        vycol = i - start
                    elif match.group('zvelocity'):
                        if verbose:
                            print "   Parsing [ zvelocity]..."
                        vzcol = i - start
            i+=1

        atom = None

        newMoleculeAtoms = []
        j = 0
        if verbose:
            print "   Parsing atoms..."

        molecules = []    
        while j < mult:
            newMolecule = copy.deepcopy(currentMolecule)
            for atom in newMolecule._atoms:
                if ':::' in lines[i]:
                    break
                else:
                    aline = shlex.split(lines[i])
                    atom.residueIndex = int(aline[rincol])
                    atom.residueName = aline[rncol].strip()
                    atom.atomIndex = int(aline[aicol])
                    atom.setPosition(float(aline[xcol]) * units.angstroms,
                                     float(aline[ycol]) * units.angstroms,
                                     float(aline[zcol]) * units.angstroms)
                    if vxcol == vycol == vzcol == None:
                        atom.setVelocity(0.0 * units.angstroms * units.picoseconds**(-1),
                                        (0.0 * units.angstroms) * units.picoseconds**(-1),
                                        (0.0 * units.angstroms) * units.picoseconds**(-1))
                    else:
                        atom.setVelocity(float(aline[vxcol]) * units.angstroms * units.picoseconds**(-1),
                                        float(aline[vycol]) * units.angstroms * units.picoseconds**(-1),
                                         float(aline[vzcol]) * units.angstroms * units.picoseconds**(-1))
                    aline1 = aline[an1col].strip()
                    aline2 = aline[an2col].strip()
                    if re.match('$^',aline1) and not re.match('$^',aline2):
                        atom.atomName = aline2
                    elif re.match('$^',aline2) and not re.match('$^',aline1):
                        atom.atomName = aline1
                    elif re.search("\d+",aline1) and not re.search("\d+",aline2):
                        if re.search("\D+",aline1) and re.search("\w+",aline1):
                            atom.atomName = aline1
                        else:
                            atom.atomName = aline2
                    elif re.search("\d+",aline2) and not re.search("\d+",aline1):
                        if re.search("\D+",aline2) and re.search("\w+",aline2):
                            atom.atomName = aline2
                        else:
                            atom.atomName = aline1
                    elif re.match('$^',aline1) and re.match('$^',aline2):
                        atom.atomName = "None"
                    else:
                        atom.atomName = aline2  #doesn't matter which we choose, so we'll go with atom name instead of pdb
                    i+=1

            molecules.append(newMolecule)        
            j+=1

        return molecules


    def loadBoxVector(self, lines, start, end, verbose = False):

#       Loading Box Vector
#       Create a Box Vector to load into the System
#        Args:
#            lines: all the lines of the file stored in an array
#            start: starting position
#            end: ending position

        i = start
        v = np.zeros([3,3])*units.angstroms
        while (i<end):
            if 'r_chorus_box_ax' in lines[i]:
                startboxlabel = i-start
            if ':::' in lines[i]:
                endlabel = i
                break
            i+=1
        startbox = startboxlabel+endlabel
        nvec = 0
        for i in range(startbox,startbox+9):
            j = (nvec)/3
            k = (nvec)%3
            v[j,k] = float(re.sub(r'\s', '', lines[i])) * units.angstrom
            nvec += 1

        System._sys.setBoxVector(v)

    def readFile(self, filename):

#        Load in data from file

#       Read data in Desmond format

#        Args:
#            filename: the name of the file to write out to

        lines = list()

        fl = open(filename, 'r')
        lines = list(fl)
        fl.close()
        i,j=0,0

        for line in lines:
            if re.search("f_m_ct",line,re.VERBOSE):
                if j > 0:
                    self.fblockpos.append(i)
                j+=1
            if re.search("m_atom",line,re.VERBOSE) and not (re.search("i_m",line)
            or re.search("s_m",line)):
                if j > 1:
                    self.a_blockpos.append(i)
                j+=1
            if re.search("m_bond",line,re.VERBOSE):
                if j > 2:
                    self.b_blockpos.append(i)
                j+=1
            if re.search("ffio_ff",line,re.VERBOSE):
                if j > 2:
                    self.ffio_blockpos.append(i)
                j+=1
            i+=1
        i-=1
        self.fblockpos.append(i)
        self.a_blockpos.append(i)
        self.b_blockpos.append(i)
        self.ffio_blockpos.append(i)
        verbose = True

        sysDirectiveTop = re.compile(r"""
          ((?P<vdwtypes>\s*ffio_vdwtypes)
          |
          (?P<sites>\s*ffio_sites)
          |
          (?P<bonds>\s*ffio_bonds)
          |
          (?P<pairs>\s*ffio_pairs)
          |
          (?P<angles>\s*ffio_angles)
          |
          (?P<dihedrals>\s*ffio_dihedrals)
          |
          (?P<constraints>\s*ffio_constraints)
          |
          (?P<exclusions>\s*ffio_exclusions)
          |
          (?P<restraints>\s*ffio_restraints))
        """, re.VERBOSE)

        sysDirectiveStr = re.compile(r"""
          ((?P<xcoord>\s*r_m_x_coord)
          |
          (?P<ycoord>\s*r_m_y_coord)
          |
          (?P<zcoord>\s*r_m_z_coord)
          |
          (?P<rindex>\s*i_m_residue_number)
          |
          (?P<rname>\s*s_m_pdb_residue_name)
          |
          (?P<aindex>\s*i_m_atomic_number)
          |
          (?P<aname1>\s*s_m_pdb_atom_name)
          |
          (?P<aname2>\s*s_m_atom_name)
          |
          (?P<xvelocity>\s*r_ffio_x_vel)
          |
          (?P<yvelocity>\s*r_ffio_y_vel)
          |
          (?P<zvelocity>\s*r_ffio_z_vel)
          |
          (?P<First>\s*[#][\s+]First[\s+]column[\s+]is[\s+]atom[\s+]index[\s+][#]))
        """, re.VERBOSE)

        #LOADING Ffio blocks

        print "Reading Ffio Block..."
        #MRS: warning -- currently no check to avoid duplicated molecule names. Investigate.
        i = 0
        j = 0
        while i < (len(self.ffio_blockpos)-1):
            j = self.fblockpos[i]
            while not re.match(r'\s*[:::]',lines[j]):
                j+=1
            self.load_ffio_block(lines, lines[j+1].strip(), self.ffio_blockpos[i], self.fblockpos[i+1]-1, sysDirectiveTop, sysDirectiveStr,  verbose)
            i+=1
        i = 0

        #LOAD RAW BOX VECTOR-Same throughout cms

        print "Reading Box Vector..."
        self.loadBoxVector(lines, self.fblockpos[0], self.a_blockpos[0], verbose)

# 
    def writeFile(self, filename, verbose=True):

#        Write this topology to file
#        Write out this topology in Desmond format
#        Args:
#            filename: the name of the file to write out to

        lines = list()
        vdwtypes = []
        sites = []
        pos = 0
        name = ''
        sig = None
        ep = None
        stemp = None
        etemp = None

        print "WARNING: MacroModel atom type is set to 7 for all cases."

        # for all CMS files
        lines.append('{\n')
        lines.append('  s_m_m2io_version\n')
        lines.append('  :::\n')
        lines.append('  2.0.0\n')
        lines.append('}\n')

        #FIRST F_M_CT BLOCK

        if verbose:
            print "Writing first f_m_ct..."
        lines.append('f_m_ct {\n')
        lines.append('  s_m_title\n')
        lines.append('  r_chorus_box_ax\n')
        lines.append('  r_chorus_box_ay\n')
        lines.append('  r_chorus_box_az\n')
        lines.append('  r_chorus_box_bx\n')
        lines.append('  r_chorus_box_by\n')
        lines.append('  r_chorus_box_bz\n')
        lines.append('  r_chorus_box_cx\n')
        lines.append('  r_chorus_box_cy\n')
        lines.append('  r_chorus_box_cz\n')
        lines.append('  s_ffio_ct_type\n')
        lines.append('  :::\n')

        #box vector
        bv = System._sys.getBoxVector()
        lines.append('  "full system"\n')
        for bi in range(3):
            for bj in range(3):
                lines.append('%22s\n'%float(bv[bi][bj].in_units_of(units.angstroms)._value))
        lines.append('  full_system\n')

        #M_ATOM
        apos = len(lines) #pos of where m_atom will be; will need to overwite later based on the number of atoms
        lines.append('m_atom\n')
        lines.append('    # First column is atom index #\n')
        lines.append('    i_m_mmod_type\n')
        lines.append('    r_m_x_coord\n')
        lines.append('    r_m_y_coord\n')
        lines.append('    r_m_z_coord\n')
        lines.append('    i_m_residue_number\n')
        lines.append('    s_m_pdb_residue_name\n')
        lines.append('    i_m_atomic_number\n')
        lines.append('    s_m_atom_name\n')
        lines.append('    r_ffio_x_vel\n')
        lines.append('    r_ffio_y_vel\n')
        lines.append('    r_ffio_z_vel\n')
        lines.append('    :::\n')

        i = 0
        nmol = 0
        totalatoms = []
        totalatoms.append(0)
        for moleculetype in System._sys._molecules.itervalues():
            for molecule in moleculetype.moleculeSet:
                for atom in molecule._atoms:
                    i += 1
                    lines.append('    %d        %d   %10.8f %10.8f %10.8f     %2d %4s    %2d  %2s    %11.8f %11.8f %11.8f\n'
                                %(i,
                                7, #NOT SURE WHAT TO PUT FOR MMOD TYPE
                                float(atom._position[0].in_units_of(units.angstroms)._value),
                                float(atom._position[1].in_units_of(units.angstroms)._value),
                                float(atom._position[2].in_units_of(units.angstroms)._value),
                                atom.residueIndex,
                                '"%s"'%atom.residueName,
                                atom.atomIndex,
                                '"%s"'%atom.atomName,
                                float(atom._velocity[0].in_units_of(units.angstroms/units.picoseconds)._value),
                                float(atom._velocity[1].in_units_of(units.angstroms/units.picoseconds)._value),
                                float(atom._velocity[2].in_units_of(units.angstroms/units.picoseconds)._value)))
            totalatoms.append(i)

        lines[apos] = '  m_atom[%d] {\n'%(i)
        lines.append('    :::\n')
        lines.append('  }\n')

        bpos = len(lines)
        i = 0

        #M_BOND
        hlines = list()
        dlines = list()
        hlines.append('  m_bond_placeholder\n')
        hlines.append('    i_m_from\n')
        hlines.append('    i_m_to\n')
        hlines.append('    i_m_order\n')
        hlines.append('    i_m_from_rep\n')
        hlines.append('    i_m_to_rep\n')
        hlines.append('    :::\n')

        i = 0
        nonecnt = 0
        nmol = 0
        for moleculetype in System._sys._molecules.itervalues():
            # sort the bondlist because Desmond requires the first time a bond is listed to have
            # the atoms in ascending order
            bondlist = sorted(moleculetype.bondForceSet.itervalues(), key=lambda x: x.atom1)
            for bond in bondlist:
                if bond and bond.order:
                    i += 1
                    dlines.append('    %d %d %d %d %d %d\n'
                                  %(i,
                                  bond.atom1 + totalatoms[nmol],
                                  bond.atom2 + totalatoms[nmol],
                                  int(bond.order),
                                  1,
                                  1))
                elif not bond:
                    nonecnt+=1
            nmol +=1
        if nonecnt > 0 and verbose:
            print 'FOUND %d BONDS THAT DO NOT EXIST'%nonecnt
        hlines[0] = '  m_bond[%d] {\n' % i    
        if (i > 0):
            lines.extend(hlines)
            lines.extend(dlines)
            lines.append('    :::\n')
            lines.append('  }\n')
            lines.append('}\n')

        solute = True
        endline = ''
        resName = ''

        #WRITE OUT ALL FFIO AND F_M_CT BLOCKS

        for moleculetype in System._sys._molecules.itervalues():
            if verbose:
                print 'Writing molecule block %s...'% moleculetype.name
            #BEGINNING BLOCK
            
            if verbose:
                print "  Writing f_m_ct..."
            lines.append('f_m_ct {\n')
            lines.append('  s_m_title\n')
            bpos = len(lines) #bpos temporarily used for position of s_m_entry_name (for TIP3)
            lines.append('  s_m_entry_name\n')
            lines.append('  i_ffio_num_component\n')
            lines.append('  r_chorus_box_ax\n')
            lines.append('  r_chorus_box_ay\n')
            lines.append('  r_chorus_box_az\n')
            lines.append('  r_chorus_box_bx\n')
            lines.append('  r_chorus_box_by\n')
            lines.append('  r_chorus_box_bz\n')
            lines.append('  r_chorus_box_cx\n')
            lines.append('  r_chorus_box_cy\n')
            lines.append('  r_chorus_box_cz\n')
            lines.append('  s_ffio_ct_type\n')
            lines.append('  :::\n')

            if solute:
                lines.append('  solute\n')
                endline = '  solute\n'
                solute = False
                del lines[bpos]
                del lines[bpos]
            else:
                for atom in molecule._atoms:
                    resName = atom.residueName
                    break
                if re.match("T3P", resName) or re.search("WAT", resName):
                    #lines[bpos] = ('  s_m_entry_name\n')
                    lines.append('  "TIP3P water box"\n')
                    lines.append('  "TIP3P water box"\n')
                    lines.append('  1\n')
                    endline = '  solvent\n'
                else:
                    lines.append('  %s\n'%(moleculetype.name))
                    endline = '  ion\n'
                    del lines[bpos]
                    del lines[bpos] #deletes line for num component (only in TIP3)

            for bi in range(3):
                for bj in range(3):
                    lines.append('%22s\n'%float(bv[bi][bj].in_units_of(units.angstroms)._value))
            lines.append(endline)

            #M_ATOMS

            if verbose:
                print "  Writing m_atoms..."
            apos = len(lines) #pos of where m_atom will be; will need to overwite later based on the number of atoms
            lines.append('m_atom\n')
            lines.append('    # First column is atom index #\n')
            lines.append('    i_m_mmod_type\n')
            lines.append('    r_m_x_coord\n')
            lines.append('    r_m_y_coord\n')
            lines.append('    r_m_z_coord\n')
            lines.append('    i_m_residue_number\n')
            lines.append('    s_m_pdb_residue_name\n')
            lines.append('    i_m_atomic_number\n')
            lines.append('    s_m_atom_name\n')
            lines.append('    r_ffio_x_vel\n')
            lines.append('    r_ffio_y_vel\n')
            lines.append('    r_ffio_z_vel\n')
            lines.append('    :::\n')

            i = 0
            for molecule in moleculetype.moleculeSet:
                for atom in molecule._atoms:
                    i += 1
                    #NOT SURE WHAT TO PUT FOR MMOD TYPE; 7 is currently used.
                    #This can't be determined currently from the information provided,
                    # unless it is stored previous, nor is it used by desmond
                    lines.append('    %d        %d   %10.8f %10.8f %10.8f     %2d %4s    %2d  %2s   %11.8f %11.8f %11.8f\n'
                                %(i,
                                7,
                                float(atom._position[0].in_units_of(units.angstroms)._value),
                                float(atom._position[1].in_units_of(units.angstroms)._value),
                                float(atom._position[2].in_units_of(units.angstroms)._value),
                                atom.residueIndex,
                                '"%s"'%atom.residueName,
                                atom.atomIndex,
                                '"%s"'%atom.atomName,
                                float(atom._velocity[0].in_units_of(units.angstroms/units.picoseconds)._value),
                                float(atom._velocity[1].in_units_of(units.angstroms/units.picoseconds)._value),
                                float(atom._velocity[2].in_units_of(units.angstroms/units.picoseconds)._value)))
            lines[apos] = '  m_atom[%d] {\n'%(i)
            lines.append('    :::\n')
            lines.append('  }\n')

            #M_BONDS
            if verbose:
                print "  Writing m_bonds..."

            hlines = list()
            dlines =  list()
            hlines.append('m_bond_placeholder\n')
            hlines.append('    i_m_from\n')
            hlines.append('    i_m_to\n')
            hlines.append('    i_m_order\n')
            hlines.append('    i_m_from_rep\n')
            hlines.append('    i_m_to_rep\n')
            hlines.append('    :::\n')

            i = 0
            nonecnt = 0
            bondlist = sorted(moleculetype.bondForceSet.itervalues(), key=lambda x: x.atom1)
            for bond in bondlist:
                if bond and bond.order:
                    i += 1
                    dlines.append('    %d %d %d %d %d %d\n'
                                 %(i,
                                   bond.atom1,
                                   bond.atom2,
                                   int(bond.order),
                                   1,
                                   1))
                else:
                    nonecnt+=1
            if nonecnt > 0 and verbose:
                print 'FOUND %d BONDS THAT DO NOT EXIST'%nonecnt
            header = '  m_bond[%d] {\n'%i

            if (i>0):
                hlines = endheadersection(False,header,hlines)
                lines.extend(hlines)
                lines.extend(dlines)
                lines.append('    :::\n')
                lines.append('  }\n')

            #FFIO
            molecule =  moleculetype.moleculeSet[0]
            if verbose:
                print "  Writing ffio..."
            lines.append('  ffio_ff {\n')
            lines.append('    s_ffio_name\n')
            lines.append('    s_ffio_comb_rule\n')
            lines.append('    i_ffio_version\n')
            lines.append('    :::\n')

            #Adding Molecule Name
            if re.search("Viparr", moleculetype.name):
                lines.append('    Generated by Viparr\n')
            else:
                #print moleculetype.name
                lines.append('    %s\n' % moleculetype.name)

            #Adding Combination Rule
            if System._sys._combinationRule == 1:
                lines.append('    ARITHMETIC\n') #NOT SURE WHAT TO PUT HERE...COME BACK TO THIS
            elif System._sys._combinationRule == 2:
                lines.append('    GEOMETRIC\n')
            elif System._sys._combinationRule == 3:
                lines.append('    ARITHMETIC/GEOMETRIC\n')

            #Adding Version
            lines.append('    1.0.0\n') #All files had this, check if version is 1.0.0

            #-ADDING VDWTYPES AND SITES
            i = 1
            vdwtypes = []
            sites = []
            sig = None
            ep = None
            stemp = None
            etemp = None
            combRule = System._sys._combinationRule
            for atom in molecule._atoms:
                if atom.residueIndex:
                    sites.append(' %3d %5s %9.8f %9.8f %2s %1d %4s\n' % (i,'atom',float(atom._charge[0].in_units_of(units.elementary_charge)._value),float(atom._mass[0].in_units_of(units.atomic_mass_unit)._value),atom._atomtype[0],atom.residueIndex,atom.residueName))
                else:
                    sites.append(' %3d %5s %9.8f %9.8f %2s\n' % (i,'atom',float(atom._charge[0].in_units_of(units.elementary_charge)._value),float(atom._mass[0].in_units_of(units.atomic_mass_unit)._value),atom._atomtype[0]))
                sig = float(atom._sigma[0].in_units_of(units.angstroms)._value)
                ep = float(atom._epsilon[0].in_units_of(units.kilocalorie_per_mole)._value)
                if combRule == 1:   #MRS: seems like this should be automated more?
                    stemp = ep * (4 * (sig**6))
                    etemp = stemp * (sig**6)
                elif combRule == 2 or combRule == 3:
                    stemp = sig
                    etemp = ep
                if ' %2s %18s %8.8f %8.8f\n' % (atom._atomtype[0],"LJ12_6_sig_epsilon",float(stemp),float(etemp)) not in vdwtypes:
                    vdwtypes.append(' %2s %18s %8.8f %8.8f\n' % (atom._atomtype[0],"LJ12_6_sig_epsilon",float(stemp),float(etemp)))
                i+=1

            if verbose:
                print "   -Writing vdwtypes..."
            lines.append("    ffio_vdwtypes[%d] {\n"%(len(vdwtypes)))
            lines.append("      s_ffio_name\n")
            lines.append("      s_ffio_funct\n")
            lines.append("      r_ffio_c1\n")
            lines.append("      r_ffio_c2\n")
            lines.append("      :::\n")
            i = 1
            for v in vdwtypes:
                lines.append('      %d%2s'%(i,v))
                i+=1
            lines.append("      :::\n")
            lines.append("    }\n")

            if verbose:
                print "   -Writing sites..."
            lines.append("    ffio_sites[%d] {\n"%(len(sites)))
            lines.append("      s_ffio_type\n")
            lines.append("      r_ffio_charge\n")
            lines.append("      r_ffio_mass\n")
            lines.append("      s_ffio_vdwtype\n")
            if len(sites[0].split()) > 5:
                lines.append("      i_ffio_resnr\n")
                lines.append("      s_ffio_residue\n")
            lines.append("      :::\n")
            for s in sites:
                lines.append('   %s'%(s))
            lines.append("      :::\n")
            lines.append("    }\n")

            #-ADDING BONDS
            if verbose:
                print "   -Writing bonds..."

            dlines = list()
            hlines = list()
            hlines.append('ffio_bonds_placeholder\n')
            hlines.append("      i_ffio_ai\n")
            hlines.append("      i_ffio_aj\n")
            hlines.append("      s_ffio_funct\n")
            hlines.append("      r_ffio_c1\n")
            hlines.append("      r_ffio_c2\n")
            hlines.append("      :::\n")

            nonecnt = 0
            name = ''
            length = None
            k = None

            i = 0
            for bond in moleculetype.bondForceSet.itervalues():
                try:
                    length = float(bond.length.in_units_of(units.angstroms)._value)   #Look at unit conversions here
                    k = float(bond.k._value)  # look at unit conversions here
                except:
                    length = None
                    k = None
                if bond and (length and not length == float(0)) and (k and not k == float(0)):  #Probably a better way to sort sites from m_bond
                    i += 1
                    if bond.c == 1:
                        name = 'Harm_constrained'
                    else:
                        name = 'Harm'
                    dlines.append('      %d %d %d %s %10.8f %10.8f\n'
                                 %(i,
                                   bond.atom1,
                                   bond.atom2,
                                   name,
                                   length,
                                   k))
                elif not bond:
                    nonecnt+=1
            if nonecnt > 0 and verbose:
                print 'FOUND %d BONDS THAT DO NOT EXIST'%nonecnt
            header = "    ffio_bonds[%d] {\n" % (i)
            hlines = endheadersection(i==0,header,hlines)
            lines.extend(hlines)
            lines.extend(dlines)
            lines.append("      :::\n")
            lines.append("    }\n")

            #-ADDING ANGLES
            if verbose:
                print "   -Writing angles..."

            dlines = list()
            hlines = list()
            hlines.append("    ffio_angles_placeholder\n")
            hlines.append("      i_ffio_ai\n")
            hlines.append("      i_ffio_aj\n")
            hlines.append("      i_ffio_ak\n")
            hlines.append("      s_ffio_funct\n")
            hlines.append("      r_ffio_c1\n")
            hlines.append("      r_ffio_c2\n")
            hlines.append("      :::\n")
            i = 1
            for angle in moleculetype.angleForceSet.itervalues():
                if angle.c == 0:
                    dlines.append('      %d %d %d %d %s %10.8f %10.8f\n' % (i, angle.atom1, angle.atom2, angle.atom3, 'Harm', float(angle.theta.in_units_of(units.degrees)._value), float(angle.k.in_units_of(units.kilocalorie_per_mole/units.radians**2)._value)))
                elif angle.c == 1:
                    dlines.append('      %d %d %d %d %s %10.8f %10.8f\n' % (i, angle.atom1, angle.atom2, angle.atom3, 'Harm_constrained', float(angle.theta.in_units_of(units.degrees)._value), float(angle.k.in_units_of(units.kilocalorie_per_mole/units.radians**2)._value)))
                i+=1

            header = "    ffio_angles[%d] {\n" % (i-1)
            hlines = endheadersection(i==1,header,hlines)

            lines.extend(hlines)
            lines.extend(dlines)

            lines.append("      :::\n")
            lines.append("    }\n")


            #-ADDING DIHEDRALS
            if verbose:
                print "   -Writing dihedrals..."
            dlines = list()
            hlines = list()
            hlines.append("    ffio_dihedrals_placeholder\n")
            hlines.append("      i_ffio_ai\n")
            hlines.append("      i_ffio_aj\n")
            hlines.append("      i_ffio_ak\n")
            hlines.append("      i_ffio_al\n")
            hlines.append("      s_ffio_funct\n")

            i = 1
            # sorting for now to make it easier to debug
            #for dihedral in moleculetype.dihedralForceSet.itervalues():
            dihedrallist = sorted(moleculetype.dihedralForceSet.itervalues(), key=lambda x: x.atom1)
            for dihedral in dihedrallist:
                if isinstance(dihedral, ProperDihedral1) or isinstance(dihedral, ProperDihedral9):
                    if i == 1:
                        hlines.append("      r_ffio_c0\n")
                        hlines.append("      r_ffio_c1\n")
                        hlines.append("      r_ffio_c2\n")
                        hlines.append("      :::\n")
                    dlines.append('      %d %d %d %d %d %s %10.8f %10.8f %10.8f\n'%(i, dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4, 'Proper_Harm', float(dihedral.phi.in_units_of(units.degrees)._value), float(dihedral.k.in_units_of(units.kilocalorie_per_mole/units.radians**2)._value), int(dihedral.multiplicity)))
                elif isinstance(dihedral, ImproperDihedral2):
                    if i == 1:
                        hlines.append("      r_ffio_c0\n")
                        hlines.append("      r_ffio_c1\n")
                        hlines.append("      :::\n")
                    dlines.append('      %d %d %d %d %d %s %10.8f %10.8f\n'%(i, dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4, 'Improper_Harm', float(dihedral.xi.in_units_of(units.radians)._value), float(dihedral.k.in_units_of(units.kilocalorie_per_mole/units.radians**2)._value)))
                elif isinstance(dihedral, RBDihedral):
                    if i == 1:
                        hlines.append("      r_ffio_c0\n")
                        hlines.append("      r_ffio_c1\n")
                        hlines.append("      r_ffio_c2\n")
                        hlines.append("      r_ffio_c3\n")
                        hlines.append("      r_ffio_c4\n")
                        hlines.append("      r_ffio_c5\n")
                        hlines.append("      r_ffio_c6\n")
                        hlines.append("      r_ffio_c7\n")
                        hlines.append("      :::\n")
                    if dihedral.i == 1:
                        name = 'Improper_Trig'
                    else:
                        name = 'Proper_Trig'
                    # convert to the     
                    c0,c1,c2,c3,c4,c5,c6 = ConvertDihedralFromRBToProperTrig(dihedral.C0,dihedral.C1,dihedral.C2,
                                                                             dihedral.C3,dihedral.C4,dihedral.C5,
                                                                             dihedral.C6)
                    dlines.append('      %d %d %d %d %d %s 0.0 %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f\n' % (i, dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4, name, float(c0.in_units_of(units.kilocalories_per_mole)._value), float(c1.in_units_of(units.kilocalories_per_mole)._value), float(c2.in_units_of(units.kilocalories_per_mole)._value), float(c3.in_units_of(units.kilocalories_per_mole)._value), float(c4.in_units_of(units.kilocalories_per_mole)._value), float(c5.in_units_of(units.kilocalories_per_mole)._value), float(c6.in_units_of(units.kilocalories_per_mole)._value)))
                else:
                    print "ERROR (writeFile): found unsupported dihedral"
                i+=1
            header = "    ffio_dihedrals[%d] {\n" % (i-1)
            hlines = endheadersection(i==1,header,hlines)

            lines.extend(hlines)
            lines.extend(dlines)
            lines.append("      :::\n")
            lines.append("    }\n")

            #ADDING EXCLUSIONS
            i = 1
            if verbose:
                print "   -Writing exclusions..."
            bpos = len(lines) #storing position for exclusions instead of bonds
            hlines = list()
            dlines = list()
            hlines.append("    ffio_exclusions_placeholder\n")
            hlines.append("      i_ffio_ai\n")
            hlines.append("      i_ffio_aj\n")
            hlines.append("      :::\n")
            for exclusion in moleculetype.exclusions.itervalues():
                dlines.append('      %d %d %d\n'%(i, int(exclusion.exclusions[0]), int(exclusion.exclusions[1])))
                i+=1
            header = "    ffio_exclusions[%d] {\n"%(i-1)
            hlines = endheadersection(i==1,header,hlines)
            lines.extend(hlines)
            lines.extend(dlines)
            lines.append("      :::\n")
            lines.append("    }\n")

            #-ADDING PAIRS
            if verbose:
                print "   -Writing pairs..."
                
            dlines = list()
            hlines = list()
            hlines.append("ffio_pairs_placeholder\n")
            hlines.append("      i_ffio_ai\n")
            hlines.append("      i_ffio_aj\n")
            hlines.append("      s_ffio_funct\n")
            hlines.append("      r_ffio_c1\n")
            hlines.append("      :::\n")
            i = 1
            for pair in moleculetype.pairForceSet.itervalues():
                if re.match("LJ", pair.type):
                    dlines.append('      %d %d %d %s %10.8f\n' % (i, pair.atom1, pair.atom2, pair.type, System._sys._ljCorrection))
                elif re.match("Coulomb", pair.type):
                    dlines.append('      %d %d %d %s %10.8f\n' % (i, pair.atom1, pair.atom2, pair.type, System._sys._coulombCorrection))
                elif re.match("Both", pair.type):
                    dlines.append('      %d %d %d %s %10.8f\n' % (i, pair.atom1, pair.atom2, "LJ", System._sys._ljCorrection))
                    i+=1
                    dlines.append('      %d %d %d %s %10.8f\n' % (i, pair.atom1, pair.atom2, "Coulomb", System._sys._coulombCorrection))
                i+=1
            header = "    ffio_pairs[%d] {\n"%(i-1)
            hlines = endheadersection(i==1,header,hlines)

            lines.extend(hlines)
            lines.extend(dlines)
            lines.append("      :::\n")
            lines.append("    }\n")


            #ADDING RESTRAINTS
            # STILL NEED TO ADD
            
            #ADDING CONSTRAINTS
            if verbose:
                print "   -Writing constraints..."
            isHOH = False
            alen = 0
            clen = 0
            alen_max = 0
            clen_max = 0
            for constraint in moleculetype.constraints.itervalues():
                if re.search('AH',constraint.type):
                    alen = int(list(constraint.type)[-1])
                    clen = alen
                elif re.match('HOH',constraint.type):
                    alen = 2
                    clen = 3
                if alen_max < alen:
                    alen_max = alen
                    clen_max = clen
            # we now know the maximum length of all constraint types

            # not sure we need to sort these, but makes it easier to debug
            i = 0
            constraintlist = sorted(moleculetype.constraints.itervalues(),key=lambda x: x.atom1)
            dlines = list()
            hlines = list()

            for constraint in constraintlist: #calculate the max number of atoms in constraint
                i+=1
                if re.search('HOH',constraint.type):
                    cline = '      %d %d %d %d ' % (i,int(constraint.atom1),int(constraint.atom2),int(constraint.atom3))
                    for j in range(alen_max-3):
                        cline += '0 '
                    cline += constraint.type
                    cline += ' %10.8f' % (float(constraint.length1.in_units_of(units.degrees)._value))
                    cline += ' %10.8f' % (float(constraint.length2.in_units_of(units.angstroms)._value))
                    cline += ' %10.8f' % (float(constraint.length2.in_units_of(units.angstroms)._value))
                    for j in range(clen_max-3):
                        cline += ' <>'
                elif re.match('AH',constraint.type):
                    alen = int(list(constraint.type)[-1])
                    cline = '      %d ' % i
                    cline += ' %d ' % int(constraint.atom1)
                    cline += ' %d ' % int(constraint.atom2)
                    if alen > 1:
                        cline += ' %d ' % int(constraint.atom3)
                    if alen > 2:
                        cline += ' %d ' % int(constraint.atom4)
                    if alen > 3:
                        cline += ' %d ' % int(constraint.atom5)
                    if alen > 4:
                        cline += ' %d ' % int(constraint.atom6)
                    if alen > 5:
                        cline += ' %d ' % int(constraint.atom7)
                    if alen > 6:
                        cline += ' %d ' % int(constraint.atom8)
                    for j in range(alen,alen_max):
                        cline += ' 0 '
                    cline += constraint.type
                    cline += ' %10.8f' % (float(constraint.length1.in_units_of(units.angstroms)._value))
                    if alen > 1:
                        cline += ' %10.8f' % (float(constraint.length2.in_units_of(units.angstroms)._value))
                    if alen > 2:
                        cline += ' %10.8f' % (float(constraint.length3.in_units_of(units.angstroms)._value))
                    if alen > 3:
                        cline += ' %10.8f' % (float(constraint.length4.in_units_of(units.angstroms)._value))
                    if alen > 4:
                        cline += ' %10.8f' % (float(constraint.length5.in_units_of(units.angstroms)._value))
                    if alen > 5:
                        cline += ' %10.8f' % (float(constraint.length6.in_units_of(units.angstroms)._value))
                    if alen > 6:
                        cline += ' %10.8f' % (float(constraint.length7.in_units_of(units.angstroms)._value))
                    for j in range(alen,alen_max):
                        cline += ' 0.0'
                cline += '\n'
                dlines.append(cline)

            hlines.append("    ffio_constraints[%d] {\n"%(i))
            if (i==0):
                hlines.append("      :::\n")
            else:
                letters = ['i','j','k','l','m','n','o','p','q']
                for j in range(alen_max+1):
                    hlines.append('      i_ffio_a%s\n'%letters[j])
                hlines.append('      s_ffio_funct\n')
                for j in range(clen_max):
                    hlines.append('      r_ffio_c%d\n' %(j+1))
                hlines.append("      :::\n")
            lines.extend(hlines)
            lines.extend(dlines)
            lines.append("      :::\n")
            lines.append("    }\n")

            lines.append("  }\n")
            lines.append("}\n")

        fout = open(filename, 'w')
        for line in lines:
            fout.write(line)
        fout.close()
