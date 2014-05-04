import sys
import os
import re
import copy
import warnings
import math
from collections import OrderedDict

# for debugging eventually remove
import pdb

import intermol.unit as units
from collections import deque
from intermol.Atom import Atom
from intermol.Molecule import Molecule
from intermol.System import System
from intermol.Types import *
from intermol.Force import *
from intermol.HashMap import *

class GromacsTopologyParser(object):
    """
    A class containing methods required to read in a Gromacs(4.5.4) Topology File
    """
    _GroTopParser = None

    def __init__(self, defines=None):
        """
        Initializes a GromacsTopologyParse object which serves to read in a Gromacs
        topology into the abstract representation.

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

    def parseTopology(self, topfile, verbose=False):
        """
        Parses a Gromacs topology reading into the abstract

        Args:
            topfile: filename of the file to be parsed
            verbose: verbose output
        """
        lines = self.preprocess(topfile, verbose)
        if lines:
            self.readTopology(lines, verbose)

    def readTopology(self, expanded, verbose=False):
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
                    System._sys._nbFunc = int(fields[0])
                    System._sys._combinationRule = int(fields[1])
                    System._sys._genpairs = fields[2]
                    System._sys._ljCorrection = float(fields[3])
                    System._sys._coulombCorrection = float(fields[4])

                    expanded.pop(i)

                elif match.group('atomtypes'):
                    if verbose:
                        print "Parsing [ atomtypes ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newAtomType = None

                        # note -- we should be able to store either C6 and C12 parameters or sigma and epsilon: not both
                        atomtype = split[0].strip()
                        if len(split) == 7:  # atom name and bond type are the same, or there is no z.
                            d = 0  #offset
                            if (split[1].isdigit()):
                                Z = int(split[1])
                                bondtype = split[0].strip()
                            else:
                                Z = -1
                                bondtype = split[1].strip()
                        elif len(split) == 8: #atom and bond name and Z
                            d = 1
                            bondtype = split[1].strip()            # bondtype
                            Z = split[2]                           # Z
                        else:
                            print "ERROR: incorrect number of points in atomtype entry (%s)" % split

                        mass =  float(split[2+d]) * units.amu
                        charge = float(split[3+d]) * units.elementary_charge
                        ptype = split[4+d]
                        if System._sys._combinationRule == 1:
                            sigma = (float(split[6+d]) / float(split[5+d])) ** (1.0/6.0)
                            epsilon = float(split[5+d]) / (4*sigma**6)

                            newAtomType = AtomCR1Type(atomtype,
                                                      bondtype,
                                                      Z,
                                                      mass,
                                                      charge,
                                                      ptype,
                                                      sigma * units.kilojoules_per_mole * units.nanometers**(6),      # C6
                                                      epsilon * units.kilojoules_per_mole * units.nanometers**(12))   # C12

                        elif (System._sys._combinationRule == 2) or (System._sys._combinationRule == 3):
                            newAtomType = AtomCR23Type(atomtype,
                                                       bondtype,
                                                       Z,
                                                       mass,
                                                       charge,
                                                       ptype,
                                                       float(split[5+d]) * units.nanometers,           # sigma
                                                       float(split[6+d]) * units.kilojoules_per_mole)  # epsilon
                        System._sys._atomtypes.add(newAtomType)

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
                                    float(split[3]) * units.nanometers,
                                    float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-2))

                        # G96Bond
                        elif int(split[2]) == 2:
                            newBondType = G96BondType(split[0],
                                        split[1],
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-4))

                        # Morse
                        elif int(split[2]) == 3:
                            newBondType = MorseBondType(split[0],
                                        split[1],
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole,
                                        float(split[5]) * units.nanometers**(-1))

                        # Cubic
                        elif int(split[2]) == 4:
                            newBondType = CubicBondType(split[0],
                                        split[1],
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-2),
                                        float(split[5]) * units.kilojoules_per_mole * units.nanometers**(-3))

                        # Harmonic
                        elif int(split[2]) == 6:
                            newBondType = HarmonicBondType(split[0],
                                        split[1],
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-2))
                        else:
                            print "%s is not a supported bond type"

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
                                if System._sys._combinationRule == 1:
                                    newPairType = LJ1PairCR1Type(split[0],
                                            split[1],
                                            split[2],
                                            float(split[3]) * units.kilojoules_per_mole * units.nanometers**(6),
                                            float(split[4]) * units.kilojoules_per_mole * units.nanometers**(12))

                                elif System._sys._combinationRule == (2 or 3):
                                    newPairType = LJ1PairCR1Type(split[0],
                                            split[1],
                                            split[2],
                                            float(split[3]) * units.nanometers,
                                            float(split[4]) * units.kilojoules_per_mole)

                            # LJ/C. pair NB
                            elif len(split) == 7:
                                if System._sys._combinationRule == 1:
                                    newPairType = LJ1PairCR1Type(split[0],
                                            split[1],
                                            split[2],
                                            float(split[3]) * units.elementary_charge,
                                            float(split[4]) * units.elementary_charge,
                                            float(split[5]) * units.kilojoules_per_mole * units.nanometers**(6),
                                            float(split[6]) * units.kilojoules_per_mole * units.nanometers**(12))

                                elif System._sys._combinationRule == (2 or 3):
                                    newPairType = LJ1PairCR1Type(split[0],
                                            split[1],
                                            split[2],
                                            float(split[3]) * units.elementary_charge,
                                            float(split[4]) * units.elementary_charge,
                                            float(split[5]) * units.nanometers,
                                            float(split[6]) * units.kilojoules_per_mole)

                        # LJ/Coul. 1-4 (Type 2)
                        elif int(split[2]) == 2:
                            if System._sys._combinationRule == 1:
                                newPairType = LJ1PairCR1Type(split[0],
                                        split[1],
                                        split[2],
                                        split[3],
                                        float(split[4]) * units.elementary_charge,
                                        float(split[5]) * units.elementary_charge,
                                        float(split[6]) * units.kilojoules_per_mole * units.nanometers**(6),
                                        float(split[7]) * units.kilojoules_per_mole * units.nanometers**(12))

                            elif System._sys._combinationRule == (2 or 3):
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
                                    float(split[5]) * units.kilojoules_per_mole * units.nanometers**(-2),
                                    float(split[6]) * units.nanometers,
                                    float(split[7]) * units.kilojoules_per_mole * units.nanometers**(-2))

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
                        split = expanded.pop(i).split()
                        newDihedralType = None
                        # first, check whether they are using 2 or 4 atom types
                        if split[2].isdigit():
                            atom1 = 'X'
                            atom2 = split[0]
                            atom3 = split[1]
                            atom4 = 'X'
                            d = 0

                        elif split[4].isdigit():
                            atom1 = split[0]
                            atom2 = split[1]
                            atom3 = split[2]
                            atom4 = split[3]
                            d = 2

                        # We can fit everything into two types of dihedrals - dihedral_trig, and improper harmonic
                        # dihedral trig is of the form fc0 + sum_i=1^6 fci (cos(nx-phi)
                        # proper dihedrals can be stored easily in this form, since they have only 1 n
                        # improper dihedrals can as well (flag as improper)
                        # RB can be stored as well, assuming phi = 0 or 180
                        # Fourier can also be stored.
                        # a full dihedral trig can be decomposied in to multiple proper dihedrals.

                        # will need to handle this a little differently, in that we will need
                        # to add multiple 9 dihedrals together into a single dihedral_trig

                        if ((int(split[2+d]) == 1 or
                            int(split[2+d]) == 9 or
                            int(split[2+d]) == 4)
                            and (len(split) == 6+d)):

                            if int(split[2+d])==4:
                                improper = True
                            else:
                                improper = False
                            fc0,fc1,fc2,fc3,fc4,fc5,fc6 = ConvertDihedralFromProperDihedralToDihedralTrig(
                                float(split[4+d])*units.kilojoules_per_mole,int(split[5+d]))
                            newDihedralType = DihedralTrigType(
                                atom1, atom2, atom3, atom4, float(split[3+d]) * units.degrees,
                                fc0, fc1, fc2, fc3, fc4, fc5, fc6, improper = improper)

                        # Improper Harmonic Dihedral: type 2. Can't be converted to any other type
                        elif (int(split[2+d]) == 2) and (len(split) == 5+d):
                            newDihedralType = ImproperHarmonicDihedralType(
                                atom1,atom2,atom3,atom4,
                                float(split[3+d]) * units.degrees,
                                float(split[4+d]) * units.kilojoules_per_mole * units.radians**(-2))

                        # RBDihedral: type 3
                        elif (int(split[2+d]) == 3) and (len(split) == 9+d):
                            fc0, fc1, fc2, fc3, fc4, fc5, fc6 = ConvertDihedralFromRBToDihedralTrig(
                                float(split[3+d]) * units.kilojoules_per_mole,
                                float(split[4+d]) * units.kilojoules_per_mole,
                                float(split[5+d]) * units.kilojoules_per_mole,
                                float(split[6+d]) * units.kilojoules_per_mole,
                                float(split[7+d]) * units.kilojoules_per_mole,
                                float(split[8+d]) * units.kilojoules_per_mole,
                                0 * units.kilojoules_per_mole)
                            newDihedralType = DihedralTrigType(
                                atom1, atom2, atom3, atom4, 0 * units.degrees,
                                fc0, fc1, fc2, fc3, fc4, fc5, fc6)  # need to look at the sign here

                        # Fourier Dihedral 5
                        elif (int(split[2+d]) == 5) and (len(split) == 7+d):
                            fc0, fc1, fc2, fc3, fc4, fc5, fc6  = ConvertDihedralFromFourierToDihedralTrig(
                                float(split[3+d]) * units.kilojoules_per_mole,
                                float(split[4+d]) * units.kilojoules_per_mole,
                                float(split[5+d]) * units.kilojoules_per_mole,
                                float(split[6+d]) * units.kilojoules_per_mole)

                            newDihedralType = DihedralTrigType(
                                atom1, atom2, atom3, atom4, 0 * units.degrees,
                                fc0, fc1, fc2, fc3, fc4, fc5, fc6)

                        elif (int(split[2+d]) == 8):
                            print "Error: Tabulated dihedrals not supported"
                        else:
                            print "could not find dihedral type"

                        if newDihedralType:
                            self.dihedraltypes.add(newDihedralType)

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
                    System._sys.addMolecule(currentMolecule)
                    currentMoleculeType = System._sys._molecules[moleculeName]
                    currentMoleculeType.nrexcl = int(split[1])
                    expanded.pop(i)

                elif match.group('atoms'):
                    if verbose:
                        print "Parsing [ atoms ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()

                        atom = Atom(int(split[0]),          # AtomNum  (index)
                                split[4].strip(),           # atomName
                                int(split[2]),              # resNum
                                split[3].strip())           # resName
                        atom.setAtomType(0, split[1].strip())
                        atom.cgnr = int(split[5])
                        atom.setCharge(0, float(split[6]) * units.elementary_charge)
                        try:
                            atom.setMass(0, float(split[7]) * units.amu)
                        except:
                            atom.setMass(0, -1 * units.amu)

                        if len(split) == 11:
                            atom.setAtomType(1, split[8].strip())
                            atom.setCharge(1, float(split[9]) * units.elementary_charge)
                            atom.setMass(1, float(split[10]) * units.amu)

                        index = 0
                        for atomType in atom._atomtype:
                            # Searching for a matching atomType to pull values from
                            tempType = AbstractAtomType(atom._atomtype[index])
                            atomType = System._sys._atomtypes.get(tempType)
                            if atomType:
                                atom.Z = atomType.Z
                                if not atom.bondtype:
                                    if atomType.bondtype:
                                        atom.bondtype = atomType.bondtype
                                    else:
                                        sys.stderr.write("Warning: A suspicious parameter was found in atom/atomtypes. Visually inspect before using.\n")
                                if atom._mass[index]._value < 0:
                                    if atomType.mass._value >= 0:
                                        atom.setMass(index, atomType.mass)
                                    else:
                                        sys.stderr.write("Warning: A suspicious parameter was found in atom/atomtypes. Visually inspect before using.\n")
                                # Assuming ptype = A
                                #atom.ptype = atomType.ptype

                                atom.setSigma(index, atomType.sigma)
                                atom.setEpsilon(index, atomType.epsilon)

                            else:
                                sys.stderr.write("Warning: A corresponding AtomType was not found. Insert missing values yourself.\n")
                            index += 1

                        currentMolecule.addAtom(atom)
                elif match.group('bonds'):
                    if verbose:
                        print "Parsing [ bonds ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newBondForce = None

                        if len(split) == 3:
                            atomtype1 = currentMolecule._atoms[int(split[0])-1].bondtype
                            atomtype2 = currentMolecule._atoms[int(split[1])-1].bondtype
                            tempType = AbstractBondType(atomtype1, atomtype2)
                            bondType = self.bondtypes.get(tempType)
                            if not bondType:
                                # we only have the reversed bond order stored, flip the atoms
                                tempType = AbstractBondType(atomtype2, atomtype1)
                                bondType = self.bondtypes.get(tempType)
                            if not bondType:
                                raise Exception("Bondtype lookup failed for '{0}'".format(" ".join(split)))

                            if isinstance(bondType, BondType):
                                split.append(bondType.length)
                                split.append(bondType.k)

                            elif isinstance(bondType, G96BondType):
                                split.append(bondType.length)
                                split.append(bondType.l)

                            elif isinstance(bondType, CubicBondType):
                                split.append(bondType.length)
                                split.append(bondType.C2)
                                split.append(bondType.C3)

                            elif isinstance(bondType, MorseBondType):
                                split.append(bondType.length)
                                split.append(bondType.D)
                                split.append(bomdType.beta)

                            elif isinstance(bondType, HarmonicBondType):
                                split.append(bondType.length)
                                split.append(bondType.k)
                            else:
                                warnings.warn("Bondtype '{0}' is unsupported or something more complicated went wrong".format(bondType.type))

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
                        System._sys._forces.add(newBondForce)

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
                                newPairForce = AbstractPair(int(split[0]), int(split[1]), "Both")

                        currentMoleculeType.pairForceSet.add(newPairForce)
                        System._sys._forces.add(newPairForce)

                elif match.group('angles'):
                    if verbose:
                        print "Parsing [ angles ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newAngleForce = None

                        if len(split) == 4:
                            #atomtype1 = currentMolecule._atoms[int(split[0])-1].getAtomType()[0]
                            atomtype1 = currentMolecule._atoms[int(split[0])-1].bondtype
                            #atomtype2 = currentMolecule._atoms[int(split[1])-1].getAtomType()[0]
                            atomtype2 = currentMolecule._atoms[int(split[1])-1].bondtype
                            #atomtype3 = currentMolecule._atoms[int(split[2])-1].getAtomType()[0]
                            atomtype3 = currentMolecule._atoms[int(split[2])-1].bondtype
                            tempType = AbstractAngleType(atomtype1, atomtype2, atomtype3, split[3])
                            angleType = self.angletypes.get(tempType)
                            if not (angleType):
                                #flip it around.
                                tempType = AbstractAngleType(atomtype3, atomtype2, atomtype1, split[3])
                                angleType = self.angletypes.get(tempType)

                            if isinstance(angleType, AngleType):
                                split.append(angleType.theta)
                                split.append(angleType.k)

                            if isinstance(angleType, G96BondType):
                                split.append(angleType.theta)
                                split.append(angleType.k)

                            if isinstance(angleType, CrossBondBondAngleType):
                                split.append(angleType.r1)
                                split.append(angleType.r2)
                                split.append(angleType.k)

                            if isinstance(angleType, CrossBondAngleAngleType):
                                split.append(angleType.r1)
                                split.append(angleType.r2)
                                split.append(angleType.r3)
                                split.append(angleType.k)
                            if isinstance(angleType, UreyBradleyAngleType):
                                split.append(angleType.theta)
                                split.append(angleType.k)
                                split.append(angleType.r)
                                split.append(angleType.kUB)

                            if isinstance(angleType, QuarticAngleType):
                                split.append(angleType.theta)
                                split.append(angleType.C0)
                                split.append(angleType.C1)
                                split.append(angleType.C2)
                                split.append(angleType.C3)
                                split.append(angleType.C4)

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
                        System._sys._forces.add(newAngleForce)

                elif match.group('dihedrals'):
                    if verbose:
                        print "Parsing [ dihedrals ]..."
                    expanded.pop(i)
                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newDihedralForce = None
                        dihedralType = None

                        if len(split) == 5:

                            atomtype1 = currentMolecule._atoms[int(split[0])-1].bondtype
                            atomtype2 = currentMolecule._atoms[int(split[1])-1].bondtype
                            atomtype3 = currentMolecule._atoms[int(split[2])-1].bondtype
                            atomtype4 = currentMolecule._atoms[int(split[3])-1].bondtype

                            # check the possible ways to match a dihedraltype
                            atomtypelists = [[atomtype1, atomtype2, atomtype3, atomtype4],  # original order
                                             [atomtype4, atomtype3, atomtype2, atomtype1],  # flip it
                                             [atomtype1, atomtype2, atomtype3, 'X'],  #single wildcard 1
                                             ['X', atomtype2, atomtype3, atomtype4],  #single wildcard 2
                                             ['X', atomtype3, atomtype2, atomtype1], # flipped single wildcard 1
                                             [atomtype4, atomtype3, atomtype2, 'X'], # flipped single wildcard 2
                                             ['X', atomtype2, atomtype3, 'X'], # double wildcard
                                             ['X', atomtype3, atomtype2, 'X'] # flipped double wildcard
                                             ]

                            for alist in atomtypelists:
                                if not dihedralType:
                                    tempType = AbstractDihedralType(alist[0], alist[1], alist[2], alist[3])
                                    dihedralType = self.dihedraltypes.get(tempType)
                                else:
                                    break

                        if isinstance(dihedralType, DihedralTrigType):
                            phi = dihedralType.phi
                            fc0 = dihedralType.fc0
                            fc1 = dihedralType.fc1
                            fc2 = dihedralType.fc2
                            fc3 = dihedralType.fc3
                            fc4 = dihedralType.fc4
                            fc5 = dihedralType.fc5
                            fc6 = dihedralType.fc6

                        if isinstance(dihedralType, ImproperHarmonicDihedralType):
                            xi = dihedralType.xi
                            k = dihedralType.k

                        atom1 = int(split[0])
                        atom2 = int(split[1])
                        atom3 = int(split[2])
                        atom4 = int(split[3])

                        # Proper Dihedral 1
                        if int(split[4]) == 1 or int(split[4]) == 9 or int(split[4]) == 4:

                            if int(split[4]) == 4:
                                improper = True
                            else:
                                improper = False

                            if len(split) > 5:
                                fc0, fc1, fc2, fc3, fc4, fc5, fc6 = ConvertDihedralFromProperDihedralToDihedralTrig(
                                    float(split[6]) * units.kilojoules_per_mole, split[7])
                                phi = float(split[5]) * units.degrees

                            newDihedralForce = DihedralTrigDihedral(
                                atom1, atom2, atom3, atom4,
                                phi, fc0, fc1, fc2, fc3, fc4, fc5, fc6, improper = improper)
                            # for dihedral type 9, there can be multiple interactions

                        # Improper Dihedral 2

                        elif int(split[4]) == 2:
                            if len(split) > 5:
                                phi = float(split[5]) * degrees
                                k = float(split[6]) * units.kilojoules_per_mole * units.radians**(-2)

                            newDihedralForce = ImproperHarmonicDihedral(
                                atom1, atom2, atom3, atom4, xi, k)

                        # RBDihedral
                        elif int(split[4]) == 3:

                            if len(split) > 5:
                                fc0, fc1, fc2, fc3, fc4, fc5, fc6 = ConvertDihedralFromRBToDihedralTrig(
                                    float(split[4]) * units.kilojoules_per_mole,
                                    float(split[5]) * units.kilojoules_per_mole,
                                    float(split[6]) * units.kilojoules_per_mole,
                                    float(split[7]) * units.kilojoules_per_mole,
                                    float(split[8]) * units.kilojoules_per_mole,
                                    float(split[9]) * units.kilojoules_per_mole,
                                    float(split[10]) * units.kilojoules_per_mole)

                            newDihedralForce = DihedralTrigDihedral(
                                atom1, atom2, atom3, atom4,
                                0 * units.degrees, fc0, fc1, fc2, fc3, fc4, fc5, fc6)  # need to look at the use of sign here

                        elif int(split[4]) == 5:

                            if len(split) > 5:
                                fc0, fc1, fc2, fc3, fc4, fc5, fc6 = ConvertDihedralFromFourierToDihedralTrig(
                                    float(split[5])*units.kilojoules_per_mole,
                                    float(split[6])*units.kilojoules_per_mole,
                                    float(split[7])*units.kilojoules_per_mole,
                                    float(split[8])*units.kilojoules_per_mole,
                                    )

                        elif int(split[4]) == 8:
                            print "Error: Cannot support tabulated dihedrals"
                        else:
                            print "ERROR: unsupported dihedral"

                        currentMoleculeType.dihedralForceSet.add(newDihedralForce)
                        System._sys._forces.add(newDihedralForce)

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
                        System._sys._forces.add(newSettlesForce)

                        # we need to add a constrainted bonded force as well between the atoms in these molecules. 
                        # we assume the gromacs default of 1. O, 2. H, 3. H
                        # reference bond strength is 900 kj/mol, but doesn't really matter since constrainted.
                        waterbondrefk = 900*units.kilojoules_per_mole * units.nanometers**(-2)
                        wateranglerefk = 400*units.kilojoules_per_mole * units.degrees**(-2)
                        angle = 2.0 * math.asin(0.5 * float(split[3]) / float(split[2])) * units.radians
                        dOH = float(split[2]) * units.nanometers
                        
                        newBondForce = Bond(1,2,dOH,waterbondrefk,c=True)
                        currentMoleculeType.bondForceSet.add(newBondForce)
                        System._sys._forces.add(newBondForce)

                        newBondForce = Bond(1,3,dOH,waterbondrefk,c=True)
                        currentMoleculeType.bondForceSet.add(newBondForce)
                        System._sys._forces.add(newBondForce)

                        newAngleForce = Angle(3,1,2,angle,wateranglerefk,c=True)
                        currentMoleculeType.angleForceSet.add(newAngleForce)
                        System._sys._forces.add(newAngleForce)

                elif match.group('exclusions'):
                    if verbose:
                        print "Parsing [ exclusions ]..."
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        for j in range(len(split)):
                            if split[0] < split[j]:
                                newExclusion = Exclusions([split[0],split[j]])
                                currentMoleculeType.exclusions.add(newExclusion)
                                System._sys._forces.add(newExclusion)

                elif match.group('molecules'):
                    if verbose:
                        print "Parsing [ molecules ]..."
                    expanded.pop(i)
                    ordered_moleculetypes = OrderedDict()
                    while i < len(expanded) and not (expanded[i].count('[')):
                        split = expanded.pop(i).split()
                        mol_name = split[0]
                        mol_num = int(split[1])
                        System._sys._components.append((mol_name, mol_num))

                        ordered_moleculetypes[mol_name] = System._sys._molecules[mol_name]
                        tempMolecule = System._sys._molecules[mol_name].moleculeSet[0]
                        if len(System._sys._molecules[mol_name].moleculeSet) > 1:
                            n = 0
                        else:
                            n = 1
                        while n < mol_num:
                            mol = copy.deepcopy(tempMolecule)
                            System._sys.addMolecule(mol)
                            n += 1
                    System._sys._molecules = ordered_moleculetypes
                else:
                    i += 1
            else:
                i += 1

    def preprocess(self, filename, verbose=False):
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
        if fd is not None:
            lines = deque(fd)   # initial read of root file
            fd.close()
            while(lines):       # while the queue isn't empty continue preprocessing
                line = lines.popleft()
                if verbose:
                    print "========================="
                    print "Original line:", line
                match = lineReg.match(line)   # ensure that the line matches what we expect!
                if match:
                    line = match.group('directive')
                    comment = match.group('comment')
                    if verbose:
                        print "Directive:", line
                        print "Comment:", comment
                    while True:     # expand the lines
                        if match.group('ignore_newline'):
                            if verbose:
                                print "Ignore newline found"
                            if lines:
                                temp = lines.popleft()
                                match = lineReg.match(temp)
                                line += match.group('directive')
                            else:
                                print "WARNING: Previous line continues yet EOF encountered"
                                break   # EOF so we can't loop again
                        else:
                            break   # line terminates
                    line = line.strip()  # just in case theres whitespace we missed

                    match = preReg.match(line)
                    if match is not None:
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
                                if fd is not None:
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

    def open(self, filename, parentfile=''):
        """
        Open a file and add to includes list

        Checks if a file exists in the same directory as the top and gro files
        or in the GMXLIB environment; then opens if found or fails if not

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
                sys.stderr.write("Local instance of topology file '%s' used\n" % (filename))
                return fd
            except IOError, (errno, strerror):
                sys.stderr.write("I/O error(%d): %s for local instance of '%s'\n" % (errno, strerror, filename))

        # if we can't include the file locally, let's look for it in
        # the directory of one of the previous files.  

        # NOTE: could be
        # dangerous if multiple files with the same name are found in
        # different directories!
        found = False
        for parentfile in self.includes:
            parentdir, oldfile = os.path.split(parentfile)
            filename = os.path.join(parentdir, temp)
            if os.path.exists(filename):
                found = True
                break

        if found:
            try:
                fd = open(filename)
                self.includes.add(temp)
                sys.stderr.write("version in %s used for topology file '%s'\n" % (parentdir, filename))
                return fd
            except:
                pass
        else:
            try:
                filename = os.path.join(os.environ['GMXLIB'], temp)
                fd = open(filename)
                self.includes.add(filename)
                sys.stderr.write("version in GMXLIB = %s used for topology file '%s'\n" % (os.environ['GMXLIB'], temp))
                return fd

            except IOError, (errno, strerror):
                sys.stderr.write("Unrecoverable I/O error(%d): can'd find %s anywhere: '%s'\n" % (errno, strerror, filename))
                sys.exit()

    def writeTopology(self, filename):
        """Write this topology in GROMACS file format.

        Args:
            filename: the name of the file to write out to
        """
        lines = list()

        # [ defaults ]
        lines.append('[ defaults ]\n')
        lines.append('; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n')
        lines.append('%d%16d%18s%20.4f%8.4f\n'
                % (System._sys._nbFunc,
                   System._sys._combinationRule,
                   System._sys._genpairs,
                   System._sys._ljCorrection,
                   System._sys._coulombCorrection))
        lines.append('\n')

        # [ atomtypes ]
        lines.append('[ atomtypes ]\n')
        lines.append(';type, bondtype, Z, mass, charge, ptype, sigma, epsilon\n')
        for atomtype in System._sys._atomtypes.itervalues():
            if atomtype.atomtype.isdigit():
                atomtype.atomtype = "LMP_" + atomtype.atomtype
            if atomtype.bondtype.isdigit():
                atomtype.bondtype = "LMP_" + atomtype.bondtype
            lines.append('%-11s%5s%6d%18.8f%18.8f%5s%18.8e%18.8e\n'
                    % (atomtype.atomtype,
                       atomtype.bondtype,
                       int(atomtype.Z),
                       atomtype.mass.in_units_of(units.atomic_mass_unit)._value,
                       atomtype.charge.in_units_of(units.elementary_charge)._value,
                       atomtype.ptype,
                       atomtype.sigma.in_units_of(units.nanometers)._value,
                       atomtype.epsilon.in_units_of(units.kilojoules_per_mole)._value))
        lines.append('\n')

        # [ moleculetype]
        for moleculeType in System._sys._molecules.itervalues():
            lines.append('[ moleculetype ]\n')
            lines.append('%s%10d\n\n'
                    % (moleculeType.name,
                       moleculeType.nrexcl))

            # [ atoms ]
            lines.append('[ atoms ]\n')
            lines.append(';num, type, resnum, resname, atomname, cgnr, q, m\n')
            molecule = moleculeType.moleculeSet[0]
            count = 1
            for atom in molecule._atoms:
                if atom.atomName.isdigit():
                    atom.atomName = "LMP_" + atom.atomName
                if atom._atomtype[0].isdigit():
                    atom._atomtype[0] = "LMP_" + atom._atomtype[0]

                try:
                    lines.append('%6d%18s%6d%8s%8s%6d%18.8f%18.8f%18s%18.8f%18.8f\n'
                            % (count,
                               atom._atomtype[0],
                               atom.residueIndex,
                               atom.residueName,
                               atom.atomName,
                               atom.cgnr,
                               atom._charge[0].in_units_of(units.elementary_charge)._value,
                               atom._mass[0].in_units_of(units.atomic_mass_unit)._value,
                               atom._atomtype[1],
                               atom._charge[1].in_units_of(units.elementary_charge)._value,
                               atom._mass[1].in_units_of(units.atomic_mass_unit)._value))
                except:
                    lines.append('%6d%18s%6d%8s%8s%6d%18.8f%18.8f\n'
                                 % (count,
                                    atom._atomtype[0],
                                    atom.residueIndex,
                                    atom.residueName,
                                    atom.atomName,
                                    atom.cgnr,
                                    atom._charge[0].in_units_of(units.elementary_charge)._value,
                                    atom._mass[0].in_units_of(units.atomic_mass_unit)._value))
                count += 1
            lines.append('\n')

            if moleculeType.bondForceSet:
                # [ bonds ]
                lines.append('[ bonds ]\n')
                lines.append(';   ai     aj funct  r               k\n')
                for bond in moleculeType.bondForceSet.itervalues():
                    if isinstance(bond, Bond):
                        b_type = 1
                        lines.append('%6d%7d%4d%18.8e%18.8e\n'
                                % (bond.atom1,
                                   bond.atom2,
                                   b_type,
                                   bond.length.in_units_of(units.nanometers)._value,
                                   bond.k.in_units_of(units.kilojoules_per_mole*units.nanometers**(-2))._value))
                    elif isinstance(bond, G96Bond):
                        b_type = 2
                        lines.append('%6d%7d%4d%5.8f%5.8f\n'
                                % (bond.atom1,
                                   bond.atom2,
                                   b_type,
                                   bond.length.in_units_of(units.nanometers)._value,
                                   bond.k.in_units_of(units.kilojoules_per_mole*units.nanometers**(-2))._value))
                    else:
                        print "ERROR (writeTopology): found unsupported bond type"
                lines.append('\n')

            if moleculeType.pairForceSet:
                # [ pair ]
                lines.append('[ pairs ]\n')
                lines.append(';  ai    aj   funct\n')
                for pair in moleculeType.pairForceSet.itervalues():

                    if isinstance(pair, AbstractPair):
                        p_type = 1
                        lines.append('%6d%7d%4d\n'
                                % (pair.atom1,
                                   pair.atom2,
                                   p_type))

                    else:
                        print "ERROR (writeTopology): found unsupported pair type"
                lines.append('\n')

            if moleculeType.angleForceSet:
                # [ angles ]
                lines.append('[ angles ]\n')
                lines.append(';   ai     aj     ak    funct   theta         cth\n')
                for angle in moleculeType.angleForceSet.itervalues():
                    if isinstance(angle, Angle):
                        a_type = 1
                        atomindex = "%6d%7d%7d" % (angle.atom1,angle.atom2,angle.atom3)
                        lines.append('%s%4d%18.8e%18.8e\n'
                                % (atomindex,
                                   a_type,
                                   angle.theta.in_units_of(units.degrees)._value,
                                   angle.k.in_units_of(units.kilojoules_per_mole*units.radians**(-2))._value))
                    elif isinstance(angle, UreyBradleyAngle):
                        lines.append('%s%4d%18.8e%18.8e%18.8e%18.8e\n'
                                     % (atomindex,
                                        a_type,
                                        angle.theta.in_units_of(units.degrees)._value,
                                        angle.k.in_units_of(units.kilojoules_per_mole*units.radians**(-2))._value,
                                        angle.r.in_units_of(units.angstroms)._value,
                                        angle.kUB.in_units_of(units.kilojoules_per_mole)._value))
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
                    lines.append('%6d%7d%4d%18.8e%18.8e\n'%(pair.atom1, pair.atom2, type, pair.V.in_units_of(units.XXX)._value, pair.W.in_units_of(units.XXX)._value)

                elif isinstance(pair, LJ2PairCR1) or isinstance(pair, LJ2PairCR23):
                    type = 2
                    lines.append('%6d%7d%4d%18.8e%18.8e\n'%(pair.atom1, pair.atom2, type, pair.V.in_units_of(units.XXX)._value, pair.W.in_units_of(units.XXX)._value))

                elif isinstance(pair, LJNBCR1) or isinstance( pair, LJNBCR23):
                    type = 1
                    lines.append('%6d%7d%4d%18.8f%18.8f%18.8f%18.8f\n'%(pair.atom1, pair.atom2, type, pair.qi.in_units_of(units.XXX)._value, pair.qj.in_units_of(units.XXX)._value, pair.V.in_units_of(units.XXX)._value, pair.W.in_units_of(units.XXX)._value))
                else:
                    print "Could not identify pair!"
            """

            if moleculeType.dihedralForceSet:
                # [ dihedrals ]
                lines.append('[ dihedrals ]\n')
                lines.append(';    i      j      k      l   func\n')

                for dihedral in moleculeType.dihedralForceSet.itervalues():
                    # this atom index will be the same for all of types.
                    atomindex = "%7d%7d%7d%7d" % (dihedral.atom1,dihedral.atom2,dihedral.atom3,dihedral.atom4)
                    if isinstance(dihedral, DihedralTrigDihedral):
                        # convienience array
                        darray = [dihedral.fc1,dihedral.fc2,dihedral.fc3,dihedral.fc4,dihedral.fc5,dihedral.fc6]
                        if (dihedral.improper):
                            for j in range(6):  # only one of these should be nonzero
                                if darray[j].in_units_of(units.kilojoules_per_mole)._value != 0:
                                    lines.append('%s%4d%18.8f%18.8f%6d\n'
                                                 % (atomindex, 4,
                                                    dihedral.phi.in_units_of(units.degrees)._value,
                                                    darray[j].in_units_of(units.kilojoules_per_mole)._value,
                                                    j+1))
                        else:
                            if (dihedral.phi==0*units.degrees or dihedral.phi==180*units.degrees):
                                d_type = 3
                                c = ConvertDihedralFromDihedralTrigToRB(
                                    math.cos(dihedral.phi.in_units_of(units.radians)._value),
                                    dihedral.phi,
                                    dihedral.fc0,
                                    dihedral.fc1,
                                    dihedral.fc2,
                                    dihedral.fc3,
                                    dihedral.fc4,
                                    dihedral.fc5,
                                    dihedral.fc6)
                                if (c[6]._value != 0):
                                    print "ERROR: Gromacs does not handle multiplicities of greater than 6"
                                lines.append('%s%4d%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f\n'
                                             % (atomindex,
                                                d_type,
                                                c[0].in_units_of(units.kilojoules_per_mole)._value,
                                                c[1].in_units_of(units.kilojoules_per_mole)._value,
                                                c[2].in_units_of(units.kilojoules_per_mole)._value,
                                                c[3].in_units_of(units.kilojoules_per_mole)._value,
                                                c[4].in_units_of(units.kilojoules_per_mole)._value,
                                                c[5].in_units_of(units.kilojoules_per_mole)._value))
                            else:
                                #print as a type 1 dihedral, or a series of type 9 dihedrals
                                ncount = 0

                                for j in range(6):
                                    if float(darray[j]._value) != 0:
                                        ncount +=1
                                if ncount > 1:
                                    dtype = 9
                                else:
                                    dtype = 1
                                for j in range(6):
                                    if float(darray[j]._value) != 0:
                                        lines.append('%s%4d%18.8f%18.8f%6d\n'
                                        % (atomindex,
                                           dtype,
                                           dihedral.phi.in_units_of(units.degrees)._value,
                                           darray[j].in_units_of(units.kilojoules_per_mole)._value,
                                           j+1))

                    elif isinstance(dihedral, ImproperHarmonicDihedral):
                        d_type = 2
                        lines.append('%s%4d%18.8f%18.8f\n'
                                     % (atomindex,
                                        d_type,
                                        dihedral.xi.in_units_of(units.degrees)._value,
                                        dihedral.k.in_units_of(units.kilojoules_per_mole*units.radians**(-2))._value))

                    else:
                        print "ERROR (writeTopology): found unsupported  dihedral type"
                lines.append('\n')

            if moleculeType.settles:
                settles = moleculeType.settles
                # [ settles ]
                lines.append('[ settles ]\n')
                lines.append('; i  funct   dOH  dHH\n')
                s_type = 1
                lines.append('%6d%6d%18.8f%18.8f\n'
                             % (settles.atom1,
                                s_type,
                                settles.dOH.in_units_of(units.nanometers)._value,
                                settles.dHH.in_units_of(units.nanometers)._value))
                lines.append('\n')

            if moleculeType.exclusions:
                # [ exclusions ]
                lines.append('[ exclusions ]\n')
                for i, exclusion in enumerate(moleculeType.exclusions.itervalues()):
                    lines.append('%6s%6s\n'
                                 % (exclusion.exclusions[0],
                                    exclusion.exclusions[1]))

        # [ system ]
        lines.append('[ system ]\n')
        lines.append('%s\n' % (System._sys._name))
        lines.append('\n')

        # [ molecules ]
        lines.append('[ molecules ]\n')
        lines.append('; Compound        nmols\n')
        #for component in System._sys._components:
        #    lines.append('%-15s%8d\n'
        #                 % (component[0],
        #                    component[1]))
        #keeping this for now, since we don't know when it might be preferable.
        # The following lines are more 'chemical'
        for molType in System._sys._molecules:
            lines.append('%-15s%8d\n'
                    % (molType,
                      len(System._sys._molecules[molType].moleculeSet)))

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
