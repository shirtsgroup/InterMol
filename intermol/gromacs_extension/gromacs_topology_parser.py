import sys
import os
import re
import copy
from warnings import warn
from collections import OrderedDict
import logging
import numpy as np

import intermol.unit as units
from collections import deque
from intermol.atom import Atom
from intermol.molecule import Molecule
from intermol.system import System
from intermol.types import *
from intermol.forces import *
from intermol.hashmap import *
import math

logger = logging.getLogger('InterMolLog')

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

    def parse_topology(self, topfile, verbose=False):
        """
        Parses a Gromacs topology reading into the abstract

        Args:
            topfile: filename of the file to be parsed
            verbose: verbose output
        """
        lines = self.preprocess(topfile, verbose)
        if lines:
            self.read_topology(lines, verbose)

    def read_topology(self, expanded, verbose=False):
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
                logger.debug(match.groups())
                if match.group('defaults'):
                    logger.debug("Parsing [ defaults ]...")
                    expanded.pop(i)

                    fields = expanded[i].split()
                    System._sys.nonbonded_function = int(fields[0])
                    System._sys.combination_rule = int(fields[1])
                    System._sys.genpairs = fields[2]
                    System._sys.lj_correction = float(fields[3])
                    System._sys.coulomb_correction = float(fields[4])

                    expanded.pop(i)

                elif match.group('atomtypes'):
                    logger.debug("Parsing [ atomtypes ]...")
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newAtomType = None

                        # note -- we should be able to store either C6 and C12 parameters or sigma and epsilon: not both
                        atomtype = split[0].strip()
                        if len(split) == 7:  # atom name and bond type are the same, or there is no z.
                            d = 0  #offset
                            if (split[1].isdigit()):
                                atomic_number = int(split[1])
                                bondtype = split[0].strip()
                            else:
                                atomic_number = -1
                                bondtype = split[1].strip()
                        elif len(split) == 8: #atom and bond name and atomic_number
                            d = 1
                            bondtype = split[1].strip()            # bondtype
                            atomic_number = split[2]                           # atomic_number
                        else:
                            raise Exception("Incorrect number of points in atomtype entry (%s)" % split)

                        mass =  float(split[2+d]) * units.amu
                        charge = float(split[3+d]) * units.elementary_charge
                        ptype = split[4+d]
                        if System._sys.combination_rule == 1:
                            sigma = (float(split[6+d]) / float(split[5+d])) ** (1.0/6.0)
                            epsilon = float(split[5+d]) / (4*sigma**6)

                            newAtomType = AtomCR1Type(atomtype,
                                                      bondtype,
                                                      atomic_number,
                                                      mass,
                                                      charge,
                                                      ptype,
                                                      sigma * units.kilojoules_per_mole * units.nanometers**(6),      # C6
                                                      epsilon * units.kilojoules_per_mole * units.nanometers**(12))   # C12

                        elif (System._sys.combination_rule == 2) or (System._sys.combination_rule == 3):
                            newAtomType = AtomCR23Type(atomtype,
                                                       bondtype,
                                                       atomic_number,
                                                       mass,
                                                       charge,
                                                       ptype,
                                                       float(split[5+d]) * units.nanometers,           # sigma
                                                       float(split[6+d]) * units.kilojoules_per_mole)  # epsilon
                        System._sys._atomtypes.add(newAtomType)

                elif match.group('bondtypes'):
                    logger.debug("Parsing [ bondtypes ]...")
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
                            newBondType = HarmonicPotentialType(split[0],
                                        split[1],
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-2))
                        else:
                            raise Exception("%s is not a supported bond type" % split[2])

                        if newBondType:
                            self.bondtypes.add(newBondType)

                elif match.group('pairtypes'):
                    logger.debug("Parsing [ pairtypes ]...")
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:

                        split = expanded.pop(i).split()
                        newPairType = None
                        if int(split[2]) == 1:
                            # LJ/Coul. 1-4 (Type 1)
                            if len(split) == 5:
                                if System._sys.combination_rule == 1:
                                    newPairType = LJ1PairCR1Type(split[0],
                                            split[1],
                                            split[2],
                                            float(split[3]) * units.kilojoules_per_mole * units.nanometers**(6),
                                            float(split[4]) * units.kilojoules_per_mole * units.nanometers**(12))

                                elif System._sys.combination_rule == (2 or 3):
                                    newPairType = LJ1PairCR1Type(split[0],
                                            split[1],
                                            split[2],
                                            float(split[3]) * units.nanometers,
                                            float(split[4]) * units.kilojoules_per_mole)

                            # LJ/C. pair NB
                            elif len(split) == 7:
                                if System._sys.combination_rule == 1:
                                    newPairType = LJ1PairCR1Type(split[0],
                                            split[1],
                                            split[2],
                                            float(split[3]) * units.elementary_charge,
                                            float(split[4]) * units.elementary_charge,
                                            float(split[5]) * units.kilojoules_per_mole * units.nanometers**(6),
                                            float(split[6]) * units.kilojoules_per_mole * units.nanometers**(12))

                                elif System._sys.combination_rule == (2 or 3):
                                    newPairType = LJ1PairCR1Type(split[0],
                                            split[1],
                                            split[2],
                                            float(split[3]) * units.elementary_charge,
                                            float(split[4]) * units.elementary_charge,
                                            float(split[5]) * units.nanometers,
                                            float(split[6]) * units.kilojoules_per_mole)

                        # LJ/Coul. 1-4 (Type 2)
                        elif int(split[2]) == 2:
                            if System._sys.combination_rule == 1:
                                newPairType = LJ1PairCR1Type(split[0],
                                        split[1],
                                        split[2],
                                        split[3],
                                        float(split[4]) * units.elementary_charge,
                                        float(split[5]) * units.elementary_charge,
                                        float(split[6]) * units.kilojoules_per_mole * units.nanometers**(6),
                                        float(split[7]) * units.kilojoules_per_mole * units.nanometers**(12))

                            elif System._sys.combination_rule == (2 or 3):
                                newPairType = LJ1PairCR1Type(split[0],
                                         split[1],
                                         split[2],
                                         split[3],
                                         float(split[4]) * units.elementary_charge,
                                         float(split[5]) * units.elementary_charge,
                                         float(split[6]) * units.nanometers,
                                         float(split[7]) * units.kilojoules_per_mole)

                        else:
                            raise Exception("Could not find pair type")

                        if newPairType:
                            self.pairtypes.add(newPairType)

                elif match.group('angletypes'):
                    logger.debug("Parsing [ angletypes ]...")
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newAngleType = None

                        # Angle
                        if int(split[3]) == 1:
                            newAngleType = AngleType(split[0],
                                    split[1],
                                    split[2],
                                    float(split[4]) * units.degrees,
                                    float(split[5]) * units.kilojoules_per_mole * units.radians**(-2))

                        # G96Angle
                        elif int(split[3]) == 2:
                            newAngleType = G96AngleType(split[0],
                                    split[1],
                                    split[2],
                                    float(split[4]) * units.degrees,
                                    float(split[5]) * units.kilojoules_per_mole)

                        # Cross bond-bond
                        elif int(split[3]) == 3:
                            newAngleType = CrossBondBondAngleType(split[0],
                                    split[1],
                                    split[2],
                                    float(split[4]) * units.nanometers,
                                    float(split[5]) * units.nanometers,
                                    float(split[6]) * units.kilojoules_per_mole * units.nanometers**(-2))

                        # Cross bond-angle
                        elif int(split[3]) == 4:
                            newAngleType = CrossBondAngleAngleType(split[0],
                                    split[1],
                                    split[2],
                                    float(split[4]) * units.nanometers,
                                    float(split[5]) * units.nanometers,
                                    float(split[6]) * units.nanometers,
                                    float(split[7]) * units.kilojoules_per_mole * units.nanometers**(-2))

                        # Urey-Bradley
                        elif int(split[3]) == 5:
                            newAngleType = UreyBradleyAngleType(split[0],
                                    split[1],
                                    split[2],
                                    float(split[4]) * units.degrees,
                                    float(split[5]) * units.kilojoules_per_mole * units.radians**(-2),
                                    float(split[6]) * units.nanometers,
                                    float(split[7]) * units.kilojoules_per_mole * units.nanometers**(-2))

                        # Quartic
                        elif int(split[3]) == 6:
                            newAngleType = QuarticAngleType(split[0],
                                    split[1],
                                    split[2],
                                    float(split[4]) * units.degrees,
                                    float(split[5]) * units.kilojoules_per_mole,
                                    float(split[6]) * units.kilojoules_per_mole * units.radians**(-1),
                                    float(split[7]) * units.kilojoules_per_mole * units.radians**(-2),
                                    float(split[8]) * units.kilojoules_per_mole * units.radians**(-3),
                                    float(split[9]) * units.kilojoules_per_mole * units.radians**(-4))

                        else:
                            raise Exception("Could not find angle type")

                        if newAngleType:
                            self.angletypes.add(newAngleType)

                elif match.group('dihedraltypes'):
                    logger.debug("Parsing [ dihedraltypes ]...")
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
                        # to add multiple 9 dihedrals together into a single dihedral_trig, as long as they
                        # have the same phi angle (seems to be always the case).

                        dtype = int(split[2+d])
                        nentries = len(split)

                        if (dtype == 1 or dtype == 4 or dtype == 9) and nentries == 6+d:

                            if dtype == 4:
                                improper = True
                            else:
                                improper = False

                            fc0,fc1,fc2,fc3,fc4,fc5,fc6 = ConvertDihedralFromProperDihedralToDihedralTrig(
                                float(split[4+d])*units.kilojoules_per_mole,int(split[5+d]))
                            newDihedralType = DihedralTrigType(
                                atom1, atom2, atom3, atom4, float(split[3+d]) * units.degrees,
                                fc0, fc1, fc2, fc3, fc4, fc5, fc6, improper = improper)

                        # Improper Harmonic Dihedral: type 2. Can't be converted to any other type
                        elif (dtype == 2) and nentries == 5+d:
                            newDihedralType = ImproperHarmonicDihedralType(
                                atom1,atom2,atom3,atom4,
                                float(split[3+d]) * units.degrees,
                                float(split[4+d]) * units.kilojoules_per_mole * units.radians**(-2))

                        # RBDihedral: type 3
                        elif (dtype == 3) and nentries == 9+d:
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
                        elif dtype == 5 and nentries == 7+d:
                            fc0, fc1, fc2, fc3, fc4, fc5, fc6  = ConvertDihedralFromFourierToDihedralTrig(
                                float(split[3+d]) * units.kilojoules_per_mole,
                                float(split[4+d]) * units.kilojoules_per_mole,
                                float(split[5+d]) * units.kilojoules_per_mole,
                                float(split[6+d]) * units.kilojoules_per_mole)

                            newDihedralType = DihedralTrigType(
                                atom1, atom2, atom3, atom4, 0 * units.degrees,
                                fc0, fc1, fc2, fc3, fc4, fc5, fc6)

                        elif dtype == 8:
                            raise Exception("Tabulated dihedrals not supported")
                        else:
                            raise Exception("Could not find dihedral type")

                        if newDihedralType:
                            if dtype == 9:
                                # we can't actually store multiple dihedral parameters in our
                                # architecture, so we add up the types into a single DihedralTrigDihedral angle
                                try:
                                    dihedralmatch = self.dihedraltypes.get(newDihedralType)
                                    dihedralmatch.sum_parameters(newDihedralType)
                                except Exception as e:
                                    logger.exception(e) # EDZ: used to be pass, now recorded but supressed
                            self.dihedraltypes.add(newDihedralType)

                elif match.group('constrainttypes'):
                    logger.debug("Parsing [ constrainttypes ]...")
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        #newPairType = PairType(expanded.pop(i).split())
                        #self.pairtypes.add(newPairType)
                        expanded.pop(i)

                elif match.group('nonbond_params'):
                    logger.debug("Parsing [ nonbond_params ]...")
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newNonbondedType = None
                        if int(split[2]) == 1:
                            if System._sys.combination_rule == 1:
                                sigma = (float(split[4]) / float(split[3])) ** (1.0/6.0)
                                sigma *= units.kilojoules_per_mole * units.nanometers**(6)
                                epsilon = float(split[3]) / (4 * sigma**6)
                                epsilon *= units.kilojoules_per_mole * units.nanometers**(12)
                                newNonbondedType = NonbondedLJCR1Type(split[0], split[1], split[2],
                                        sigma, epsilon)
                            elif System._sys.combination_rule in (2, 3):
                                sigma = float(split[3]) * units.nanometers
                                epsilon = float(split[4]) * units.kilojoules_per_mole
                                newNonbondedType = NonbondedLJCR23Type(split[0], split[1], split[2],
                                        sigma, epsilon)

                        elif int(split[2]) == 2:
                            # TODO
                            warn("Found Buckingham entry in [ nonbond_param ]. Not yet implemented!")
                        else:
                            warn("Found unknown entry type in [ nonbond_param ]. Ignoring.")
                        System._sys._nonbonded.add(newNonbondedType)

                elif match.group('system'):
                    logger.debug("Parsing [ system ]...")
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
                    logger.debug("Parsing [ moleculetype ]...")
                    expanded.pop(i)
                    split = expanded[i].split()

                    moleculeName = split[0]
                    currentMolecule = Molecule(moleculeName)
                    System._sys.add_molecule(currentMolecule)
                    currentMoleculeType = System._sys._molecules[moleculeName]
                    currentMoleculeType.nrexcl = int(split[1])
                    expanded.pop(i)

                elif match.group('atoms'):
                    logger.debug("Parsing [ atoms ]...")
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()

                        atom = Atom(int(split[0]),          # AtomNum  (index)
                                split[4].strip(),           # name
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
                                atom.atomic_number = atomType.atomic_number
                                if not atom.bondtype:
                                    if atomType.bondtype:
                                        atom.bondtype = atomType.bondtype
                                    else:
                                        logger.warn("A suspicious parameter was found in atom/atomtypes. Visually inspect before using.\n")
                                if atom._mass[index]._value < 0:
                                    if atomType.mass._value >= 0:
                                        atom.setMass(index, atomType.mass)
                                    else:
                                        logger.warn("A suspicious parameter was found in atom/atomtypes. Visually inspect before using.\n")
                                # Assuming ptype = A
                                #atom.ptype = atomType.ptype

                                atom.setSigma(index, atomType.sigma)
                                atom.setEpsilon(index, atomType.epsilon)

                            else:
                                logger.warn("A corresponding AtomType was not found. Insert missing values yourself.\n")
                            index += 1

                        currentMolecule.addAtom(atom)
                elif match.group('bonds'):
                    logger.debug("Parsing [ bonds ]...")
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
                                split.append(bondType.k)

                            elif isinstance(bondType, CubicBondType):
                                split.append(bondType.length)
                                split.append(bondType.C2)
                                split.append(bondType.C3)

                            elif isinstance(bondType, MorseBondType):
                                split.append(bondType.length)
                                split.append(bondType.D)
                                split.append(bondType.beta)

                            elif isinstance(bondType, HarmonicPotentialType):
                                split.append(bondType.length)
                                split.append(bondType.k)
                            else:
                                warn("Bondtype '{0}' is unsupported or something more complicated went wrong".format(bondType.type))

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

                        if int(split[2]) == 6:
                            try:
                                newBondForce = HarmonicPotential(int(split[0]),
                                        int(split[1]),
                                        float(split[3]) * units.nanometers,
                                        float(split[4]) * units.kilojoules_per_mole * units.nanometers**(-2))
                            except:
                                newBondForce = HarmonicPotential(int(split[0]),
                                        int(split[1]),
                                        split[3],
                                        split[4])

                        currentMoleculeType.bondForceSet.add(newBondForce)

                elif match.group('pairs'):
                    logger.debug("Parsing [ pairs ]...")
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newPairForce = None

                        if int(split[2]) == 1:
                            if len(split) == 3:
                                # this probably won't work due to units
                                newPairForce = AbstractPair(int(split[0]), int(split[1]), "Both")

                        currentMoleculeType.pairForceSet.add(newPairForce)

                elif match.group('angles'):
                    logger.debug("Parsing [ angles ]...")
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newAngleForce = None

                        if len(split) == 4:
                            atomtype1 = currentMolecule._atoms[int(split[0])-1].bondtype
                            atomtype2 = currentMolecule._atoms[int(split[1])-1].bondtype
                            atomtype3 = currentMolecule._atoms[int(split[2])-1].bondtype
                            tempType = AbstractAngleType(atomtype1, atomtype2, atomtype3)
                            angleType = self.angletypes.get(tempType)
                            if not (angleType):
                                #flip it around.
                                tempType = AbstractAngleType(atomtype3, atomtype2, atomtype1)
                                angleType = self.angletypes.get(tempType)

                            if isinstance(angleType, AngleType):
                                split.append(angleType.theta)
                                split.append(angleType.k)

                            if isinstance(angleType, G96AngleType):
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
                                        float(split[5]) * units.kilojoules_per_mole * units.radians**(-2),
                                        float(split[6]) * units.nanometers,
                                        float(split[7]) * units.kilojoules_per_mole * units.nanometers**(-2))
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
                                newAngleForce = QuarticAngle(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        float(split[4]) * units.degrees,
                                        float(split[5]) * units.kilojoules_per_mole,
                                        float(split[6]) * units.kilojoules_per_mole * units.radians**(-1),
                                        float(split[7]) * units.kilojoules_per_mole * units.radians**(-2),
                                        float(split[8]) * units.kilojoules_per_mole * units.radians**(-3),
                                        float(split[9]) * units.kilojoules_per_mole * units.radians**(-4))
                            except:
                                newAngleForce = QuarticAngle(int(split[0]),
                                        int(split[1]),
                                        int(split[2]),
                                        split[4],
                                        split[5],
                                        split[6],
                                        split[7],
                                        split[8],
                                        split[9])

                        currentMoleculeType.angleForceSet.add(newAngleForce)

                elif match.group('dihedrals'):
                    logger.debug("Parsing [ dihedrals ]...")
                    expanded.pop(i)
                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newDihedralForce = None
                        dihedralType = None

                        dtype = int(split[4])
                        improper = (dtype == 4) or (dtype == 2)
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
                                             ['X', atomtype3, atomtype2, 'X'], # flipped double wildcard
                                             ['X', 'X', atomtype3, atomtype4], # end double wildcard
                                             [atomtype1, atomtype2,'X', 'X'] # flipped end double wildcard
                                             ]

                            for alist in atomtypelists:
                                if not dihedralType:
                                    tempType = AbstractDihedralType(alist[0], alist[1], alist[2], alist[3], improper)
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

                            elif isinstance(dihedralType, ImproperHarmonicDihedralType):
                                phi = dihedralType.xi
                                k = dihedralType.k
                            else:
                                warn("Did not find matching dihedraltype!")

                        atom1 = int(split[0])
                        atom2 = int(split[1])
                        atom3 = int(split[2])
                        atom4 = int(split[3])
                        nentries = len(split)

                        # Proper Dihedral 1
                        if int(split[4]) == 1 or int(split[4]) == 9 or int(split[4]) == 4:

                            if nentries > 5:
                                fc0, fc1, fc2, fc3, fc4, fc5, fc6 = ConvertDihedralFromProperDihedralToDihedralTrig(
                                        float(split[6]) * units.kilojoules_per_mole, int(split[7]))
                                phi = float(split[5]) * units.degrees
                            newDihedralForce = DihedralTrigDihedral(
                                atom1, atom2, atom3, atom4, phi,
                                fc0, fc1, fc2, fc3, fc4, fc5, fc6, improper=improper)

                            # for dihedral type 9, there can be multiple interactions

                        # Improper Dihedral 2
                        elif dtype == 2:
                            if nentries > 5:
                                phi = float(split[5]) * units.degrees
                                k = float(split[6]) * units.kilojoules_per_mole * units.radians**(-2)

                            newDihedralForce = ImproperHarmonicDihedral(
                                atom1, atom2, atom3, atom4, phi, k)

                        # RBDihedral
                        elif dtype == 3:

                            if nentries > 5:
                                fc0, fc1, fc2, fc3, fc4, fc5, fc6 = ConvertDihedralFromRBToDihedralTrig(
                                    float(split[5]) * units.kilojoules_per_mole,
                                    float(split[6]) * units.kilojoules_per_mole,
                                    float(split[7]) * units.kilojoules_per_mole,
                                    float(split[8]) * units.kilojoules_per_mole,
                                    float(split[9]) * units.kilojoules_per_mole,
                                    float(split[10]) * units.kilojoules_per_mole,
                                    0 * units.kilojoules_per_mole)

                            newDihedralForce = DihedralTrigDihedral(
                                atom1, atom2, atom3, atom4,
                                0 * units.degrees, fc0, fc1, fc2, fc3, fc4, fc5, fc6)  # need to look at the use of sign here

                        elif dtype == 5:

                            if nentries > 5:
                                fc0, fc1, fc2, fc3, fc4, fc5, fc6 = ConvertDihedralFromFourierToDihedralTrig(
                                    float(split[5])*units.kilojoules_per_mole,
                                    float(split[6])*units.kilojoules_per_mole,
                                    float(split[7])*units.kilojoules_per_mole,
                                    float(split[8])*units.kilojoules_per_mole,
                                    )

                            newDihedralForce = DihedralTrigDihedral(
                                atom1, atom2, atom3, atom4,
                                0 * units.degrees, fc0, fc1, fc2, fc3, fc4, fc5, fc6)  # need to look at the use of sign here

                        elif dtype == 8:
                            raise Exception("Cannot support tabulated dihedrals")
                        else:
                            raise Exception("Unsupported dihedral")

                        if dtype == 9:
                            # we need to retrive the information add to the dihedral,
                            # rather than overwrite it
                            # warning: right now, I don't think we can have BOTH duplicate dihedrals and
                            # duplicate dihedral types.  Will have to investigate.
                            try:
                                # retrive a dihedral with the same atoms, add to it if it exists
                                dihedralmatch = currentMoleculeType.dihedralForceSet.map[newDihedralForce]
                                dihedralmatch.sum_parameters(newDihedralForce)
                            except Exception as e:
                                logger.exception(e) # used to be pass, now recorded but surpressed
                        currentMoleculeType.dihedralForceSet.add(newDihedralForce)

                elif match.group('constraints'):
                    logger.debug("Parsing [ constraints ]...")
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        expanded.pop(i)

                elif match.group('settles'):
                    logger.debug("Parsing [ settles ]...")
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        newSettlesForce = None

                        if len(split) == 4:
                            newSettlesForce = Settles(int(split[0]),
                                    float(split[2]) * units.nanometers,
                                    float(split[3]) * units.nanometers)

                        currentMoleculeType.settles = newSettlesForce

                        # we need to add a constrainted bonded forces as well between the atoms in these molecules.
                        # we assume the gromacs default of 1. O, 2. H, 3. H
                        # reference bond strength is 900 kj/mol, but doesn't really matter since constrainted.
                        #waterbondrefk = 900*units.kilojoules_per_mole * units.nanometers**(-2)
                        #wateranglerefk = 400*units.kilojoules_per_mole * units.degrees**(-2)
                        
                        # From oplsaa.ff/tip3p.itp - JPT
                        waterbondrefk = 502416.0 * units.kilojoules_per_mole * units.nanometers**(-2)
                        wateranglerefk = 628.02 * units.kilojoules_per_mole * units.radians**(-2)

                        angle = 2.0 * math.asin(0.5 * float(split[3]) / float(split[2])) * units.radians
                        dOH = float(split[2]) * units.nanometers

                        newBondForce = Bond(1,2,dOH,waterbondrefk,c=True)
                        currentMoleculeType.bondForceSet.add(newBondForce)

                        newBondForce = Bond(1,3,dOH,waterbondrefk,c=True)
                        currentMoleculeType.bondForceSet.add(newBondForce)

                        newAngleForce = Angle(3,1,2,angle,wateranglerefk,c=True)
                        currentMoleculeType.angleForceSet.add(newAngleForce)

                elif match.group('exclusions'):
                    logger.debug("Parsing [ exclusions ]...")
                    expanded.pop(i)

                    while not (expanded[i].count('[')) and i < len(expanded)-1:
                        split = expanded.pop(i).split()
                        for j in range(len(split)):
                            if split[0] < split[j]:
                                newExclusion = Exclusions([int(split[0]),int(split[j])])
                                currentMoleculeType.exclusions.add(newExclusion)


                elif match.group('molecules'):
                    logger.debug("Parsing [ molecules ]...")
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
                            System._sys.add_molecule(mol)
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
                logger.debug("=========================")
                logger.debug("Original line: %s" % line)
                match = lineReg.match(line)   # ensure that the line matches what we expect!
                if match:
                    line = match.group('directive')
                    comment = match.group('comment')
                    logger.debug("Directive: %s" % line)
                    logger.debug("Comment: %s" % comment)
                    while True:     # expand the lines
                        if match.group('ignore_newline'):
                            logger.debug("Ignore newline found")
                            if lines:
                                temp = lines.popleft()
                                match = lineReg.match(temp)
                                line += match.group('directive')
                            else:
                                logger.warn("Previous line continues yet EOF encountered")
                                break   # EOF so we can't loop again
                        else:
                            break   # line terminates
                    line = line.strip()  # just in case theres whitespace we missed

                    match = preReg.match(line)
                    if match is not None:
                        if match.group('ifdef') and write[condDepth]:
                            logger.debug('Found an ifdef: %s' % match.group('ifdef'))
                            condDepth += 1
                            ifdef = match.group('ifdef')
                            if ifdef in self.defines:
                                write.append(True)
                            else:
                                write.append(False)
                        elif match.group('ifndef') and write[condDepth]:
                            logger.debug('Found an ifndef: %s' % match.group('ifndef'))
                            condDepth += 1
                            ifndef = match.group('ifndef')
                            if ifndef not in self.defines:
                                write.append(True)
                            else:
                                write.append(False)
                        elif match.group('else'):
                            logger.debug('Found an else')
                            if condDepth > 0:
                                write[condDepth] = not write[condDepth]
                            else:
                                logger.warn("Found an else not associated with a conditional")
                        elif match.group('endif'):
                            logger.debug('Found an endif')
                            if condDepth > 0:
                                condDepth -= 1
                                write.pop()
                            else:
                                raise Exception("Found an endif not associated with a conditional")
                        elif write[condDepth]:
                            if match.group('include'):
                                logger.debug("Found a include: %s" % line)
                                include = match.group('include')
                                fd = self.open(include)
                                if fd is not None:
                                    tempList = list(fd)
                                    tempList.reverse()
                                    lines.extendleft(tempList)
                            elif match.group('defineLiteral'):
                                logger.debug("Found a define: %s" % line)
                                define = match.group('defineLiteral')

                                if define not in self.defines:
                                    self.defines[define] = match.group('defineData')
                                else:
                                    logger.warn("WARNING: Overriding define: %s" % define)
                                    self.define[define] = match.group('defineData')
                            elif match.group('undef'):
                                logger.debug("Found a undefine: %s" % line)
                                undef = match.group('undef')
                                if undef in self.defines:
                                    self.defines.pop(undef)
                    elif write[condDepth]:
                        if line != '':
                            for define in self.defines:
                                if define in line:
                                    line = line.replace(define, self.defines[define])
                            logger.debug("Writing: %s" % line)
                            expanded.append(line)
                        if comment:
                            self.comments.append(comment)
                else:
                    raise Exception("Unreadable line!")
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
            logger.warn("Omitting file %s. It has already been included!" % filename)
            return None
        temp = filename

        if os.path.exists(filename):
            try:
                fd = open(filename)
                self.includes.add(temp)
                logger.info("Local instance of topology file '%s' used\n" % (filename))
                return fd
            except IOError, (errno, strerror):
                logger.exception("I/O error(%d): %s for local instance of '%s'\n" % (errno, strerror, filename))

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
                logger.info("version in %s used for topology file '%s'\n" % (parentdir, filename))
                return fd
            except Exception as e:
                logger.exception(e) # EDZ: used to be pass, now recorded but supressed
        else:
            try:
                filename = os.path.join(os.environ['GMXLIB'], temp)
                fd = open(filename)
                self.includes.add(filename)
                logger.info("version in GMXLIB = %s used for topology file '%s'\n" % (os.environ['GMXLIB'], temp))
                return fd

            except IOError, (errno, strerror):
                logger.exception("Unrecoverable I/O error(%d): can'd find %s anywhere: '%s'\n" % (errno, strerror, filename))
                sys.exit()

    def write_topology(self, filename):
        """Write this topology in GROMACS file format.

        Args:
            filename: the name of the file to write out to
        """
        lines = list()

        # [ defaults ]
        lines.append('[ defaults ]\n')
        lines.append('; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n')
        lines.append('%d%16d%18s%20.4f%8.4f\n'
                % (System._sys.nonbonded_function,
                   System._sys.combination_rule,
                   System._sys.genpairs,
                   System._sys.lj_correction,
                   System._sys.coulomb_correction))
        lines.append('\n')

        # [ atomtypes ]
        lines.append('[ atomtypes ]\n')
        lines.append(';type, bondtype, atomic_number, mass, charge, ptype, sigma, epsilon\n')
        atomtypelist = sorted(System._sys._atomtypes.itervalues(), key=lambda x: x.atomtype)
        for atomtype in atomtypelist:
#            if atomtype.atomtype.isdigit():
#                atomtype.atomtype = "LMP_" + atomtype.atomtype
#            if atomtype.bondtype.isdigit():
#                atomtype.bondtype = "LMP_" + atomtype.bondtype
            if System._sys.combination_rule == 1:
                lines.append('%-11s%5s%6d%18.8f%18.8f%5s%18.8e%18.8e\n'
                        % (atomtype.atomtype,
                           atomtype.bondtype,
                           int(atomtype.atomic_number),
                           atomtype.mass.in_units_of(units.atomic_mass_unit)._value,
                           atomtype.charge.in_units_of(units.elementary_charge)._value,
                           atomtype.ptype,
                           atomtype.sigma.in_units_of(units.kilojoules_per_mole * units.nanometers**(6))._value,
                           atomtype.epsilon.in_units_of(units.kilojoules_per_mole * units.nanometers**(12))._value))
            elif System._sys.combination_rule in (2, 3):
                lines.append('%-11s%5s%6d%18.8f%18.8f%5s%18.8e%18.8e\n'
                        % (atomtype.atomtype,
                           atomtype.bondtype,
                           int(atomtype.atomic_number),
                           atomtype.mass.in_units_of(units.atomic_mass_unit)._value,
                           atomtype.charge.in_units_of(units.elementary_charge)._value,
                           atomtype.ptype,
                           atomtype.sigma.in_units_of(units.nanometers)._value,
                           atomtype.epsilon.in_units_of(units.kilojoules_per_mole)._value))
        lines.append('\n')

        if System._sys._nonbonded:
            # [ nonbond_params ]
            lines.append('[ nonbond_params ]\n')
            nonbondedlist = sorted(moleculeType.bondForceSet.itervalues(), key=lambda x: (x.atom1, x.atom2))
            for nonbonded in nonbondedlist:
                if System._sys.combination_rule == 1:
                    lines.append('{0:6s} {1:6s} {2:3d} {3:18.8e} {4:18.8e}\n'.format(
                            nonbonded.atom1, nonbonded.atom2, nonbonded.type,
                            nonbonded.sigma.in_units_of(units.kilojoules_per_mole * units.nanometers**(6))._value,
                            nonbonded.epsilon.in_units_of(units.kilojoules_per_mole * units.nanometers**(12))._value))
                elif System._sys.combination_rule in (2, 3):
                    lines.append('{0:6s} {1:6s} {2:3s} {3:18.8e} {4:18.8e}\n'.format(
                            nonbonded.atom1, nonbonded.atom2, nonbonded.type,
                            nonbonded.sigma.in_units_of(units.nanometers)._value,
                            nonbonded.epsilon.in_units_of(units.kilojoules_per_mole)._value))
        lines.append('\n')

        # [ moleculetype]
        moleculeTypelist = sorted(System._sys._molecules.itervalues(), key=lambda x: x.name)        
        for moleculeType in moleculeTypelist:
            lines.append('[ moleculetype ]\n')
            # gromacs can't handle spaces in the molecule name
            printname = moleculeType.name
            printname = printname.replace(' ','_')
            printname = printname.replace('"','')
            lines.append('%s%10d\n\n'
                    % (printname,
                       moleculeType.nrexcl))

            # [ atoms ]
            lines.append('[ atoms ]\n')
            lines.append(';num, type, resnum, resname, atomname, cgnr, q, m\n')
            molecule = moleculeType.moleculeSet[0]
            count = 1
            for atom in molecule._atoms:
#                if atom.name.isdigit():
#                    atom.name = "LMP_" + atom.name
#                if atom._atomtype[0].isdigit():
#                    atom._atomtype[0] = "LMP_" + atom._atomtype[0]

                try:
                    lines.append('%6d%18s%6d%8s%8s%6d%18.8f%18.8f%18s%18.8f%18.8f\n'
                            % (count,
                               atom._atomtype[0],
                               atom.residue_index,
                               atom.residue_name,
                               atom.name,
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
                                    atom.residue_index,
                                    atom.residue_name,
                                    atom.name,
                                    atom.cgnr,
                                    atom._charge[0].in_units_of(units.elementary_charge)._value,
                                    atom._mass[0].in_units_of(units.atomic_mass_unit)._value))
                count += 1
            lines.append('\n')

            if moleculeType.bondForceSet and not moleculeType.settles:
                # [ bonds ]
                lines.append('[ bonds ]\n')
                lines.append(';   ai     aj funct  r               k\n')
                bondlist = sorted(moleculeType.bondForceSet.itervalues(), key=lambda x: (x.atom1,x.atom2))
                for bond in bondlist:
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
                                   bond.k.in_units_of(units.kilojoules_per_mole*units.nanometers**(-4))._value))
                    elif isinstance(bond, MorseBond):
                        b_type = 3
                        lines.append('%6d%7d%4d%5.8f%5.8f%5.8f\n'
                                % (bond.atom1,
                                   bond.atom2,
                                   b_type,
                                   bond.length.in_units_of(units.nanometers)._value,
                                   bond.D.in_units_of(units.kilojoules_per_mole)._value,
                                   bond.beta.in_units_of(units.nanometers**(-1))._value))
                    elif isinstance(bond, CubicBond):
                        b_type = 4
                        lines.append('%6d%7d%4d%5.8f%5.8f%5.8f\n'
                                % (bond.atom1,
                                   bond.atom2,
                                   b_type,
                                   bond.length.in_units_of(units.nanometers)._value,
                                   bond.C2.in_units_of(units.kilojoules_per_mole * units.nanometers**(-2))._value,
                                   bond.C3.in_units_of(units.kilojoules_per_mole * units.nanometers**(-3))._value))
                    else:
                        raise Exception("WriteError: found unsupported bond type")
                lines.append('\n')

            #MRS: why are there two pairs sections?
            if moleculeType.pairForceSet:
                # [ pair ]
                lines.append('[ pairs ]\n')
                lines.append(';  ai    aj   funct\n')
                pairlist = sorted(moleculeType.pairForceSet.itervalues(), key=lambda x: (x.atom1, x.atom2))
                for pair in pairlist:

                    if isinstance(pair, AbstractPair):
                        p_type = 1
                        lines.append('%6d%7d%4d\n'
                                % (pair.atom1,
                                   pair.atom2,
                                   p_type))

                    else:
                        raise Exception("WriteError: found unsupported pair type")
                lines.append('\n')

            if moleculeType.angleForceSet and not moleculeType.settles:
                # [ angles ]
                lines.append('[ angles ]\n')
                lines.append(';   ai     aj     ak    funct   theta         cth\n')

                anglelist = sorted(moleculeType.angleForceSet.itervalues(), key=lambda x: (x.atom1, x.atom2, x.atom3))
                for angle in anglelist:
                    atomindex = "%6d%7d%7d" % (angle.atom1,angle.atom2,angle.atom3)
                    if isinstance(angle, Angle):
                        a_type = 1
                        lines.append('%s%4d%18.8e%18.8e\n'
                                % (atomindex,
                                   a_type,
                                   angle.theta.in_units_of(units.degrees)._value,
                                   angle.k.in_units_of(units.kilojoules_per_mole*units.radians**(-2))._value))
                    elif isinstance(angle, G96Angle):
                        a_type = 2
                        lines.append('%s%4d%18.8e%18.8e\n'
                                     % (atomindex,
                                        a_type,
                                        angle.theta.in_units_of(units.degrees)._value,
                                        angle.k.in_units_of(units.kilojoules_per_mole)._value))
                    elif isinstance(angle, UreyBradleyAngle):
                        a_type = 5
                        lines.append('%s%4d%18.8e%18.8e%18.8e%18.8e\n'
                                     % (atomindex,
                                        a_type,
                                        angle.theta.in_units_of(units.degrees)._value,
                                        angle.k.in_units_of(units.kilojoules_per_mole*units.radians**(-2))._value,
                                        angle.r.in_units_of(units.nanometers)._value,
                                        angle.kUB.in_units_of(units.kilojoules_per_mole*units.nanometers**(-2))._value))
                    else:
                        raise Exception("WriteError: found unsupported angle type")
                lines.append('\n')

            """
            # [ pairs]
            lines.append('[ pairs ]\n')
            lines.append(';   ai     aj    funct')

            pairlist = sorted(moleculeType.pairForceSet.itervalues(), key=lambda x: x.atom1)
            for pair in pairlist:
                if isinstance(pair, LJ1PairCR1) or isinstance(pair, LJ1PairCR23)
                    type = 1
                    lines.append('%6d%7d%4d%18.8e%18.8e\n'%(pair.atom1, pair.atom2, type,
                            pair.V.in_units_of(units.XXX)._value,
                            pair.W.in_units_of(units.XXX)._value)

                elif isinstance(pair, LJ2PairCR1) or isinstance(pair, LJ2PairCR23):
                    type = 2
                    lines.append('%6d%7d%4d%18.8e%18.8e\n'%(pair.atom1, pair.atom2, type,
                            pair.V.in_units_of(units.XXX)._value,
                            pair.W.in_units_of(units.XXX)._value))

                elif isinstance(pair, LJNBCR1) or isinstance( pair, LJNBCR23):
                    type = 1
                    lines.append('%6d%7d%4d%18.8f%18.8f%18.8f%18.8f\n'%(pair.atom1, pair.atom2, type,
                            pair.qi.in_units_of(units.XXX)._value,
                            pair.qj.in_units_of(units.XXX)._value,
                            pair.V.in_units_of(units.XXX)._value,
                            pair.W.in_units_of(units.XXX)._value))
                else:
                    print "Could not identify pair!"
            """

            if moleculeType.dihedralForceSet:
                # [ dihedrals ]
                lines.append('[ dihedrals ]\n')
                lines.append(';    i      j      k      l   func\n')
                dihedrallist = sorted(moleculeType.dihedralForceSet.itervalues(),
                        key=lambda x: (x.atom1, x.atom2, x.atom3, x.atom4))
                for dihedral in dihedrallist:
                    # this atom index will be the same for all of types.
                    atomindex = "%7d%7d%7d%7d" % (dihedral.atom1, dihedral.atom2,
                            dihedral.atom3, dihedral.atom4)
                    if isinstance(dihedral, DihedralTrigDihedral):
                        # convienience array
                        coefficients = [dihedral.fc1, dihedral.fc2, dihedral.fc3,
                                dihedral.fc4, dihedral.fc5, dihedral.fc6]
                        if dihedral.improper:
                            found_nonzero = False
                            for n, coeff in enumerate(coefficients):  # only one of these should be nonzero
                                if coeff._value != 0.0:
                                    if found_nonzero == False:
                                        found_nonzero = True
                                    else:
                                        raise ValueError("Found more than one nonzero "
                                                "coefficient in improper trigonal dihedral!")
                                    lines.append('%s%4d%18.8f%18.8f%6d\n'
                                                 % (atomindex, 4,
                                                    dihedral.phi.in_units_of(units.degrees)._value,
                                                    coeff.in_units_of(units.kilojoules_per_mole)._value,
                                                    n + 1))
                        else:
                            rb_coeffs = ConvertDihedralFromDihedralTrigToRB(
                                np.cos(dihedral.phi.in_units_of(units.radians)._value),
                                dihedral.phi, dihedral.fc0, * coefficients)
                            # there are some cases where some dihedrals will have c[6] values, which gromacs
                            # can't yet handle.  We need a workaround, and will go through the route of multiple 9s.
                            if (dihedral.phi in [0*units.degrees, 180*units.degrees] and
                                    rb_coeffs[6]._value == 0):
                                d_type = 3
                                lines.append('%s%4d%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f\n'
                                             % (atomindex,
                                                d_type,
                                                rb_coeffs[0].in_units_of(units.kilojoules_per_mole)._value,
                                                rb_coeffs[1].in_units_of(units.kilojoules_per_mole)._value,
                                                rb_coeffs[2].in_units_of(units.kilojoules_per_mole)._value,
                                                rb_coeffs[3].in_units_of(units.kilojoules_per_mole)._value,
                                                rb_coeffs[4].in_units_of(units.kilojoules_per_mole)._value,
                                                rb_coeffs[5].in_units_of(units.kilojoules_per_mole)._value))
                            else:
                                ncount = sum(coeff._value != 0.0 for coeff in coefficients)
                                if ncount > 1:
                                    dtype = 9
                                else:
                                    dtype = 1
                                # all of the terms should have consistent phi now
                                for n, coeff in enumerate(coefficients):
                                    if coeff._value < 0:
                                        # kludge for different definitions of trigonometric
                                        # dihedral with multiplicity in the central representation
                                        # and in GROMACS
                                        coefficients[n] *= -1
                                        if dihedral.phi.value_in_unit(units.degrees)==180:
                                            printphi = 0*units.degrees
                                        else:
                                            printphi = 180*units.degrees
                                    else:
                                        printphi = dihedral.phi
                                    lines.append('%s%4d%18.8f%18.8f%6d\n'
                                                 % (atomindex,
                                                    dtype,
                                                    printphi.in_units_of(units.degrees)._value,
                                                    coeff.in_units_of(units.kilojoules_per_mole)._value,
                                                    n + 1))

                    elif isinstance(dihedral, ImproperHarmonicDihedral):
                        d_type = 2
                        lines.append('%s%4d%18.8f%18.8f\n'
                                     % (atomindex,
                                        d_type,
                                        dihedral.xi.in_units_of(units.degrees)._value,
                                        dihedral.k.in_units_of(units.kilojoules_per_mole*units.radians**(-2))._value))

                    else:
                        raise Exception("WriteError: found unsupported dihedral type")
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
                exclusionlist = sorted(moleculeType.exclusions.itervalues(), key=lambda x: (x.exclusions[0], x.exclusions[1]))
                for exclusion in exclusionlist:
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
            printname = molType
            printname = printname.replace(' ','_')
            printname = printname.replace('"','')            
            lines.append('%-15s%8d\n'
                    % (printname,
                      len(System._sys._molecules[molType].moleculeSet)))

        fout = open(filename, 'w')
        for line in lines:
            fout.write(line)
        fout.close()

