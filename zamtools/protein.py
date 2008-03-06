#!/usr/bin/env python

#LAST MODIFIED: 04-13-07

from numpy import *
from proteinconst import *
import proteinfunc as pfunc

import os, sys, copy, gzip, time
import sequence, coords
from geometry import *


#try to load the compiled functions
USELIB = True
try:
  import proteinlib
except ImportError:
  USELIB = False


#GLOBALS

#pdb file with amino acids
TemplateFile = "template.pdb"

#defaults
ResAtomDflt = "*"         #default atom for residue-residue distances (* is centroid)
ResRadiusDflt = 8.0       #default residue-residue cutoff distance
AtomRadiusDflt = 4.0      #default atom-atom cutoff distance
SaltAtomRadiusDflt = 4.0  #default atom-atom cutoff distance for salt atoms
IntMinCoordDflt = 9       #minimum coordination number of interior residues
                          #9 is for use with an 8A cutoff; 7 for a 7A cutoff
MinCODflt = 3             #minimum contact order for contacts calculation
AssignBondsDflt = True    #whether or not to automatically assign bonds to pdbfiles without explicit bonds


     
#======== CLASSES =========

class AtomClass:
  def __init__(self, Element = "", Name = "", BFactor = 0.,
    Mass = 0., FormalCharge = None, Occupancy = 1., ResName = ""):
    self.Element = Element
    self.FullName = AtomConvertAliases.get(Name, Name)
    self.Name = AtomConvertAliases.get(Name.strip(), Name.strip())
    #move numbers to the end of the name
    if self.Name[0] in "0123456789": self.Name = self.Name[1:] + self.Name[0]
    self.BFactor = BFactor
    self.Mass = Mass
    self.FormalCharge = FormalCharge
    self.Occupancy = Occupancy
    #get force field parameters
    self.Radius, self.SqrtEps, self.Charge = pfunc.GetFFData(self.Element, self.Name, ResName)
  def __repr__(self):
    return self.Name

class ResClass:
  def __init__(self, Name = "", Chain = "", Atoms = None,
               StartAtom = 0, Bonds = None):
    self.FullName = Name
    self.Name = Name.strip()
    self.Chain = Chain
    self.Atoms = Atoms
    self.Bonds = Bonds
    if self.Atoms is None: self.Atoms = []
    if self.Bonds is None: self.Bonds = []
    #StartAtom and StopAtom are the only members of ResClass
    #that are updated by ProteinClass and not ResClass routines
    self.StartAtom = StartAtom
    self.StopAtom = self.StartAtom + len(self.Atoms)
    self.ChiAtoms = []
    if len(self.Bonds) > 0: self.__GenerateChi()
  def __repr__(self):
    return self.Name
  def __len__(self):
    return len(self.Atoms)
  def HasAtom(self, AtomName, Aliases = []):
    if type(Aliases) == dict:
      Aliases = Aliases.get(AtomName, Aliases["*"])
      Aliases = Aliases.get(self.Name, Aliases["*"])
    AtomNames = [a.Name for a in self.Atoms]
    for Alias in [AtomName] + Aliases:
      if Alias in AtomNames: return True
    return False
  def AtomNum(self, AtomName, Aliases = [],
              NotFoundError = True, Absolute = False):
    if type(Aliases) == dict:
      Aliases = Aliases.get(AtomName, Aliases["*"])
      Aliases = Aliases.get(self.Name, Aliases["*"])
    AtomNames = [a.Name for a in self.Atoms]
    for Alias in [AtomName] + Aliases:
      if Alias in AtomNames:
        if Absolute:
          return AtomNames.index(Alias) + self.StartAtom
        else:
          return AtomNames.index(Alias)
    if NotFoundError:
      raise IndexError, "Atom name %s not found in res %s." % (AtomName, self.Name)
    else:
      return -1
  def AddAtom(self, Atom, Ind = None):
    if Ind is None:
      self.Atoms.append(Atom)
    else:
      self.Atoms = self.Atoms[:Ind] + [Atom] + self.Atoms[Ind:]
  def DelAtom(self, AtomNum):
    del self.Atoms[AtomNum]
    NewBonds = []
    for (a,b) in self.Bonds:
      if a == AtomNum or b == AtomNum: continue
      if a > AtomNum: a -= 1
      if b > AtomNum: b -= 1
      NewBonds.append((a,b))
    self.Bonds = NewBonds
    self.__GenerateChi()
  def AddBond(self, Bond, UpdateChi = True):
    a, b = min(Bond), max(Bond)
    if not (a,b) in self.Bonds: self.Bonds.append((a,b))
    if UpdateChi: self.__GenerateChi()
  def DelBond(self, Bond, UpdateChi = True):
    a, b = min(Bond), max(Bond)
    self.Bonds = [x for x in self.Bonds if not x == (a,b)]
    if UpdateChi: self.__GenerateChi()
  def ClearBonds(self):
    self.Bonds = []
    self.ChiAtoms = []
  def GetBonded(self, AtomNum, Depth = 1, Exclude = []):
    l = []
    if Depth < 1: return l
    for (a,b) in self.Bonds:
      if a == AtomNum and not b in Exclude:
        l.append(b)
      elif b == AtomNum and not a in Exclude:
        l.append(a)
    for a in copy.deepcopy(l):
      m = self.GetBonded(a, Depth - 1, Exclude + [AtomNum] + l)
      l = l + [x for x in m if not x in l and not x == AtomNum]
    return l
  def GetLoops(self, Branch = [-1,0]):
    Loops = []
    for x in self.GetBonded(Branch[-1], Exclude = [Branch[-2]]):
      if x in Branch:
        Loop = sorted(Branch[Branch.index(x):])
        if not Loop in Loops: Loops.append(Loop)
      else:
        for Loop in self.GetLoops(Branch + [x]):
          if not Loop in Loops: Loops.append(Loop)
    return Loops
  def GetDihedrals(self, Exclude = []):
    Dihedrals = []
    Bonds = [(a,b) for (a,b) in self.Bonds
             if not a in Exclude and not b in Exclude]
    for (a,b) in Bonds:
      aBonds = self.GetBonded(a, Exclude = Exclude + [b])
      if len(aBonds) == 0: continue
      bBonds = self.GetBonded(b, Exclude = Exclude + [a])
      if len(bBonds) == 0: continue
      Dih = [min(aBonds)] + [a,b] + [min(bBonds)]
      if Dih[-1] < Dih[0]: Dih.reverse()
      Dihedrals.append(Dih)
    return Dihedrals
  def __GenerateChi(self):
    #for this to work properly, atoms farther from CA must have
    #higher indices than bonded atoms that are closer
    def InLoop(a, b, Loops):
      return any([a in l and b in l for l in Loops])
    Loops = self.GetLoops()
    Exclude = [self.Atoms.index(a) for a in self.Atoms
               if a.Name in ["O","C"] or a.Element == "H"]
    Dihedrals = self.GetDihedrals(Exclude = Exclude)
    Dihedrals = [d for d in Dihedrals if not InLoop(d[1], d[2], Loops)]
    self.ChiAtoms = []
    for l in Dihedrals:
      Bonded = sorted(self.GetBonded(l[2], len(self.Atoms), l[0:2]))
      l.extend(Bonded)
      self.ChiAtoms.append(array(l,int))
  def GenerateBonds(self, Pos):
    Radii = [CovalentRadii.get(a.Element, CovalentRadii["default"])
             for a in self.Atoms]
    for i in range(0, len(self.Atoms)):
      for j in range(i+1, len(self.Atoms)):
        Dist = Length(Pos[i] - Pos[j])
        if Dist < Radii[i] + Radii[j]: self.AddBond((i,j), UpdateChi = False)
    self.__GenerateChi()
   

class ProteinClass:

  #list of settings to be copied from proteinclass to another
  SettingsList = ["ResAtom", "ResRadius", "AtomRadius", "SaltAtomRadius", "MinCO"]

  #items to exclude when copying or pickling
  GetStateExcludeVars = ["Trj", "CoordsObj"]
  
  def __init__(self, Seq = None, Pdb = None, Res = None, Pos = None):
    #default variables
    self.Res = []
    self.Pos = zeros((0,3), float)
    self.Chains = []
    self.ResAtom = ResAtomDflt
    self.ResRadius = ResRadiusDflt
    self.AtomRadius = AtomRadiusDflt
    self.SaltAtomRadius = SaltAtomRadiusDflt
    self.MinCO = MinCODflt
    self.StartChain = ""
    self.UseBonds = False
    #figure out what to do
    #check for filename in Seq
    if not Seq is None:
      if "/" in Seq or "\\" in Seq or "." in Seq:
        Pdb = Seq
        Seq = None
    if not Pdb is None:
      self.ReadPdb(Pdb)
    elif not Seq is None:
      self.Res = []
      self.Pos = zeros((0,3), float)
      for r in sequence.SeqToList(Seq):
        if r in Templates:
          self.Join(Templates[r])
        else:
          print "Could not find residue %s." % r
    else:
      if not Res is None: self.Res = Res
      if not Pos is None: self.Pos = Pos
    self.Update()

  def __setattr__(self, key, val):
    """Sets a class attribute."""
    #check if we have to recalculate the indices
    if key == "ResAtom":
      if not "ResAtom" in self.__dict__:
        self.__dict__[key] = val
      elif not self.__dict__[key] == val:
        self.__dict__[key] = val
        self.__UpdateResAtoms()
    elif key == "StartChain":
      if not "StartChain" in self.__dict__:
        self.__dict__[key] = val
      elif val == "" or val in ChainNames:
        self.__dict__[key] = val
        self.__UpdateChains()
    else:
      self.__dict__[key] = val

  def __repr__(self):
    return " ".join([str(x) for x in self.Res])

  def __len__(self):
    return len(self.Res)

  def __getitem__(self, Ind):
    "Returns a new ProteinClass that is a slice of the current residues."
    if type(Ind) is int:
      start, stop, step = Ind, Ind+1, 1
    elif type(Ind) is slice:
      start, stop = Ind.start, Ind.stop
    else:
      raise IndexError, "Syntax error in slice index."
      return
    if start < 0: start += len(self.Res)
    start = max(0, start)
    if stop < 0: stop += len(self.Res)
    stop = min(len(self.Res), stop)
    NewRes = copy.deepcopy(self.Res[start:stop])
    if stop == start:
      NewPos = zeros((0,3), float)
    else:
      StartAtom = self.Res[start].StartAtom
      StopAtom = self.Res[stop-1].StopAtom
      NewPos = self.Pos[StartAtom:StopAtom,:]
    p = ProteinClass()
    p.Res = NewRes
    p.Pos = NewPos
    p.GetSettings(self)
    p.Update()
    return p

  def __getstate__(self):
    """Returns a dictionary for pickling."""
    NewDict = self.__dict__.copy()
    for v in ProteinClass.GetStateExcludeVars:
      if v in NewDict: del NewDict[v]
    return NewDict  

  def __add__(self, Other):
    "Joins two proteins."
    return self.Copy().Join(Other)

  def __sub__(self, Other):
    "Joins two proteins and optimizes for non-overlap."
    n1 = len(self.Res)
    n2 = len(Other.Res)
    p = self.Copy().Join(Other)
    if n1 > 0 and n2 > 0:
      p.OptimizePep(n1)
    return p

  def __mul__(self, NCopy):
    "Duplicates a protein a specified number of times."
    p = self.Copy()
    for i in range(1, NCopy):
      p += self
    return p

  def __div__(self, NCopy):
    "Duplicates a protein a specified number of times and optimizes for non-overlap."
    p = self.Copy()
    for i in range(1, NCopy):
      p -= self
    return p
    
  def __UpdateResAtoms(self, ResNum = None):
    #find the residue atoms and maintain in an array for speed
    if self.ResAtom == "*":
      self.__ResAtomNum = zeros(len(self.Res), int)
    else:
      if ResNum is None:
        self.__ResAtomNum = [self.AtomNum(i, self.ResAtom, AtomAliases)
                             for (i,r) in enumerate(self.Res)]
        self.__ResAtomNum = array(self.__ResAtomNum, int)
      else:
        self.__ResAtomNum[ResNum] = self.AtomNum(ResNum, self.ResAtom,
                                                 AtomAliases) 

  def __UpdateQuickAtoms(self, ResNum = None, NextStartAtom = None):
    #put certain atoms in an array for speed;
    #nonexistent atoms are given the atom number -1.
    if ResNum is None:
      self.__AtomNums = {}
      for AtomName in QuickAtoms:
        self.__AtomNums[AtomName] = array([r.AtomNum(AtomName,
            NotFoundError = False, Absolute = True) for r in self.Res], int)
        self.__AtomNums[AtomName].flags["WRITEABLE"] = False
    else:
      #calculate the difference in atoms
      r = self.Res[ResNum]
      Diff = r.StartAtom + len(r.Atoms) - NextStartAtom
      for AtomName in self.__AtomNums.iterkeys():
        self.__AtomNums[AtomName].flags["WRITEABLE"] = True
        #update this residue
        self.__AtomNums[AtomName][ResNum] = self.Res[ResNum].AtomNum(
            AtomName, NotFoundError = False, Absolute = True)
        #update following residues; make sure to keep -1 values the same
        Mask = self.__AtomNums[AtomName] < 0
        self.__AtomNums[AtomName][ResNum+1:] += Diff
        self.__AtomNums[AtomName][Mask] = -1
        self.__AtomNums[AtomName].flags["WRITEABLE"] = False

  def __UpdateFFData(self, ResNum = None, NextStartAtom = None):
    #update an array of VDW radii
    if ResNum is None:
      self.Radius = array([a.Radius for a in self.Atoms], float)
      self.SqrtEps = array([a.SqrtEps for a in self.Atoms], float)
      self.Charge = array([a.Charge for a in self.Atoms], float)
    else:
      r = self.Res[ResNum]
      self.Radius = concatenate((self.Radius[:r.StartAtom],
                                 array([a.Radius for a in r.Atoms], float),
                                 self.Radius[NextStartAtom:]))
      self.SqrtEps = concatenate((self.SqrtEps[:r.StartAtom],
                                  array([a.SqrtEps for a in r.Atoms], float),
                                  self.SqrtEps[NextStartAtom:]))
      self.Charge = concatenate((self.Charge[:r.StartAtom],
                                 array([a.Charge for a in r.Atoms], float),
                                 self.Charge[NextStartAtom:]))
    #update the cutoffs (corresponding to 2.5sigma for lj)
    MaxSigma = 0.
    if len(self.Radius) > 0: MaxSigma = 2. * self.Radius.max()
    self.StericDistCut = 1.5 * MaxSigma
    self.LJDistCut = 2.2272 * MaxSigma

  def __UpdateChains(self):
    #updates the chain information
    NRes, NAtom = len(self.Res), len(self.Atoms)
    self.ChainResNums = [i for i in range(NRes) if i==0 or
                         self.Res[i].Chain != self.Res[i-1].Chain] + [NRes]
    self.ChainAtomNums = [self.Res[i].StartAtom for i in self.ChainResNums[:-1]] + [NAtom]
    #rename the chains
    NChain = len(self.ChainResNums) - 1
    if NChain < 2:
      for r in self.Res: r.Chain = self.StartChain
    else:
      s = copy.deepcopy(ChainNames)
      if not self.StartChain == "": s = s[s.index(self.StartChain):]
      for (rn, r) in enumerate(self.Res):
        if rn in self.ChainResNums: Chain = s.pop(0)
        r.Chain = Chain
    #update thelist of chains
    self.Chains = [self.Res[i].Chain for i in self.ChainResNums[:-1]]

  def __UpdateBase(self, ResNum = None, NextStartAtom = None):
    #updates all the internal numbering items and lists:
    #self.Atoms, self.AtomResNum, self.Res[:].StartAtom, self.Seq
    if ResNum is None:
      #rebuild AtomResNum and Atoms
      self.AtomResNum = []
      self.Atoms = []
      for i in range(0, len(self.Res)):
        self.AtomResNum.extend([i] * len(self.Res[i].Atoms))
        self.Atoms.extend(self.Res[i].Atoms)
      self.AtomResNum = array(self.AtomResNum, int)
      #rebuild the StartRes 
      Ind = 0
      for r in self.Res:
        r.StartAtom = Ind
        Ind += len(r.Atoms)
        r.StopAtom = Ind
      #rebuild the sequence
      self.Seq = [r.Name.strip() for r in self.Res]
    else:
      r = self.Res[ResNum]
      Ind = r.StartAtom
      n = len(r.Atoms)
      #update AtomResNum and Atoms
      self.AtomResNum = concatenate((self.AtomResNum[:Ind], [ResNum]*n,
                                     self.AtomResNum[NextStartAtom:]))
      self.Atoms = self.Atoms[:Ind] + r.Atoms + self.Atoms[NextStartAtom:]
      #update the startatoms
      for r in self.Res[ResNum:]:
        r.StartAtom = Ind
        Ind += len(r.Atoms)
        r.StopAtom = Ind
      #update the seq
      self.Seq[ResNum] = self.Res[ResNum].Name
    #reset the bond holders;
    #this will be automatically filled in when the scoring function
    #is calculated
    self.__Bonds23, self.__Bonds4 = None, None

  def Update(self, ResNum = None):
    """Updates internal class properties; not to be invoked by user."""
    #get the next startatom if we are updating just one residue;
    #this is the number corresponding to the first atom in the
    #next residue prior to the most recent change
    if ResNum is None:
      NextStartAtom = None
    else:
      NextStartAtom = self.Res[ResNum].StopAtom
    #first fix all internal lists;
    #this must be done first!!
    self.__UpdateBase(ResNum, NextStartAtom)
    #other updates
    self.__UpdateChains()
    self.__UpdateResAtoms(ResNum)
    self.__UpdateQuickAtoms(ResNum, NextStartAtom)
    self.__UpdateFFData(ResNum, NextStartAtom)

  def GetSettings(self, p):
    "Copies settings from another ProteinClass."
    for itm in ProteinClass.SettingsList:
      self.__dict__[itm] = p.__dict__[itm]

  def Copy(self):
    "Returns a deep copy of self."
    return copy.deepcopy(self)


#======== BASIC FUNCTIONS AND INPUT/OUTPUT ======== 

  def Seq(self):
    "Returns a list of 3-letter codes."
    return copy.deepcopy(self.Seq)

  def AtomNames(self):
    "Returns a list of atom names."
    return [a.Name for a in self.Atoms]

  def AtomRes(self):
    "Returns a list of atom residues."
    return [self.Res[i].Name for i in r.AtomResNum]
  
  def AddRes(self, Res, Pos = None):
    "Adds a residue class to the chain."
    self.Res.append(Res)
    if Pos is None: Pos = zeros((len(Res.Atoms), 3), float)
    self.Pos = concatenate((self.Pos, Pos))
    self.Update()

  def DelRes(self, ResNum):
    "Deletes a residue from the chain."
    r = self.Res[ResNum]
    self.Res = self.Res[:ResNum] + self.Res[ResNum+1:]
    TakeAtoms = range(0, r.StartAtom) + range(r.StopAtom, len(self.Pos))
    self.Pos = self.Pos.take(TakeAtoms, axis=0)
    self.Update()

  def AddAtom(self, ResNum, Atom, Pos = None, Ind = None):
    "Adds an atom class to the chain."
    if Pos is None:
      Pos = array((1,3), float)
    elif not Pos.shape == (1,3):
      Pos = Pos.reshape((1,3))
    if Ind is None: Ind = len(self.Res[ResNum].Atoms)
    a = self.Res[ResNum].StartAtom + Ind
    self.Res[ResNum].AddAtom(Atom, Ind)
    self.Pos = concatenate((self.Pos[:a], Pos, self.Pos[a:]))
    self.Update(ResNum = ResNum)

  def DelResAtom(self, ResNum, ResAtomNum):
    "Removes an atom from the chain, specified by residue."
    if ResAtomNum >= len(self.Res[ResNum].Atoms):
      raise IndexError, "ResAtomNum %d exceeds atom list in residue." % ResAtomNum
    a = self.Res[ResNum].StartAtom + ResAtomNum
    self.Res[ResNum].DelAtom(ResAtomNum)
    self.Pos = concatenate((self.Pos[:a], self.Pos[a+1:]))
    self.Update(ResNum = ResNum)

  def DelAtom(self, AtomNum):
    "Removes an atom specified by absolute number along the chain."
    ResNum = self.AtomResNum[AtomNum]
    ResAtomNum = AtomNum - self.Res[ResNum].StartAtom
    self.Res[ResNum].DelAtom(ResAtomNum)
    self.Pos = concatenate((self.Pos[:AtomNum], self.Pos[AtomNum+1:]))
    self.Update(ResNum = ResNum)    

  def HasAtom(self, ResNum, AtomName, Aliases = []):
    "Indicates if a residue has a given atom."
    if ResNum >=0 and ResNum < len(self.Res):
      if AtomName in self.__AtomNums:
        if self.__AtomNums[AtomName][ResNum] >= 0:
          return True
        else:
          return self.Res[ResNum].HasAtom(AtomName, Aliases)
      else:
        return self.Res[ResNum].HasAtom(AtomName, Aliases)
    else:
      raise IndexError, "Residue number %d out of range." % ResNum  

  def AtomNum(self, ResNum, AtomName, Aliases = [], NotFoundError = True):
    "Gives the number of the atom."
    if ResNum >=0 and ResNum < len(self.Res):
      if AtomName in self.__AtomNums:
        n = self.__AtomNums[AtomName][ResNum]
        if n < 0:
          m = self.Res[ResNum].AtomNum(AtomName, Aliases, NotFoundError)
          if m < 0:
            return m
          else:
            return m + self.Res[ResNum].StartAtom
        else:
          return n
      else:
        m = self.Res[ResNum].AtomNum(AtomName, Aliases, NotFoundError)
        if m < 0:
          return m
        else:
          return m + self.Res[ResNum].StartAtom
    else:
      if NotFoundError:
        raise IndexError, "Residue number %d out of range." % ResNum
      else:
        return -1

  def AtomList(self, AtomName):
    """Gives the number of the first AtomName in each residue,
    or -1 if not found.  This method is much faster than AtomInd
    if QuickAtoms are set."""
    if AtomName in self.__AtomNums:
      return self.__AtomNums[AtomName]
    else:
      return array([r.AtomNum(AtomName, NotFoundError = False, Absolute = True)
                    for r in self.Res], int)

  def ResNum(self, AtomNum):
    "Gives the residue number corresponding to an atom number."
    return self.AtomResNum[AtomNum]

  def Bonds(self):
    "Gives all bonds."
    l = []
    for (i,r) in enumerate(self.Res):
      l.extend([(a + r.StartAtom, b + r.StartAtom) for (a,b) in r.Bonds])
      if not i in self.ChainResNums:
        a = self.Res[i-1].AtomNum("C", NotFoundError = False, Absolute = True)
        b = self.Res[i].AtomNum("N", NotFoundError = False, Absolute = True)
        if a >= 0 and b >= 0: l.append((a,b))
    l.sort()
    return l

  def ResBonds(self, ResNum, IncludePep = True):
    "Gives all intraresidue bonds and the peptide links for a res number."
    r = self.Res[ResNum]
    Bonds = [(a + r.StartAtom, b + r.StartAtom) for (a,b) in r.Bonds]
    if ResNum > 0 and IncludePep:
      if self.ResChain(ResNum) == self.ResChain(ResNum-1):
        r2 = self.Res[ResNum - 1]
        a = r2.AtomNum("C", NotFoundError = False, Absolute = True)
        b = r.AtomNum("N", NotFoundError = False, Absolute = True)
        if a >= 0 and b >= 0: Bonds.append((a,b))
    if ResNum < len(self.Res) - 1 and IncludePep:
      if self.ResChain(ResNum) == self.ResChain(ResNum+1):
        r2 = self.Res[ResNum + 1]
        a = r.AtomNum("C", NotFoundError = False, Absolute = True)
        b = r2.AtomNum("N", NotFoundError = False, Absolute = True)
        if a >= 0 and b >= 0: Bonds.append((a,b))
    return Bonds

  def AtomBonds(self, AtomNum, Depth = 1, Exclude = [],
                Bonds = None):
    """"Returns all intraresidue bonds from AtomNum to other atoms.
Does not search across the peptide bond."""
    l = []
    if Depth < 1: return l
    if Bonds is None: Bonds = self.ResBonds(AtomResNum[AtomNum])
    for (a,b) in Bonds:
      if a == AtomNum and not b in Exclude:
        l.append(b)
      elif b == AtomNum and not a in Exclude:
        l.append(a)
    for a in copy.deepcopy(l):
      m = self.AtomBonds(a, Depth - 1, [AtomNum] + l, Bonds)
      l = l + [x for x in m if not x in l and not x == AtomNum]
    return l  

  def ClearBonds(self, ResInd = None):
    "Clears out all bond information."
    if ResInd is None: ResInd = range(len(self.Res))
    for i in ResInd: self.Res[i].ClearBonds()
    
  def GenerateBonds(self, ResInd = None):
    """Generates intraresidue bonds based on covalent radii."""
    if ResInd is None: ResInd = range(len(self.Res))
    for i in ResInd:
      r = self.Res[i]
      r.GenerateBonds(self.Pos[r.StartAtom:r.StopAtom])

  def NewChain(self, ResInd = None):
    """Makes a new chain for a set of residues."""
    if ResInd is None: ResInd = range(len(self.Res))
    #use a new character that will be updated
    for i in ResInd: self.Res[i].Chain = "!"
    #update the chain residues
    self.__UpdateChains()

  def ResChain(self, ResNum):
    """Returns the chain number for a residue."""
    for i in range(1, len(self.ChainResNums)):
      if ResNum < self.ChainResNums[i]:
        return i-1
    return 0

  def ChainRange(self, ResNum):
    "Gives the residue range and atom range for this chain, (r1,r2),(a1,a2)."
    i = self.ResChain(ResNum)
    return [self.ChainResNums[i:i+2], self.ChainAtomNums[i:i+2]]

  def GenerateChains(self):
    """Generates chains based on peptide bonds detected."""
    n = len(self.Res)
    PepLen = CovalentRadii["C"] + CovalentRadii["N"]
    #clear old chains
    self.NewChain()
    #detect peptide bonds
    for i in range(n-1):
      a = self.AtomNum(i, "C", NotFoundError = False)
      b = self.AtomNum(i+1, "N", NotFoundError = False)
      if a < 0 or b < 0:
        self.NewChain(ResInd = range(i+1,n))
      elif Length(self.Pos[a] - self.Pos[b]) > PepLen:
        self.NewChain(ResInd = range(i+1,n))


#======== MASKS ========

  def AtomInd(self, AtomName = None, ResName = None, ResNum = None,
             AtomCharged = None, ResCharged = None,
             ChargeSign = None, Hydrophobic = None, ChainNum = None,
             UserFunc = None, ReturnMask = False, OnePerRes = False):
    """Returns a list of array indices defaulting to all atoms,
    to which the specified property filters are applied. Names are a string
    or a list of strings. UserFunc is a user-defined function taking string
    variables ResName and AtomName, and returning True or False or None."""
    Mask = ones(len(self.Pos), bool)
    #check for names
    if not AtomName is None:
      if type(AtomName) is str: AtomName = [AtomName]
      #remove spaces
      AtomName = [x.strip() for x in AtomName]
      #check for wildcard
      if "*" in AtomName or "" in AtomName: AtomName = None
    if not ResName is None:
      if type(ResName) is str: ResName = [ResName]
      #remove spaces
      ResName = [x.strip() for x in ResName]
      #check for wildcard
      if "*" in ResName or "" in ResName: ResName = None
    if not ResNum is None:
      if type(ResNum) is int: ResNum = [ResNum]
    if not ChainNum is None:
      if type(ChainNum) is int: ChainNum = [ChainNum]
    #loop over atoms
    i = 0
    rNum = 0
    for r in self.Res:
      ThisResDone = False
      rInst = sequence.AAInst(r.Name)
      rCharge = sequence.Charge(r.Name)
      rCharged = sequence.Charged(r.Name)
      rPhobic = sequence.Hydrophobic(r.Name)
      rChainNum = self.Chains.index(r.Chain)
      for a in r.Atoms:
        aName = a.Name.strip()
        #names
        if not AtomName is None:
          Mask[i] = Mask[i] and (aName in AtomName)
        if not ResName is None:
          Mask[i] = Mask[i] and (r.Name in ResName)
        if not ResNum is None:
          Mask[i] = Mask[i] and (rNum in ResNum)
        #charged, by residue
        if not ResCharged is None:
          Mask[i] = Mask[i] and (rCharged == ResCharged)
        if not ChargeSign is None:
          Mask[i] = Mask[i] and (rCharge * ChargeSign > 0)
        #charged, by atom
        if a.Charge is None:
          if rInst is None:
            aCharge = 0
          else:
            if aName in rInst.SaltAtoms:
              aCharge = rInst.Charge
            else:
              aCharge = 0
        else:
          aCharge = a.Charge
        if not AtomCharged is None:
          Mask[i] = Mask[i] and ((aCharge != 0) == AtomCharged)
        if not ChargeSign is None:
          Mask[i] = Mask[i] and (aCharge * ChargeSign > 0)
        #hydrophobicity
        if not Hydrophobic is None:
          Mask[i] = Mask[i] and (rPhobic == Hydrophobic)
        #chain number
        if not ChainNum is None:
          Mask[i] = Mask[i] and (rChainNum in ChainNum)
        #user function
        if not UserFunc is None:
          Val = UserFunc(r.Name, aName)
          if not Val is None:
            Mask[i] = Mask[i] and Val
        i += 1
      rNum += 1
    if ReturnMask:
      return Mask
    else:
      return nonzero(Mask)[0]
            

  def ResInd(self, ResName = None, ResCharged = None, ChargeSign = None,
             Hydrophobic = None, MinCoord = None, MaxCoord = None,
             ChainNum = None, UserFunc = None, ReturnMask = False):
    """Returns a list of array indices defaulting to all residues,
    to which the specified property filters are applied. ResName is a string
    or a list of strings. UserFunc is a user-defined function taking string
    variable ResName, and returning True or False or None.
    MinCoord and MaxCoord specify limits on the coordination number."""
    Mask = ones(len(self.Res), bool)
    #check for names
    if not ResName is None:
      if type(ResName) is str: ResName = [ResName]
      #remove spaces
      ResName = [x.strip() for x in ResName]
      #check for wildcard
      if "*" in ResName or "" in ResName: ResName = None
    if not ChainNum is None:
      if type(ChainNum) is int: ChainNum = [ChainNum]
    #loop over atoms
    i = 0
    #check if any coordination numbers are specified
    if not MinCoord is None or not MaxCoord is None:
      ResCoord = self.ResCoordination()
    for r in self.Res:
      rCharge = sequence.Charge(r.Name)
      rCharged = sequence.Charged(r.Name)
      rPhobic = sequence.Hydrophobic(r.Name)
      rChainNum = self.Chains.index(r.Chain)
      #name
      if not ResName is None:
        Mask[i] = Mask[i] and (r.Name in ResName)
      #charged, by residue
      if not ResCharged is None:
        Mask[i] = Mask[i] and (rCharged == ResCharged)
      if not ChargeSign is None:
        Mask[i] = Mask[i] and (rCharge * ChargeSign > 0)
      #hydrophobicity
      if not Hydrophobic is None:
        Mask[i] = Mask[i] and (rPhobic == Hydrophobic)
      #coordination
      if not MinCoord is None:
        Mask[i] = Mask[i] and (ResCoord[i] >= MinCoord)
      if not MaxCoord is None:
        Mask[i] = Mask[i] and (ResCoord[i] <= MaxCoord)
      #chain number
      if not ChainNum is None:
        Mask[i] = Mask[i] and (rChainNum in ChainNum)
      #user function
      if not UserFunc is None:
        Val = UserFunc(r.Name)
        if not Val is None:
          Mask[i] = Mask[i] and Val
      i += 1
    if ReturnMask:
      return Mask
    else:
      return nonzero(Mask)[0]
    

#======== MONTE CARLO BASE ========

  def GetMCMoves(self, ResIndPhiPsi = None, ResIndChi = None, ChainInd = None, 
                 ResProbPhiPsi = None, ResProbChi = None, ChainProb = None,
                 WeightPhiPsi = 1.0, WeightChi = 1.0, WeightChain = 0.,
                 CenterChain = 0, Point = None):
    """Returns standard MC moves in a list of MCMoveClass objects.
ResInd is the list of residues to apply moves to; defaults to all.
ChainInd is the list of chains to apply moves to; defaults to all but CenterChain.
ResProb are probabilities for each residue; overrides ResInd.
ChainProb are probabilities for each chain; overrides ChainInd.
Weight is the probability weight of each type of move.
Point is the point about which to translate and rotate chains; default is
  the centroid of CenterChain, or overall centroid if CenterChain is None.
"""
    #check for indices and probabilities
    if ResProbPhiPsi is None:
      ResProbPhiPsi = zeros(len(self.Res), float)
      if ResIndPhiPsi is None:
        ResProbPhiPsi += 1.
      else:
        ResProbPhiPsi[ResIndPhiPsi] = 1.
    if ResProbChi is None:
      ResProbChi = zeros(len(self.Res), float)
      if ResIndChi is None:
        ResProbChi += 1.
      else:
        ResProbChi[ResIndChi] = 1.
    if ChainProb is None:
      ChainProb = zeros(len(self.Chains), float)
      if ChainInd is None:
        ChainProb += 1.
      else:
        ChainProb[ChainInd] = 1.
      if not CenterChain is None: ChainProb[CenterChain] = 0.
    if Point is None:
      self.Center(ChainNum = CenterChain)
      Point = self.Centroid(ChainNum = CenterChain)
    #make the move classes
    MCPhiPsi = pfunc.MCPhiPsiMoveClass(self, ResProb = ResProbPhiPsi, CenterChain = CenterChain)
    MCChi = pfunc.MCChiMoveClass(self, ResProb = ResProbChi)
    MCChainTrans = pfunc.MCChainTranslateClass(self, ChainProb = ChainProb, Point = Point)
    MCChainRot1 = pfunc.MCChainRotateClass(self, ChainProb = ChainProb, Point = Point)
    MCChainRot2 = pfunc.MCChainRotateClass(self, ChainProb = ChainProb, Point = None)
    MCMoves = [MCPhiPsi, MCChi, MCChainTrans, MCChainRot1, MCChainRot2]
    #adjust move weights as specified
    MCPhiPsi.W = WeightPhiPsi
    MCChi.W = WeightChi
    MCChainTrans.W = 0.33333 * WeightChain
    MCChainRot1.W = 0.33333 * WeightChain
    MCChainRot2.W = 0.33333 * WeightChain
    #prune for dead moves
    MCMoves = [m for m in MCMoves if not m.W * m.Active() == 0]
    return MCMoves

  def RunMC(self, ScoreFn, NSteps, MCMoves,
            Ti = 0.01, Tf = None, TargetAcc = 0.5, 
            StepFn = None, StepsUpdateMax = None, StepsUpdateTemp = None,
            KeepMin = False, Verbose = False):
    """Runs a MC algorithm using the current protein.
ScoreFn must be a function of two variables, a ProteinClass object and
  a residue number; residue number will not be None only if the side chain
  of the given residue is changed and nothing else.
NSteps is the total number of move attempts.
MCMoves is a list of MCMoveClass objects, initialized to the current protein.
Ti,Tf are the initial and final temperatures.
TargetAcc is the overall target acceptance ratio, for temperature adjustment.
StepFn will be run after each step and will be sent the following arguments:
  the ProteinClass, the step number, the current T, and a list of the current
  move classes.
StepsUpdateMax gives the interval for updating all of the maximum deflections.
StepsUpdateTemp gives the interval for updating the temperature.
KeepMin will remember the minimum scoring structure, substituting it at the end.
Verbose will provide detailed stats about the simulation.
"""
    #setup temperatures
    if Tf is None: Tf = Ti
    lnTi = log(Ti)
    dlnT = (log(Tf) - log(Ti)) / float(NSteps - 1)
    T = Ti
    #make the move probabilities
    x = sum([m.Active() * m.W for m in MCMoves])
    if x == 0:
      raise IndexError, "No available MC moves."
    MCMoveProbs = array([m.Active() * m.W / x for m in MCMoves], float)
    #prepare and perform the steps
    Ret = []
    OldE = ScoreFn(self, None)
    MinE = OldE
    MinEPos = self.Pos.copy()
    NAcc, NAtt = 0., 0.
    for m in MCMoves: m.Reset()
    for n in range(NSteps):
      #fraction done
      FracDone = float(n) / float(NSteps - 1)
      NAtt += 1.
      #get temp
      if not Tf == Ti: T = exp(lnTi + n*dlnT)
      #random move type
      m = MCMoves[pfunc.GetRandomInd(MCMoveProbs)]
      NewE, AccFact = m.Make(self, ScoreFn, OldE, FracDone = FracDone)
      Acc = exp(min(0., (OldE - NewE) / T + AccFact)) > random.random()
      m.Accept(self, Acc)
      if Acc:
        NAcc += 1.
        OldE = NewE
      #update parameters
      if not StepsUpdateMax is None and (n+1) % StepsUpdateMax == 0:
        for m in MCMoves:
          m.Update()
          m.Reset()
        if Verbose: print "Ds updated to [" + ",".join(["%.4f" % m.D for m in MCMoves]) + "]"
      if not StepsUpdateTemp is None and (n+1) % StepsUpdateTemp == 0:
        T = T * (TargetAcc + 0.1)/(NAcc/NAtt + 0.1)
        NAcc, NAtt = 0., 0.
        if Verbose: print "Temp updated to %f" % T
      #update minimum energy structure
      if KeepMin and OldE < MinE:
        MinE, MinEPos = OldE, self.Pos.copy()
      #run the step function
      if not StepFn is None:
        v = StepFn(self, n, T, MCMoves)
        if not v is None: Ret.append(v)
    #update structure
    if KeepMin: self.Pos = MinEPos
    #return the data
    if len(Ret) > 0: return Ret

    
#======== SPLICING, DICING, AND OPTIMIZING ========

  def Concat(self, Prot, NewChain = True):
    "Modifies self in place with a copy of Prot added without aligning bonds"
    #remember where new res start
    ResNum = len(self.Res)
    AtomNum = len(self.Pos)
    #first add residues
    self.Res.extend(copy.deepcopy(Prot.Res))
    #now add positions
    self.Pos = concatenate((self.Pos, Prot.Pos))
    #renumber
    self.Update()
    if NewChain: self.NewChain(ResInd = range(ResNum, len(self.Res)))
    return self

  def AlignPep(self, ResNum):
    """Modifies self in place to align the peptide bond to the ideal
    geometry between residues ResNum-1 and ResNum."""
    #rotate and align
    AtomNum = self.Res[ResNum].StartAtom
    PosC = self.Pos[self.AtomNum(ResNum-1, "C", AtomAliases), :]
    PosN = self.Pos[self.AtomNum(ResNum, "N", AtomAliases), :]
    PosCA0 = self.Pos[self.AtomNum(ResNum-1, "CA", AtomAliases), :]
    PosCA1 = self.Pos[self.AtomNum(ResNum, "CA", AtomAliases), :]
    PosNproj = self.ProjectN(ResNum-1)
    PosCproj = self.ProjectC(ResNum)
    #first align the N and C to make the peptide bond
    Vec, Ang = GetVecMapping(PosNproj - PosC, PosN - PosCproj)
    rm = RotMat(Vec, Ang)
    self.Pos[AtomNum:,:] = RotateAboutPoint(self.Pos[AtomNum:,:], rm, PosN) + PosNproj - PosN
    #next rotate the peptide bond
    Omega = Dihedral(PosCA0, PosC, PosN, PosCA1)
    rm = RotMat(PosC - PosN, 180. - Omega)
    self.Pos[AtomNum:,:] = RotateAboutPoint(self.Pos[AtomNum:,:], rm, PosN)

  def Join(self, Prot):
    "Modifies self in place with a copy of Prot spliced onto its end."
    #remember where new res start
    ResNum = len(self.Res)
    #update chains of prot
    if ResNum > 0:
      OldStartChainSelf = self.StartChain 
      OldStartChainProt = Prot.StartChain
      #have to change the start chain to A to get numbering right
      self.StartChain = "A"
      Prot.StartChain = self.Chains[-1]
    #first add residues
    self.Res.extend(copy.deepcopy(Prot.Res))
    #now add positions
    self.Pos = concatenate((self.Pos, Prot.Pos))
    #renumber
    self.Update()
    #update terminal atoms
    if ResNum > 0: self.RemoveTermini(ResInd = [ResNum - 1])
    if len(self.Res) > ResNum: self.RemoveTermini(ResInd = [ResNum])
    #put chains back
    if ResNum > 0:
      Prot.StartChain = OldStartChainProt
      self.StartChain = OldStartChainSelf
    #check if this was it
    if ResNum == 0 or len(Prot.Res) == 0:
      return self
    #now rotate and align
    self.AlignPep(ResNum)
    return self

  def Prejoin(self, Prot):
    "Modifies self in place with a copy of Prot spliced onto its beginning."
    #remember where old res start
    ResNum = len(Prot.Res)
    #renumber chains
    if ResNum > 0:
      OldStartChainSelf = self.StartChain 
      OldStartChainProt = Prot.StartChain
      #have to change the start chain to A to get numbering right
      Prot.StartChain = "A"
      self.StartChain = Prot.Chains[-1]
    #first add residues
    self.Res = copy.deepcopy(Prot.Res) + self.Res
    #now add positions
    self.Pos = concatenate((Prot.Pos, self.Pos))
    #renumber
    self.Update()
    #update terminal atoms
    if ResNum > 0: self.RemoveTermini(ResInd = [ResNum - 1])
    if len(self.Res) > ResNum: self.RemoveTermini(ResInd = [ResNum])
    #put chains back
    if ResNum > 0:
      Prot.StartChain = OldStartChainProt
      self.StartChain = OldStartChainSelf
    #check if this was it
    if ResNum == 0 or ResNum == len(self.Res):
      return self
    #now rotate and align
    self.AlignPep(ResNum)
    return self  

  def AddTermini(self, TermN = True, TermC = True):
    "Modifies self in place to add terminal atoms."
    self = self.Copy()
    r = self.Res[0]
    if TermN and not sequence.Cap(r.Name):
      #find the atoms
      i = r.AtomNum("H", NotFoundError = False, Absolute = True)
      j = r.AtomNum("N", NotFoundError = False, Absolute = True)
      k = r.AtomNum("CA", NotFoundError = False, Absolute = True)
      #modify nitrogen
      if i >= 0:
        self.Atoms[j].FormalCharge = 1.
      #modify hydrogens
      if i >= 0 and j >= 0 and k >= 0:
        #convert hydrogen
        self.Atoms[i].Name = "H1"
        self.Atoms[i].FullName = " H1 "
        #calculate vectors
        VecCAN = self.Pos[j] - self.Pos[k]
        VecNH = self.Pos[i] - self.Pos[j]
        PosN = self.Pos[j].copy()
        #make the other atoms
        Atom1 = AtomClass(Name = " H2 ", Element = "H", ResName = r.Name)
        rm = RotMat(VecCAN, -120.)
        Pos1 = RotateAboutPoint(self.Pos[i,:], rm, self.Pos[j,:])
        Atom2 = AtomClass(Name = " H3 ", Element = "H", ResName = r.Name)
        rm = RotMat(VecCAN, 120.)
        Pos2 = RotateAboutPoint(self.Pos[i,:], rm, self.Pos[j,:])
        self.AddAtom(0, Atom1, Pos = Pos1, Ind = i+1)
        self.AddAtom(0, Atom2, Pos = Pos2, Ind = i+2)
        #add displacements to make tetrahedral
        Vec = UnitVec(VecCAN)*Length(VecNH) - sum(self.Pos[i:i+3,:] - PosN, axis=0)
        self.Pos[i:i+3,:] = self.Pos[i:i+3,:] + Vec / 3.
    r = self.Res[-1]
    if TermC and not sequence.Cap(r.Name):
      #find the atoms
      i = r.AtomNum("O", NotFoundError = False, Absolute = True)
      j = r.AtomNum("C", NotFoundError = False, Absolute = True)
      k = r.AtomNum("CA", NotFoundError = False, Absolute = True)
      if i >= 0 and j >= 0 and k >= 0:
        #calculate vectors
        VecCAC = self.Pos[j] - self.Pos[k]
        #add an oxygen atom
        Atom1 = AtomClass(Name = " OXT", Element = "O", ResName = r.Name)
        Atom1.FormalCharge = -1
        rm = RotMat(VecCAC, -180.)
        Pos1 = RotateAboutPoint(self.Pos[i,:], rm, self.Pos[j,:])
        self.AddAtom(len(self.Res) - 1, Atom1, Pos = Pos1, Ind = i+1)
    return self

  def RemoveTermini(self, ResInd = None):
    "Modifies self in place to remove any terminal atoms."
    if ResInd is None: ResInd = [0, len(self) - 1]
    #remove termini atoms
    for r in ResInd:
      Res = self.Res[r]
      i = 0
      while i < len(Res.Atoms):
        AtomName = Res.Atoms[i].Name.strip()
        if AtomName in TerminiAtoms:
          self.DelResAtom(r, i)
        else:
          #convert one terminal H to the amide H
          if AtomName == "H1":
            Res.Atoms[i].Name = "H"
            Res.Atoms[i].FullName = " H  "
          i += 1
    return self

  def Decap(self):
    "Returns a copy of self with caps and termini atoms removed."
    p = self.Copy()
    #remove caps
    aa = sequence.AAInst(p.Res[0].Name)
    if not aa is None and len(p) > 0 and aa.Cap:
      p = p[1:]
    aa = sequence.AAInst(p.Res[-1].Name)
    if not aa is None and len(p) > 0 and aa.Cap:
      p = p[:-1]
    #remove termini atoms
    for r in [0, len(p) - 1]:
      Res = p.Res[r]
      i = 0
      while i < len(Res.Atoms):
        AtomName = Res.Atoms[i].Name.strip()
        if AtomName in TerminiAtoms:
          p.DelResAtom(r, i)
        else:
          #convert one terminal H to the amide H
          if AtomName == "H1":
            Res.Atoms[i].Name = "H"
            Res.Atoms[i].FullName = " H  "
          i += 1          
    return p

  def Cap(self, CapN = True, CapC = True,
          CapNRes = "ACE", CapCRes = "NME"):
    "Returns a copy of self with capping residues added."
    p = self.Copy()
    r = p.Res[0].Name
    if CapN and not (sequence.Cap(r) or sequence.Terminal(r)):
      p.Prejoin(ProteinClass(Seq = [CapNRes]))
    r = p.Res[-1].Name
    if CapC and not (sequence.Cap(r) or sequence.Terminal(r)):
      p.Join(ProteinClass(Seq = [CapCRes]))
    return p

  def MutateRes(self, ResNum, ResName, InPlace = False):
    "Returns a copy where residue number ResNum is mutated to ResName."
    #make new and old versions of the residue
    rn = Templates[ResName].Copy()
    #get the old and new bb atom positions
    BBAtoms = ["N", "CA", "C"]
    Poso = self.Pos.take([self.AtomNum(ResNum, a, AtomAliases)
                          for a in BBAtoms], axis=0)
    Posn = rn.Pos.take([rn.AtomNum(0, a, AtomAliases)
                        for a in BBAtoms], axis=0)
    #get CA position
    PosnCen = Posn.mean(axis=0)
    PosoCen = Poso.mean(axis=0) 
    #find the rotation matrix to take new to old
    rm = RotMatRMSD(Poso, Posn, True)
    #now rotate and translate the new residue, aligning at CA
    rn.Pos = dot(rn.Pos - PosnCen, rm) + PosoCen
    #make the phi, psi angles the same
    Phi, Psi = self.PhiPsi(ResNum)
    rn.RotateToPhiPsi(0, Phi, Psi)
    #replace the backbone atoms
    for a in BackboneAtoms:
      ano = self.AtomNum(ResNum, a, NotFoundError = False)
      ann = rn.AtomNum(0, a, NotFoundError = False)
      if ano >= 0 and ann >= 0:
        rn.Pos[ann] = self.Pos[ano]
    #now put the protein together
    if InPlace:
      p = self
    else:
      p = self.Copy()
    a1 = self.Res[ResNum].StartAtom
    a2 = self.Res[ResNum].StopAtom
    p.Pos = concatenate((self.Pos[:a1,:], rn.Pos, self.Pos[a2:,:]))
    rn.Res[0].Chain = p.Res[ResNum].Chain
    p.Res[ResNum] = rn.Res[0]
    p.Update()
    return p

  def MutateSeq(self, Seq):
    "Returns a copy where the sequence has been changed to Seq."
    Seq = sequence.SeqToList(Seq)
    n = min(len(Seq), len(self.Res))
    p = self.Copy()
    for i in range(0,n):
      p = p.MutateRes(i, Seq[i], InPlace = True)
    return p

  def Pep(self, ResNum, N = 5):
    "Rotates psi(ResNum-1) and phi(ResNum) to minimize overlaps."
    MinE = None
    OptPhi, OptPsi = 0., 0.
    for ((PhiStart,PsiStart),(PhiStop,PsiStop)) in PepAngleRanges:
      for i in range(0,N):
        Phi = PhiStart + (PhiStop-PhiStart)*(float(i)+0.5)/float(N)
        self.RotateToPhiPsi(ResNum, Phi = Phi)
        for j in range(0,N):
          Psi = PsiStart + (PsiStop-PsiStart)*(float(j)+0.5)/float(N)
          self.RotateToPhiPsi(ResNum-1, Psi = Psi)
          E = self.StericScore()
          if MinE is None or MinE > E:
            MinE = E
            OptPhi, OptPsi = Phi, Psi
    self.RotateToPhiPsi(ResNum-1, Psi = OptPsi)
    self.RotateToPhiPsi(ResNum, Phi = OptPhi)

  def OptimizeRes(self, ResNum, N = 5):
    "Rotates phi(ResNum) and psi(ResNum) to minimize overlaps."
    MinE = None
    OptPhi, OptPsi = 0., 0.
    for ((PhiStart,PsiStart),(PhiStop,PsiStop)) in ResAngleRanges:
      for i in range(0,N):
        Phi = PhiStart + (PhiStop-PhiStart)*(float(i)+0.5)/float(N)
        self.RotateToPhiPsi(ResNum, Phi = Phi)
        for j in range(0,N):
          Psi = PsiStart + (PsiStop-PsiStart)*(float(j)+0.5)/float(N)
          self.RotateToPhiPsi(ResNum, Psi = Psi)
          E = self.StericScore()
          if MinE is None or MinE > E:
            MinE = E
            OptPhi, OptPsi = Phi, Psi
    self.RotateToPhiPsi(ResNum, Phi = OptPhi, Psi = OptPsi)

  def OptimizeSC(self, ResInd = None, N = 100, Ti = 100.0, Tf = 0.01):
    """Optimizes chi angles to remove overlaps, using a MC algorithm."""
    if ResInd is None: ResInd = range(len(self))
    #make a scoring function
    def ScoreFn(q, rn): return q.LJScore(rn)
    #run the MC
    MCMoves = self.GetMCMoves(ResIndChi = ResInd, WeightPhiPsi = 0.,
                              WeightChi = 1., WeightChain = 0.)
    self.RunMC(ScoreFn, N*len(ResInd), MCMoves, Ti = Ti, Tf = Tf)

  def Optimize(self, ResInd = None, N = 100, Ti = 0.001, Tf = 0.0001,
               FracPhiPsiMove = 0.1, FracChainMove = 0., BBScale = 0.1):
    """Optimizes dihedrals and side chains to remove overlaps,
    using a MC algorithm.  Tries to maintain backbone configuration."""
    if ResInd is None: ResInd = range(len(self))
    #make a scoring function;
    #add harmonic restraint to current backbone
    ResPos = self.ResPos(ResAtom="CA")
    def ScoreFn(p, rn):
      if rn is None:
        BBScore = sum((ResPos - p.ResPos(ResAtom="CA"))**2, axis = None)
        return p.LJScore() + BBScale * BBScore
      else:
        return p.LJScore(rn)
    #run the MC
    MCMoves = self.GetMCMoves(ResIndPhiPsi = ResInd, ResIndChi = ResInd,
                              WeightPhiPsi = (1.-FracChainMove)*FracPhiPsiMove,
                              WeightChi =(1.-FracChainMove)*(1.-FracPhiPsiMove),
                              WeightChain = FracChainMove)
    self.RunMC(ScoreFn, N*len(ResInd), MCMoves, Ti = Ti, Tf = Tf)
    

#======== ALIGNING TO TEMPLATES ========

  def TemplateBonds(self, TemplateProtein = None):
    """Generates intraresidue bonds based on template residues."""
    for (rn, r) in enumerate(self.Res):
      if TemplateProtein is None:
        if not r.Name in Templates: continue
        tr = Templates[r.Name].Res[0]
      else:
        tr = TemplateProtein.Res[rn]
      r.ClearBonds()
      for (a,b) in tr.Bonds:
        Name1 = tr.Atoms[a].Name
        Name2 = tr.Atoms[b].Name
        i = r.AtomNum(Name1, NotFoundError = False)
        j = r.AtomNum(Name2, NotFoundError = False)
        if i >= 0 and j >= 0: r.AddBond((i,j))

  def TemplateAtoms(self, ResInd = None, AtomNames = [], Elements = [],
                    DelAtoms = True, AddAtoms = True, Verbose = False,
                    TemplateProtein = None):
    """Generates missing atoms and deletes extra atoms based on the
template library of residues.  Bond information in templates
is necessary to add atoms."""
    if ResInd is None: ResInd = range(len(self.Res))
    #delete any extra atoms
    if DelAtoms:
      for ResNum in ResInd:
        r = self.Res[ResNum]
        if TemplateProtein is None:
          if not r.Name in Templates: continue
          tr = Templates[r.Name].Res[0]
        else:
          tr = TemplateProtein.Res[ResNum]
        DelList = []
        for i in range(len(r.Atoms)):
          Name, Element = r.Atoms[i].Name.strip(), r.Atoms[i].Element
          if tr.HasAtom(Name): continue
          #check if specific names or elements were specified
          if not Name in AtomNames and len(AtomNames) > 0: continue
          if not Element in Elements and len(Elements) > 0: continue
          #add the extra atom to the delete list
          DelList.append(i - len(DelList))
          if Verbose:
            print "Deleting atom %s in residue %s%d." % (Name, r.Name, ResNum)
        #delete the extra atoms
        for i in DelList:
          self.DelResAtom(ResNum, i)
    #now add any missing atoms
    if AddAtoms:
      for ResNum in ResInd:
        r = self.Res[ResNum]
        if TemplateProtein is None:
          if not r.Name in Templates: continue
          tp = Templates[r.Name]
          tr = tp.Res[0]
        else:
          tp = TemplateProtein
          tr = tp.Res[ResNum]
        Name = ""
        for i in range(len(tr.Atoms)):
          LastName = Name
          Name, Element = tr.Atoms[i].Name.strip(), tr.Atoms[i].Element
          if r.HasAtom(Name): continue
          #check if specific names or elements were specified
          if not Name in AtomNames and len(AtomNames) > 0: continue
          if not Element in Elements and len(Elements) > 0: continue
          #check for special case backbone atoms that need a peptide bond to align
          if Name == "H" and not r.Name in NoDihedrals:
            Pos = self.ProjectH(ResNum)
          elif Name == "O" and not r.Name in NoDihedrals:
            Pos = self.ProjectO(ResNum) 
          else:
            #get at least three bonded atoms
            Depth = 1
            Bonded = []
            while True:
              #get bonds of this depth and check if we ran out of atoms
              OldBonded = Bonded
              Bonded = tr.GetBonded(i, Depth = Depth)
              if OldBonded == Bonded: break
              an, tan = [], []
              for j in Bonded:
                k = r.AtomNum(tr.Atoms[j].Name, NotFoundError = False)
                if k < 0: continue
                an.append(k + r.StartAtom)
                tan.append(j + tr.StartAtom)
              #try to get at least three bonded atoms;
              #increase the depth if necessary
              if len(an) < 3:
                Depth += 1
              else:
                break
            #cant continue if less than one atom to map
            if len(an) < 1: continue
            #now make the position arrays for the atoms
            Pos = self.Pos.take(an, axis=0)
            tPos = tp.Pos.take(tan, axis=0)
            tCen = tPos.mean(axis = 0)
            #get the rotation matrix which minimizes rmsd for bonded atoms
            RM = RotMatRMSD(Pos, tPos, Center = True)
            #rotate the missing atom position and the nearest bonded
            tPos = tp.Pos.take([i + tr.StartAtom, tan[0]], axis=0)
            tPos = RotateAboutPoint(tPos, RM, tCen)
            #make the correct dist from the first bonded atom
            Pos = tPos[0,:] + (Pos[0,:] - tPos[1,:])
          #find where to put the new atom
          if LastName == "":
            Ind = 0
          else:
            Ind = r.AtomNum(LastName) + 1
          #add the new atom
          Atom = copy.deepcopy(tr.Atoms[i])
          self.AddAtom(ResNum, Atom, Pos = Pos, Ind = Ind)
          if Verbose:
            print "Added atom %s to residue %s%d." % (Atom.Name, r.Name, ResNum)
        
  def Template(self, TemplateProtein = None):
    """Matches any residues to templates, fills in missing bond
information, and adds missing atoms."""
    self.TemplateAtoms(TemplateProtein = TemplateProtein)
    self.TemplateBonds(TemplateProtein = TemplateProtein)

  def Hydrogenate(self, BackboneOnly = False):
    """Adds hydrogens if they are missing."""
    if BackboneOnly:
      self.TemplateAtoms(AtomNames = ["H"])
    else:
      self.TemplateAtoms(Elements = ["H"])
    self.TemplateBonds()

  def Dehydrogen(self):
    """Removes all hydrogens."""
    Ind = []
    for i in range(len(self.Atoms)):
      if self.Atoms[i].Element == "H":
        Ind.append(i - len(Ind))
    for i in Ind:
      self.DelAtom(i) 
    

#======== COARSE GRAINING ========

  def CoarseGrained(self, AlignToCB = False, CentroidDist = None,
                    RgScale = None):
    """Returns a copy of self in which side chains are replaced
with a single large atom."""
    if RgScale is None: RgScale = CoarseGrainRgScaleDflt
    p = self.Copy()
    ResList = [i for i in range(len(p)) if p.HasAtom(i, "CB")]
    for i in ResList:
      r = p.Res[i]
      j = 0
      Pos, PosSq = zeros(3, float), zeros(3, float)
      m = 0
      while j < len(r.Atoms):
        if r.Atoms[j].Name == "CB":
          CBind = j
          j += 1
        elif r.Atoms[j].Name in BackboneAtoms:
          j += 1
        else:
          Pos = Pos + p.Pos[r.StartAtom + j]
          PosSq = PosSq + p.Pos[r.StartAtom + j]**2
          m += 1
          p.DelAtom(r.StartAtom + j)
      if m > 0:
        Pos, PosSq = Pos / float(m), PosSq / float(m)
        Rg = sqrt(sum(PosSq - Pos*Pos))
        r.Atoms[CBind].Name = "cen"
        r.Atoms[CBind].FullName = "cen"
        r.Atoms[CBind].Element = "cen"
        r.Atoms[CBind].Radius = Rg * RgScale
        PosCA = p.Pos[r.AtomNum("CA", Absolute = True)]
        l = Length(Pos - PosCA)
        v = (Pos - PosCA) / l
        if AlignToCB:
          v = UnitVec(p.Pos[r.StartAtom + CBind] - PosCA)
        if not CentroidDist is None:
          l = CentroidDist
        p.Pos[r.StartAtom + CBind] = PosCA + l*v
    p.Update()
    return p

  def Reduced(self, SideChainAtoms = ["CB"]):
    """Returns a copy of self in which all side chain atoms
are removed except those specified (default is CB)."""
    p = self.Copy()
    ResList = [i for i in range(len(p)) if p.HasAtom(i, "CB")]
    for i in ResList:
      r = p.Res[i]
      j = 0
      while j < len(r.Atoms):
        if r.Atoms[j].Name in BackboneAtoms + SideChainAtoms:
          j += 1
        else:
          p.DelAtom(r.StartAtom + j)
    p.Update()
    return p



#======== FILE READING/WRITING ========

  def ReadPdb(self, Pdb, AssignBonds = None):
    "Reads informaton from a Pdb file or data string."
    if AssignBonds is None: AssignBonds = AssignBondsDflt
    self.Res = []
    self.Pos = []
    #get the data
    if not "\n" in Pdb:
      if os.path.isfile(Pdb):
        #check for gzipped
        if Pdb.endswith(".gz"):
          Pdb = gzip.GzipFile(Pdb,"r").read()
        else:
          Pdb = file(Pdb,"r").read()
      else:
        raise IOError, "Pdb file %s not found." % Pdb
        return
    #parse the data into residues first
    ResID, Chain = None, None
    cn = copy.deepcopy(ChainNames)
    #the following two dictionaries give the atom number
    #in the residue and the residue class:
    PdbAtomNums = {}
    PdbAtomRes = {}
    for l in Pdb.split("\n"):
      if l.startswith("ATOM"):
        #check if atom name is too long;
        #just a fix for baker's decoy set
        if l[16] in "0123456789" and not " " in l[13:16] and l[12] == " ":
          l = l[:12] + l[16] + l[13:16] + l[17:]
        #first get the atom data
        Name = l[12:16]
        Element = Name[1]
        try:
          BFactor = float(l[60:66])
        except StandardError:
          BFactor = 0.
        try:
          Occupancy = float(l[54:60])
        except StandardError:
          Occupancy = 1.
        try:
          Charge = l[78:80]
          if Charge[1] in "+-": Charge = Charge[1] + Charge[0]
          Charge = int(Charge)
        except:
          Charge = None
        #get the residue information
        ResID, LastResID = l[22:27], ResID
        Chain, LastChain = l[21], Chain
        ResName = l[17:20]
        if not Chain == LastChain:
          CurChain = cn.pop(0)
        if not ResID == LastResID:
          Res = ResClass(Name = l[17:20], Chain = CurChain)
          self.Res.append(Res)
        #append the atom
        Atom = AtomClass(Name = Name, Element = Element, BFactor = BFactor,
                         FormalCharge = Charge, Occupancy = Occupancy,
                         ResName = ResName)
        self.Res[-1].Atoms.append(Atom)
        #now get the coordinates
        try:
          Pos = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
        except StandardError:
          Pos = [0., 0., 0.]
        self.Pos.append(Pos)
        #add the atom number and residue object to the dictionary
        PdbAtomNums[l[6:11]] = len(self.Res[-1].Atoms) - 1
        PdbAtomRes[l[6:11]] = self.Res[-1]
      elif l.startswith("TER"):
        #chain termination
        ResID, Chain = None, None
      elif l.startswith("CONECT"):
        #get a bond
        Atom1 = PdbAtomNums.get(l[6:11], -1)
        Res1 = PdbAtomRes.get(l[6:11], None)
        if Atom1 >= 0:
          Atom2List = [PdbAtomNums.get(l[i:i+5], -1) for i in [11,16,21,26]]
          Res2List = [PdbAtomRes.get(l[i:i+5], None) for i in [11,16,21,26]]
          for i in range(4):
            Atom2, Res2 = Atom2List[i], Res2List[i]
            if Atom2 < 0 or not Res1 is Res2: continue
            Res1.AddBond((Atom1,Atom2))
      elif l.startswith("ENDMDL"):
        break
    #now renumber and convert Pos to an array
    self.Pos = array(self.Pos, float)
    self.Update()
    #check if we need to find bonds
    if AssignBonds: self.GenerateBonds()
    

  def GetPdb(self, UseCharge = False, UseBonds = False):
    "Returns a string of pdb-formatted data."
    Pdb = "MODEL\n"
    an = -1
    for (rn, r) in enumerate(self.Res):
      if rn in self.ChainResNums[1:]:
        Pdb += "TER\n"
      for a in r.Atoms:
        an += 1
        (x, y, z) = self.Pos[an,:]
        s = "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f %5.2f%6.2f" \
          % (an + 1, a.FullName, r.FullName, r.Chain, rn + 1, x, y, z,
             a.Occupancy, a.BFactor)
        if UseCharge:
          if not a.FormalCharge is None:
            s += "            %1d" % abs(a.FormalCharge)
            if a.FormalCharge > 0:
              s += "+"
            elif a.FormalCharge < 0:
              s += "-"
            else:
              s += " "
        Pdb += s + "\n"
    Pdb += "TER\n"
    if UseBonds:
      for r in self.Res:
        for (a,b) in r.Bonds:
          Pdb += "%6s%5d%5d\n" % ('CONECT', a + r.StartAtom + 1, b + r.StartAtom + 1)
    Pdb += "ENDMDL\n"
    return Pdb
  
  def WritePdb(self, PdbFile, UseCharge = False, UseBonds = False):
    "Writes data to a pdb file."
    file(PdbFile, "w").write(self.GetPdb(UseCharge, UseBonds))

  def ReadPrmtop(self, PrmtopFile):
    """Reads from a prmtop file. Coordinates must be supplied 
    by a subsequent statement self.Pos = Pos."""
    #get names and convert names to pdb style
    AtomNames = coords.GetPrmtopAtomNames(PrmtopFile)
    AtomNames = coords.AmbToPdbAtomNames(AtomNames)
    AtomRes = coords.GetPrmtopAtomRes(PrmtopFile)
    Seq = coords.GetPrmtopSeq(PrmtopFile)
    NAtom = len(AtomNames)
    #make the residues
    self.Res = []
    for r in Seq:
      self.Res.append(ResClass(Name = r))
    #parse atoms into residues
    for i in range(0, NAtom):
      Name = AtomNames[i]
      Element = Name[1]
      ResName = Seq[AtomRes[i]]
      Atom = AtomClass(Name = Name, Element = Element, ResName = ResName)
      self.Res[AtomRes[i]].Atoms.append(Atom)
    #initialize the positions to zero
    self.Pos = zeros((NAtom,3), float)
    self.Update()
    #check for chains
    self.GenerateChains()       

  def LinkTrj(self, Trj):
    """Links to an amber trajectory object which will
    automatically update self.Pos.  Note that assignments of
    the form self.Pos = x will break the link."""
    #forge the links
    self.Trj = Trj
    #read in the atom names, etc
    self.ReadPrmtop(Trj.PrmtopFile)
    #link the positions
    Trj.LinkPos = self.Pos

  def LinkCoordsObj(self, CoordsObj):
    """Links to a generic coords object which will
    automatically update self.Pos.  Note that assignments of
    the form self.Pos = x will break the link.  The
    coords object must supply member arrays Seq, AtomNames,
    and AtomRes."""
    NAtom = len(CoordsObj.AtomNames)
    #make the residues
    self.Res = []
    for r in CoordsObj.Seq:
      self.Res.append(ResClass(Name = r))
    #parse atoms into residues
    for i in range(0, NAtom):
      Name = CoordsObj.AtomNames[i]
      ResName = CoordsObj.Seq[CoordsObj.AtomRes[i]]
      Atom = AtomClass(Name = Name, ResName = ResName)
      self.Res[CoordsObj.AtomRes[i]].Atoms.append(Atom)
    #initialize the positions to zero
    self.Pos = zeros((NAtom,3), float)
    self.Update()
    #link the positions and trajectory
    CoordsObj.LinkPos = self.Pos
    self.CoordsObj = CoordsObj

  def UnlinkTrj(self):
    "Unlinks from an amber trajectory class."
    self.Trj.LinkPos = None
    del self.Trj

  def UnlinkCoordsObj(self):
    "Unlinks from an generic coords obj class."
    self.CoordsObj.LinkPos = None
    del self.CoordsObj
    

#======== GEOMETRY AND ANGLES ========

  def Centroid(self, ResNum = None, ChainNum = None):
    "Returns the centroid of the whole structure, a residue, or a chain."
    if not ResNum is None:
      return proteinlib.centroidrange(self.Pos, self.Res[ResNum].StartAtom,
                                      self.Res[ResNum].StopAtom)
    elif not ChainNum is None:
      return proteinlib.centroidrange(self.Pos, self.ChainAtomNums[ChainNum],
                                      self.ChainAtomNums[ChainNum+1])
    else:
      return proteinlib.centroidall(self.Pos)    

  def ResPos(self, ResInd = None, ResAtom = None):
    """Returns a compressed position array with one entry for each residue.
    Can use ResAtom = '*' to do centroid.
    ResInd can be a list of indices for residues."""
    #check for empty specification
    if len(self.Res) == 0: return zeros((0,3), float)
    if not ResInd is None and len(ResInd) == 0: return zeros((0,3), float)
    #check for library versus normal routines
    if USELIB:
      if ResAtom is None: ResAtom = self.ResAtom
      if ResAtom == "*":
        if ResInd is None:
          N = len(self.Res)
          return proteinlib.resposcentroid(self.Pos, self.AtomResNum, N)
        else:
          ResInd = array(ResInd, int)
          N = len(self.Res)
          return proteinlib.resposcentroidind(self.Pos, self.AtomResNum, N, ResInd)
      elif ResAtom == self.ResAtom:
        Pos = self.Pos.take(self.__ResAtomNum, axis=0)
        if not ResInd is None: Pos = Pos.take(ResInd, axis=0)
        return Pos
      else:
        if ResInd is None: ResInd = range(0, len(self.Res))
        Pos = zeros((len(ResInd), 3), float)
        for i in range(0, len(ResInd)):
          an = self.AtomNum(ResInd[i], ResAtom, AtomAliases)
          Pos[i,:] = self.Pos[an,:]
        return Pos
    else:
      if ResAtom is None: ResAtom = self.ResAtom
      if ResAtom == "*":
        if ResInd is None: ResInd = range(0, len(self.Res))
        Pos = zeros((len(ResInd), 3), float)
        for i in range(0, len(ResInd)):
          Pos[i,:] = self.Centroid(ResInd[i])
      elif ResAtom == self.ResAtom:
        Pos = self.Pos.take(self.__ResAtomNum, axis=0)
        if not ResInd is None: Pos = Pos.take(ResInd, axis=0)
      else:
        if ResInd is None: ResInd = range(0, len(self.Res))
        Pos = zeros((len(ResInd), 3), float)
        for i in range(0, len(ResInd)):
          an = self.AtomNum(ResInd[i], ResAtom, AtomAliases)
          Pos[i,:] = self.Pos[an,:]
      return Pos

  def Center(self, ResNum = None, ChainNum = None):
    "Places the centroid (optionally of a residue or chain) at the origin."
    self.Pos = proteinlib.center(self.Pos, self.Centroid(ResNum = ResNum, ChainNum = ChainNum))
      
  def Rotate(self, Vec, Ang, Point = None, ChainNum = None):
    "Rotates the whole structure or a chain."
    if Point is None: Point = self.Centroid()
    rm = RotMat(Vec, Ang)
    if not ChainNum is None:
      a1, a2 = self.ChainAtomNums[ChainNum:ChainNum+2]
      self.Pos[a1:a2,:] = RotateAboutPoint(self.Pos[a1:a2,:], rm, Point)
    else:
      self.Pos = RotateAboutPoint(self.Pos, rm, Point)

  def Translate(self, Vec, ChainNum = None):
    "Tranlates the whole structure or a chain."
    if not ChainNum is None:
      a1, a2 = self.ChainAtomNums[ChainNum:ChainNum+2]
      self.Pos[a1:a2,:] = self.Pos[a1:a2,:] + Vec
    else:
      self.Pos = self.Pos + Vec

  def ProjectC(self, ResNum):
    "Projects the position of the preceding C."
    PosCA = self.Pos[self.AtomNum(ResNum, "CA", AtomAliases), :]
    PosN = self.Pos[self.AtomNum(ResNum, "N", AtomAliases), :]
    #check if there are hydrogens
    if self.HasAtom(ResNum, "H", AtomAliases):
      PosH = self.Pos[self.AtomNum(ResNum, "H", AtomAliases), :]
      return PosN + UnitVec((1-NCFracHN)*UnitVec(PosN-PosCA) +
        NCFracHN*UnitVec(PosN-PosH)) * CNBondLen
    else:
      #fudge it based on C
      PosC = self.Pos[self.AtomNum(ResNum, "C", AtomAliases), :]
      return PosN + UnitVec(PosCA - PosC) * CNBondLen

  def ProjectN(self, ResNum):
    "Projects the position of the following N."
    PosCA = self.Pos[self.AtomNum(ResNum, "CA", AtomAliases), :]
    PosC = self.Pos[self.AtomNum(ResNum, "C", AtomAliases), :]
    if self.HasAtom(ResNum, "O"):
      PosO = self.Pos[self.AtomNum(ResNum, "O", AtomAliases), :]
      return PosC + UnitVec((1-CNFracOC)*UnitVec(PosC-PosCA) +
        CNFracOC*UnitVec(PosC-PosO)) * CNBondLen
    else:
      #fudge it based on N
      PosN = self.Pos[self.AtomNum(ResNum, "N", AtomAliases), :]
      return PosC + UnitVec(PosCA - PosN) * CNBondLen 

  def ProjectH(self, ResNum):
    "Projects the position of H based on the peptide bond."
    PosCA = self.Pos[self.AtomNum(ResNum, "CA", AtomAliases), :]
    PosN = self.Pos[self.AtomNum(ResNum, "N", AtomAliases), :]
    n = self.AtomNum(ResNum - 1, "C", AtomAliases, NotFoundError = False)
    if n < 0:
      #fudge it based on C
      PosC = self.Pos[self.AtomNum(ResNum, "C", AtomAliases), :]
      return PosN + NHBondLen * UnitVec(NHFracCAC * UnitVec(PosC - PosCA)
             + (1 - NHFracCAC) * UnitVec(PosN - PosCA))
    else:
      PosC = self.Pos[n, :]
      return PosN + NHBondLen * UnitVec(NHFracCN * UnitVec(PosN - PosC)
             + (1 - NHFracCN) * UnitVec(PosN - PosCA))

  def ProjectO(self, ResNum):
    "Projects the position of O based on the peptide bond."
    PosCA = self.Pos[self.AtomNum(ResNum, "CA", AtomAliases), :]
    PosC = self.Pos[self.AtomNum(ResNum, "C", AtomAliases), :]
    n = self.AtomNum(ResNum + 1, "N", AtomAliases, NotFoundError = False)
    if n < 0:
      #fudge it based on N
      PosN = self.Pos[self.AtomNum(ResNum, "N", AtomAliases), :]
      return PosC + COBondLen * UnitVec(COFracCAC * UnitVec(PosC - PosCA)
             + (1 - COFracCAC) * UnitVec(PosN - PosCA))
    else:
      PosN = self.Pos[n, :]
      return PosC + COBondLen * UnitVec(COFracNC * UnitVec(PosC - PosN)
             + (1 - COFracNC) * UnitVec(PosC - PosCA))

  def PhiPsi(self, ResNum):
    "Returns the phi and psi angles."
    if self.Res[ResNum].Name in NoDihedrals: return None, None
    PosN = self.Pos[self.AtomNum(ResNum, "N", AtomAliases), :]
    PosCA = self.Pos[self.AtomNum(ResNum, "CA", AtomAliases), :]
    PosC = self.Pos[self.AtomNum(ResNum, "C", AtomAliases), :]
    #get starting and stopping atom numbers for this chain
    (r1,r2), (a1,a2) = self.ChainRange(ResNum)
    #get the preceding C
    if ResNum > r1:
      i = self.AtomNum(ResNum-1, "C", AtomAliases, NotFoundError = False)
      if i < 0:
        PosC0 = self.ProjectC(ResNum)
      else:
        PosC0 = self.Pos[i, :]
    else:
      #calculate where the preceeding C would be
      PosC0 = self.ProjectC(ResNum)
    #get the following N
    if ResNum < r2 - 1:
      i = self.AtomNum(ResNum+1, "N", AtomAliases, NotFoundError = False)
      if i < 0:
        PosN2 = self.ProjectN(ResNum)
      else:
        PosN2 = self.Pos[i, :]
    else:
      #calculate where the following N would be
      PosN2 = self.ProjectN(ResNum)
    #calculate the angles
    Phi = Dihedral(PosC0, PosN, PosCA, PosC)
    Psi = Dihedral(PosN, PosCA, PosC, PosN2)
    return Phi, Psi

  def RotatePhi(self, ResNum, Ang):
    "Rotates the structure Ang degrees along the phi angle."
    Res = self.Res[ResNum]
    if Res.Name in NoDihedrals or Res.Name in FixedPhi: return False
    CA = self.AtomNum(ResNum, "CA", AtomAliases)
    N = self.AtomNum(ResNum, "N", AtomAliases)
    Vec = self.Pos[CA,:] - self.Pos[N,:]
    Point = self.Pos[CA,:]
    #get starting and stopping atom numbers for this chain
    (r1,r2), (a1,a2) = self.ChainRange(ResNum)
    #rotate preceding residues and then current
    rm = RotMat(Vec, Ang)
    if Res.StartAtom > a1:
      self.Pos[a1:Res.StartAtom,:] = RotateAboutPoint(self.Pos[a1:Res.StartAtom,:], rm, Point)
    for Atom in ["N", "H"]:
      if self.HasAtom(ResNum, Atom, AtomAliases):
        AtomNum = self.AtomNum(ResNum, Atom, AtomAliases)
        self.Pos[AtomNum,:] = RotateAboutPoint(self.Pos[AtomNum,:], rm, Point)
    return True

  def RotatePsi(self, ResNum, Ang):
    "Rotates the structure Ang degrees along the psi angle."
    Res = self.Res[ResNum]
    if Res.Name in NoDihedrals or Res.Name in FixedPsi: return False
    CA = self.AtomNum(ResNum, "CA", AtomAliases)
    C = self.AtomNum(ResNum, "C", AtomAliases)
    Vec = self.Pos[CA,:] - self.Pos[C,:]
    Point = self.Pos[CA,:]
    #get starting and stopping atom numbers for this chain
    (r1,r2), (a1,a2) = self.ChainRange(ResNum)
    #rotate following residues and then current
    n = Res.StopAtom
    rm = RotMat(Vec, Ang)
    if n < a2 - 1:
      self.Pos[n:a2,:] = RotateAboutPoint(self.Pos[n:a2,:], rm, Point)
    for Atom in ["C", "O"]:
      if self.HasAtom(ResNum, Atom, AtomAliases):
        AtomNum = self.AtomNum(ResNum, Atom, AtomAliases)
        self.Pos[AtomNum,:] = RotateAboutPoint(self.Pos[AtomNum,:], rm, Point)
    return True

  def RotateToPhiPsi(self, ResNum, Phi = None, Psi = None):
    """Rotates the Phi and Psi angles to a specified value.
       Returns the Phi, Psi pair after the rotation."""
    CurPhi, CurPsi = self.PhiPsi(ResNum)
    if not Phi is None and not CurPhi is None:
      if self.RotatePhi(ResNum, Phi - CurPhi): CurPhi = Phi
    if not Psi is None and not CurPsi is None:
      if self.RotatePsi(ResNum, Psi - CurPsi): CurPsi = Psi
    return Phi, Psi

  def Chi(self, ResNum):
    "Returns the chi angles."
    Angs = []
    r = self.Res[ResNum]
    for ChiAtoms in r.ChiAtoms:
      Pos = self.Pos.take(ChiAtoms[:4] + r.StartAtom, axis=0)
      Angs.append(Dihedral(Pos[0], Pos[1], Pos[2], Pos[3]))
    return Angs

  def RotateChi(self, ResNum, ChiNum, Ang):
    "Rotates the structure Ang degrees along the specified chi angle."
    r = self.Res[ResNum]
    if ChiNum in FixedChi.get(r.Name, []): return False
    ChiAtoms = r.ChiAtoms[ChiNum] + r.StartAtom
    Vec = self.Pos[ChiAtoms[1],:] - self.Pos[ChiAtoms[2],:]
    Point = self.Pos[ChiAtoms[1],:]
    #rotate associated atoms
    Pos = self.Pos[ChiAtoms[4:]]
    rm = RotMat(Vec, Ang)
    Pos = RotateAboutPoint(Pos, rm, Point)
    self.Pos[ChiAtoms[4:]] = Pos
    return True

  def RotateToChi(self, ResNum, ChiNum, Ang):
    """Rotates the chi angle to a specified value.
       Returns the angle after the rotation."""
    CurChi = self.Chi(ResNum)[ChiNum]
    if self.RotateChi(ResNum, ChiNum, Ang - CurChi):
      return Ang
    else:
      return CurChi

  def HasOverlap(self, OverlapDist = 0.65):
    """Returns True if there is an atomic overlap."""
    for i in range(len(self.Pos)):
      for j in range(i+1, len(self.Pos)):
        if Length(self.Pos[i] - self.Pos[j]) <= OverlapDist:
          return True
    return False
        

#======== SCORES AND ENERGIES ========

  def __PrepBonds(self):
    """This prepares the proteinclass for energies by building bond information."""
    #updates bond info... must be done AFTER chain update
    if self.UseBonds and self.__Bonds23 is None:
      N = len(self.Atoms)
      self.__Bonds23, self.__Bonds4 = pfunc.GetBonds(self.Bonds())
      self.__Bonds23 = pfunc.ConvertBonds(self.__Bonds23, N)
      self.__Bonds4 = pfunc.ConvertBonds(self.__Bonds4, N)

  def StericScore(self, ResNum = None):
    """Returns an energy based on sum_ij (sig_ij/r_ij)^12 where
    sig_ij is based on the radii of the atoms involved.
    If supplied, ResNum will return the contribution from atoms in
    residue number ResNum with themselves and other residues."""
    if self.UseBonds: self.__PrepBonds()
    if USELIB:
      if self.UseBonds:
        if ResNum is None:
          return proteinlib.stericscorebond(self.Pos, self.Radius, self.StericDistCut,
                                            self.__Bonds23, self.__Bonds4, FFSteric14Scale)
        else:
          return proteinlib.stericscoreresbond(self.Pos, self.Radius, self.StericDistCut,
                                               ResNum, self.AtomResNum,
                                               self.__Bonds23, self.__Bonds4, FFSteric14Scale)  
      else:
        if ResNum is None:
          return proteinlib.stericscore(self.Pos, self.Radius, self.StericDistCut)
        else:
          return proteinlib.stericscoreres(self.Pos, self.Radius, self.StericDistCut,
                                           ResNum, self.AtomResNum)    
    else:
      N = len(self.Pos)
      E = 0.
      DistCutSq = self.StericDistCut**2
      if ResNum is None:
        AtomRange = range(0, N-1)
      else:
        r = self.Res[ResNum]
        a1, a2 = r.StartAtom, r.StopAtom
        AtomRange = range(a1, a2)
      MaskA = ones(N,bool)
      for i in AtomRange:
        MaskA[i] = False
        Vecs = self.Pos[MaskA] - self.Pos[i,:]
        Vecs = (Vecs*Vecs).sum(axis=1)
        if self.UseBonds:
          MaskNot23, Mask4 = pfunc.GetBondMasks(i, N, self.__Bonds23, self.__Bonds4)
          MaskB = logical_and(Vecs < DistCutSq, MaskNot23[MaskA])
        else:
          MaskB = Vecs < DistCutSq
        Rad = self.Radius[MaskA][MaskB] + self.Radius[i]
        Rad = Rad*Rad
        Vecs = Rad / Vecs.compress(MaskB)
        Vecs = Vecs*Vecs
        Vecs = Vecs*Vecs*Vecs
        if self.UseBonds: Vecs[Mask4[MaskA][MaskB]] *= FFSteric14Scale
        E += Vecs.sum()          
      return E

  def LJScore(self, ResNum = None):
    """Returns an energy based on LJ repulsion/dispersion interactions.
    If supplied, ResNum will return the contribution from atoms in
    residue number ResNum with themselves and other residues."""
    if self.UseBonds: self.__PrepBonds()
    if USELIB:
      if self.UseBonds:
        if ResNum is None:
          return proteinlib.ljscorebond(self.Pos, self.Radius, self.SqrtEps, self.LJDistCut,
                                        self.__Bonds23, self.__Bonds4, FFLJ14Scale)
        else:
          return proteinlib.ljscoreresbond(self.Pos, self.Radius, self.SqrtEps,
                                           self.LJDistCut, ResNum, self.AtomResNum,
                                           self.__Bonds23, self.__Bonds4, FFLJ14Scale)
      else:
        if ResNum is None:
          return proteinlib.ljscore(self.Pos, self.Radius, self.SqrtEps, self.LJDistCut)
        else:
          return proteinlib.ljscoreres(self.Pos, self.Radius, self.SqrtEps,
                                       self.LJDistCut, ResNum, self.AtomResNum)        
    else:
      N = len(self.Pos)
      E = 0.
      DistCutSq = self.LJDistCut**2
      if ResNum is None:
        AtomRange = range(0, N-1)
      else:
        r = self.Res[ResNum]
        a1, a2 = r.StartAtom, r.StopAtom
        AtomRange = range(a1, a2)
      MaskA = ones(N,bool)
      for i in AtomRange:
        MaskA[i] = False
        Vecs = self.Pos[MaskA] - self.Pos[i,:]
        Vecs = (Vecs*Vecs).sum(axis=1)
        if self.UseBonds:
          MaskNot23, Mask4 = pfunc.GetBondMasks(i, N, self.__Bonds23, self.__Bonds4)
          MaskB = logical_and(Vecs < DistCutSq, MaskNot23[MaskA])
        else:
          MaskB = Vecs < DistCutSq
        Rad = self.Radius[MaskA][MaskB] + self.Radius[i]
        Rad = Rad*Rad
        Eps = self.SqrtEps[MaskA][MaskB] * self.SqrtEps[i]
        if self.UseBonds: Eps[Mask4[MaskA][MaskB]] *= FFLJ14Scale
        Vecs = Rad / Vecs.compress(MaskB)
        Vecs = Vecs*Vecs*Vecs
        E += sum(Eps * Vecs*(Vecs - 2.))
      return E

  def ChargeScore(self, ResNum = None):
    """Returns an energy based on electrostatic interactions.
    If supplied, ResNum will return the contribution from atoms in
    residue number ResNum with themselves and other residues."""
    if self.UseBonds: self.__PrepBonds()
    if USELIB:
      if self.UseBonds:
        if ResNum is None:
          return CoulFact * proteinlib.chargescorebond(self.Pos, self.Charge, ChargeDistCut,
                                                       self.__Bonds23, self.__Bonds4, FFCharge14Scale)
        else:
          return CoulFact * proteinlib.chargescoreresbond(self.Pos, self.Charge, ChargeDistCut,
                                                          ResNum, self.AtomResNum,
                                                          self.__Bonds23, self.__Bonds4, FFCharge14Scale)
      else:
        if ResNum is None:
          return CoulFact * proteinlib.chargescore(self.Pos, self.Charge, ChargeDistCut)
        else:
          return CoulFact * proteinlib.chargescoreres(self.Pos, self.Charge, ChargeDistCut,
                                                      ResNum, self.AtomResNum)
    else:
      N = len(self.Pos)
      E = 0.
      DistCutSq = ChargeDistCut**2
      if ResNum is None:
        AtomRange = range(0, N-1)
      else:
        r = self.Res[ResNum]
        a1, a2 = r.StartAtom, r.StopAtom
        AtomRange = range(a1, a2)
      MaskA = ones(N,bool)
      for i in AtomRange:
        MaskA[i] = False
        Vecs = self.Pos[MaskA] - self.Pos[i,:]
        Vecs = (Vecs*Vecs).sum(axis=1)
        if self.UseBonds:
          MaskNot23, Mask4 = pfunc.GetBondMasks(i, N, self.__Bonds23, self.__Bonds4)
          MaskB = logical_and(Vecs < DistCutSq, MaskNot23[MaskA])
        else:
          MaskB = Vecs < DistCutSq
        Vecs = 1. / sqrt(Vecs.compress(MaskB))
        if self.UseBonds: Vecs[Mask4[MaskA][MaskB]] *= FFCharge14Scale
        E += self.Charge[i] * sum(self.Charge[MaskA][MaskB] * Vecs)
      return CoulFact * E    

  def Energy(self, ResNum = None):
    """Returns an energy based on LJ + electrostatic interactions.
    If supplied, ResNum will return the contribution from atoms in
    residue number ResNum with themselves and other residues."""
    if self.UseBonds: self.__PrepBonds()
    if USELIB:
      if self.UseBonds:
        if ResNum is None:
          return proteinlib.energybond(self.Pos, self.Radius, self.SqrtEps, self.Charge,
                                       self.LJDistCut, ChargeDistCut, CoulFact,
                                       self.__Bonds23, self.__Bonds4,
                                       FFLJ14Scale, FFCharge14Scale)
        else:
          return proteinlib.energyresbond(self.Pos, self.Radius, self.SqrtEps,
                                          self.Charge, self.LJDistCut, ChargeDistCut,
                                          CoulFact, ResNum, self.AtomResNum,
                                          self.__Bonds23, self.__Bonds4,
                                          FFLJ14Scale, FFCharge14Scale)
      else:
        if ResNum is None:
          return proteinlib.energy(self.Pos, self.Radius, self.SqrtEps, self.Charge,
                                   self.LJDistCut, ChargeDistCut, CoulFact)
        else:
          return proteinlib.energyres(self.Pos, self.Radius, self.SqrtEps,
                                      self.Charge, self.LJDistCut, ChargeDistCut,
                                      CoulFact, ResNum, self.AtomResNum)
    else:
      N = len(self.Pos)
      E = 0.
      DistCutSq = max(self.LJDistCut, ChargeDistCut)**2
      if ResNum is None:
        AtomRange = range(0, N-1)
      else:
        r = self.Res[ResNum]
        a1, a2 = r.StartAtom, r.StopAtom
        AtomRange = range(a1, a2)
      MaskA = ones(N,bool)
      for i in AtomRange:
        MaskA[i] = False
        Vecs = self.Pos[MaskA] - self.Pos[i,:]
        Vecs = (Vecs*Vecs).sum(axis=1)
        if self.UseBonds:
          MaskNot23, Mask4 = pfunc.GetBondMasks(i, N, self.__Bonds23, self.__Bonds4)
          MaskB = logical_and(Vecs < DistCutSq, MaskNot23[MaskA])
        else:
          MaskB = Vecs < DistCutSq
        Rad = self.Radius[MaskA][MaskB] + self.Radius[i]
        Rad = Rad*Rad
        Eps = self.SqrtEps[MaskA][MaskB] * self.SqrtEps[i]
        Vecs = Vecs.compress(MaskB)
        InvDist = 1. / sqrt(Vecs)
        if self.UseBonds:
          Eps[Mask4[MaskA][MaskB]] *= FFLJ14Scale
          InvDist[Mask4[MaskA][MaskB]] *= FFCharge14Scale
        E += CoulFact * self.Charge[i] * sum(self.Charge[MaskA][MaskB] * InvDist)
        Vecs = Rad / Vecs
        Vecs = Vecs*Vecs*Vecs
        E += sum(Eps * Vecs*(Vecs - 2.))        
      return E    

  def ResContactScore(self, ResID, EpsMat, Sigma = 6.):
    """Returns an energy based on the attractive portion of the LJ potential
    and allows residues to have specific interactions in a matrix.
    ResID are indices {ind} assigned to each residue.  EpsMat[ind1,ind2]
    gives the energy scale between residues with the given indices.
    EpsMat must be symmetric."""
    if USELIB:
      Pos = self.ResPos()
      return proteinlib.rescontactscore(Pos, ResID, EpsMat, Sigma)
    else:
      #get compressed positions
      Ind = where(ResID >= 0)[0]
      Pos = self.ResPos().take(Ind, axis=0)
      NewResID = ResID.take(Ind)
      N = len(Pos)
      E = 0.
      SigmaSq = Sigma*Sigma
      for i in range(0, N-1):
        Vecs = Pos[i+1:] - Pos[i,:]
        Vecs = SigmaSq / (Vecs*Vecs).sum(axis=1)
        Vecs = where(Vecs > 0.793700526, 0.793700526, Vecs)
        Vecs = Vecs*Vecs*Vecs
        Eps = EpsMat[NewResID[i],].take(NewResID[i+1:])
        E += sum(Eps*(Vecs - Vecs*Vecs))
      return 4.*E

  def HBondLJScore(self, ResInd = None, CosPower = 0, MinCO = 3):
    """Returns a dimensionless LJ 12-10 backbone hydrogen bonding energy.
    If supplied, ResInd will only use the specified residue numbers."""
    if USELIB:
      if ResInd is None:
        return proteinlib.hbondljscore(self.Pos,
               self.AtomList("N"), self.AtomList("H"), self.AtomList("O"),
               HBondLJA, HBondLJB, HBondLJDistCut, CosPower, MinCO)
      else:
        return proteinlib.hbondljscoreind(self.Pos,
               self.AtomList("N"), self.AtomList("H"), self.AtomList("O"),
               HBondLJA, HBondLJB, HBondLJDistCut, CosPower, MinCO, ResInd)
    else:
      E = 0.
      DistCutSq = HBondLJDistCut**2
      HInd, NInd, NHResInd = pfunc.GetCommonInd(self.AtomList("H"), self.AtomList("N"), ResInd)
      OInd, OResInd = pfunc.GetInd(self.AtomList("O"), ResInd)
      PosH = self.Pos.take(HInd, axis=0)
      PosO = self.Pos.take(OInd, axis=0)
      VecHN = self.Pos.take(NInd, axis=0) - PosH
      VecHNInvSq = 1./sum(VecHN*VecHN, axis=1)
      for i in range(0, len(PosO)):
        VecOH = PosH - PosO[i]
        Vecs1 = (VecOH*VecOH).sum(axis=1)
        Mask = logical_and(Vecs1 < DistCutSq, abs(OResInd[i] - NHResInd) >= MinCO)
        Vecs1 = 1. / Vecs1.compress(Mask, axis=0)
        Vecs3 = Vecs1*Vecs1*Vecs1
        if CosPower == 0:
          E += ((HBondLJA * Vecs3 - HBondLJB * Vecs1*Vecs1) * Vecs3).sum()
        else:
          Vecs2 = sum(VecHN[Mask] * VecOH[Mask], axis=1)
          Vecs2[Vecs2 < 0.] = 0.
          Vecs2 = Vecs2*Vecs2 * VecHNInvSq[Mask] * Vecs1
          if CosPower != 2: Vecs2 = Vecs2**(CosPower/2.)
          E += ((HBondLJA * Vecs3 - HBondLJB * Vecs1*Vecs1) * Vecs3 * Vecs2).sum()
      return E

  def HBondDipoleScore(self, ResInd = None, MinCO = 3):
    """Returns a dimensionless dipole-dipole backbone hydrogen bonding energy.
    If supplied, ResInd will only use the specified residue numbers."""
    if USELIB:
      if ResInd is None:
        return proteinlib.hbonddipolescore(self.Pos,
               self.AtomList("N"), self.AtomList("H"),
               self.AtomList("C"), self.AtomList("O"),
               HBondDipoleDistCut, MinCO)
      else:
        return proteinlib.hbonddipolescoreind(self.Pos,
               self.AtomList("N"), self.AtomList("H"),
               self.AtomList("C"), self.AtomList("O"),
               HBondDipoleDistCut, MinCO, ResInd)
    else:
      E = 0.
      DistCutSq = HBondDipoleDistCut**2
      HInd, NInd, NHResInd = pfunc.GetCommonInd(self.AtomList("H"), self.AtomList("N"), ResInd)
      CInd, OInd, COResInd = pfunc.GetCommonInd(self.AtomList("C"), self.AtomList("O"), ResInd)
      PosH = self.Pos.take(HInd, axis=0)   
      PosO = self.Pos.take(OInd, axis=0)
      VecNH = PosH - self.Pos.take(NInd, axis=0)
      VecCO = PosO - self.Pos.take(CInd, axis=0)
      for i in range(0, len(PosO)):
        Den = PosH - PosO[i]
        Den = (Den*Den).sum(axis=1)
        Mask = logical_and(Den < DistCutSq, abs(COResInd[i] - NHResInd) >= MinCO)
        Den = Den.compress(Mask)
        Den = Den**(-1.5)
        Num = dot(VecNH.compress(Mask, axis=0), VecCO[i])
        E += (Num * Den).sum()
      return E

  def HBondChargeScore(self, ResInd = None, CosPower1 = 2, CosPower2 = 2,
                       Coef = 1., MinCO = 3):
    """Returns a dimensionless charge-based backbone hydrogen bonding energy.
    If supplied, ResInd will only use the specified residue numbers."""
    if USELIB:
      if ResInd is None:
        return HBondChargeScale * proteinlib.hbondchargescore(self.Pos,
              self.AtomList("N"), self.AtomList("H"),
              self.AtomList("C"), self.AtomList("O"),
              HBondChargeDistCut, CosPower1, CosPower2, Coef, MinCO)
      else:
        return HBondChargeScale * proteinlib.hbondchargescoreind(self.Pos,
              self.AtomList("N"), self.AtomList("H"),
              self.AtomList("C"), self.AtomList("O"),
              HBondChargeDistCut, CosPower1, CosPower2, Coef, MinCO, ResInd)
    else:
      E = 0.
      DistCutSq = HBondChargeDistCut**2
      HInd, NInd, NHResInd = pfunc.GetCommonInd(self.AtomList("H"), self.AtomList("N"), ResInd)
      CInd, OInd, COResInd = pfunc.GetCommonInd(self.AtomList("C"), self.AtomList("O"), ResInd)
      PosH = self.Pos.take(HInd, axis=0)   
      PosO = self.Pos.take(OInd, axis=0)
      PosC = self.Pos.take(CInd, axis=0)   
      PosN = self.Pos.take(NInd, axis=0)
      VecHN = (PosN - PosH)
      DistHN = sqrt((VecHN*VecHN).sum(axis=1))
      VecHN = VecHN / DistHN[:,newaxis]
      VecCO = (PosO - PosC)
      DistCO = sqrt((VecCO*VecCO).sum(axis=1))
      VecCO = VecCO / DistCO[:,newaxis]
      for i in range(0, len(PosO)):
        VecOH = PosH - PosO[i]
        DistOH = (VecOH*VecOH).sum(axis=1)
        Mask = logical_and(DistOH < DistCutSq, abs(COResInd[i] - NHResInd) >= MinCO)
        DistOH = DistOH.compress(Mask)**(-0.5)
        VecOH = VecOH[Mask] * DistOH[:,newaxis]
        DistON = PosN[Mask] - PosO[i]
        DistCH = PosH[Mask] - PosC[i]
        DistCN = PosN[Mask] - PosC[i]
        DistON = (DistON*DistON).sum(axis=1)**(-0.5)
        DistCH = (DistCH*DistCH).sum(axis=1)**(-0.5)
        DistCN = (DistCN*DistCN).sum(axis=1)**(-0.5)
        CosVec1 = (VecOH * VecHN[Mask]).sum(axis=1)
        CosVec2 = (VecOH * VecCO[i]).sum(axis=1)
        CosVec1[CosVec1 < 0.] = 0.
        CosVec2[CosVec2 < 0.] = 0.
        CosVec = (CosVec1**CosPower1) * (CosVec2**CosPower2)
        E += ((Coef*DistON + Coef*DistCH - DistOH - Coef*DistCN) * CosVec).sum()
      return HBondChargeScale * E    

  def RadiusOfGyration(self, ResInd = None, ResPos = None):
    """Calculates residue-based radius of gyration."""
    #get residue positions
    if USELIB:
      #get residue positions
      if ResPos is None:
        Pos = self.ResPos(ResInd = ResInd)
        return proteinlib.rg(Pos)
      elif ResInd is None:
        return proteinlib.rg(ResPos)
      else:
        return proteinlib.rgatomind(ResPos, ResInd)      
    else:
      #get residue positions
      if ResPos is None:
        Pos = self.ResPos(ResInd = ResInd)
      elif ResInd is None:
        Pos = ResPos
      else:
        Pos = ResPos.take(ResInd, axis=0)
      #center positions
      Pos = Pos - Pos.sum(axis=0) / len(Pos)
      #calculate Rg
      return sqrt((Pos * Pos).sum() / len(Pos))


#======== PAIRWISE CONTACTS ========

  def ResContactList(self, ResInd = None, MinCO = None):
    """Returns a list of residue-residue contacts.
    ResInd is a list of residue indices to consider.
    MinCO is minimum contact order."""
    if ResInd is None: ResInd = range(0, len(self.Res))
    if MinCO is None: MinCO = self.MinCO
    #positions
    Pos = self.ResPos(ResInd = ResInd)
    #get the contact list
    return pfunc.GetContactList(Pos, ResInd, Radius = self.ResRadius, MinCO = MinCO)

  def ResContactMap(self, ResInd = None, MinCO = None):
    """Returns a list of residue-residue contacts.
    ResInd is a list of residue indices to consider.
    MinCO is minimum contact order."""
    if ResInd is None: ResInd = range(0, len(self.Res))
    if MinCO is None: MinCO = self.MinCO
    #positions
    Pos = self.ResPos(ResInd = ResInd)
    #get the contact list
    return pfunc.GetContactMap(Pos, ResInd, len(self.Res),
                         Radius = self.ResRadius, MinCO = MinCO)

  def PhobicContactList(self, MinCO = None):
    """Gives a list of hydrophobic-hydrophobic contacts in residue numbers.
    MinCO is minimum contact order."""
    ResInd = self.ResInd(Hydrophobic = True)
    return self.ResContactList(ResInd = ResInd, MinCO = MinCO)

  def NPhobicContact(self, MinCO = None):
    "Returns the number of hydrophobic contacts."
    return len(self.PhobicContactList(MinCO = MinCO))

  def SaltContactList(self, UseResInd = True):
    """Gives a list of salt bridge contacts in atom numbers.
    UseResInd = False will return contacts in terms of atom numbers
    rather than residue numbers."""
    #first get just those positions that are salt briges
    IndP = self.AtomInd(AtomCharged = True, ChargeSign = 1)
    IndN = self.AtomInd(AtomCharged = True, ChargeSign = -1)
    PosP = self.Pos.take(IndP, axis=0)
    PosN = self.Pos.take(IndN, axis=0)
    l = pfunc.GetContactList(PosP, IndP, PosN, IndN, Radius = self.SaltAtomRadius)
    if UseResInd:
      l = [(self.AtomResNum[a], self.AtomResNum[b]) for (a,b) in l]
      #remove duplicates
      l2 = []
      for v in l:
        if not v in l2: l2.append(v)
      return l2
    else:
      return l

  def NSaltContact(self):
    "Gives the number of salt bridge contacts."
    return len(self.SaltContactList(UseResInd = False))
  
  def ResCoordination(self, ResInd = None, MinCO = 1):
    """Returns the coordination number of each residue from all residues
    to those specified in ResInd (defaulting to all)."""
    N = len(self.Res)
    CM = self.ResContactMap(MinCO = MinCO)
    if ResInd is None:
      RC = CM.sum(axis=1)
    else:
      RC = zeros(N, float)
      #now filter along one direction
      for i in range(0,N):
        RC[i] = CM[i,:].take(ResInd).sum()
    return RC

  def MeanCoordination(self, FromResInd = None, ToResInd = None):
    """Computes the mean coordination number from central residues in
    FromResInd to other residues in ToResInd."""
    if FromResInd is None:
      return mean(self.ResCoordination(ResInd = ToResInd))
    else:
      return mean(self.ResCoordination(ResInd = ToResInd).take(FromResInd))

  def ResPartitioning(self, IntMinCoord = IntMinCoordDflt, ResInd = None):
    """Computes the number of exterior and interior residues.
    Returns NExteriorRes, NInteriorRes.
    IntMinCoord = minimum coordination number for internal residues."""
    if ResInd is None: ResInd = range(0, len(self.Res))
    N = len(ResInd)
    RC = self.ResCoordination().take(ResInd)
    NIntRes = where(RC >= IntMinCoord, 1, 0).sum()
    NExtRes = N - NIntRes
    return NExtRes, NIntRes

  def ResFracInterior(self, IntMinCoord = IntMinCoordDflt, ResInd = None):
    """Returns the fraction of residues in the interior.
    IntMinCoord = minimum coordination number for internal residues."""
    NExt, NInt = self.ResPartitioning(IntMinCoord = IntMinCoord, ResInd = ResInd)
    return float(NInt) / float(max(NExt + NInt, 1))

  def ResContactOrder(self, ResInd = None):
    """Returns a contact order matrix of the protein using residues in
    ResInd (defaulting to all)."""
    N = len(self.Res)
    ContactMap = self.ResContactMap(ResInd = ResInd)
    Order = abs(subtract.outer(arange(N), arange(N)))
    return where(ContactMap > 0, Order, ContactMap)

  def MeanContactOrder(self, ResInd = None):
    """Returns the mean contact order of the protein using
    residues in ResInd (defaulting to all)."""
    CO = self.ResContactOrder(ResInd = ResInd)
    N = where(CO > 0, 1, 0).sum()
    return CO.sum() / max(N, 1)
    


#======== HYDROGEN BONDS ========

  def HBondContactMap(self, ResInd = None):
    """Returns a contact map of backbone hydrogen bonding;
       the first residue (index) is for N-H, the second for O.
       ResInd is a list of residue indices to consider.
       Will use hydrogens to evaluate the bonds if there; otherwise
       the evaluation is made based on N--O-C."""
    N = len(self.Res)
    if ResInd is None: ResInd = range(0, N)
    NoH = all(self.AtomList("H") < 0)
    if USELIB:
      if NoH:
        return proteinlib.bondcontactmap(self.Pos,
          self.AtomList("N"), self.AtomList("O"), self.AtomList("C"),
          ResInd, HBondNODist, HBondNOCAng)  
      elif HBondDSSP:
        return proteinlib.hbondcontactmap(self.Pos,
          self.AtomList("N"), self.AtomList("H"), self.AtomList("C"),
          self.AtomList("O"), ResInd, HBondDSSPDistCut)
      else:
        return proteinlib.bondcontactmap(self.Pos,
          self.AtomList("O"), self.AtomList("H"), self.AtomList("N"),
          ResInd, HBondHODist, HBondNHOAng).transpose()
    else:
      CM = zeros((N, N), int)
      for i in ResInd:
        AtomN = self.AtomList("N")[i]
        AtomH = self.AtomList("H")[i]
        if AtomN < 0 or (AtomH < 0 and not NoH): continue
        for j in ResInd:
          if i == j: continue
          AtomO = self.AtomList("O")[j]
          AtomC = self.AtomList("C")[j]
          if AtomO < 0 or AtomC < 0: continue
          if NoH:
            v = pfunc.IsHBond(self.Pos[AtomN], self.Pos[AtomO], self.Pos[AtomC])
          else:
            v = pfunc.IsHBond(self.Pos[AtomN], self.Pos[AtomO],
                        self.Pos[AtomC], self.Pos[AtomH])
          CM[i,j] = int(v)
      return CM

  def HBondContactList(self, ResInd = None):
    """Returns a list of backbone hydrogen bonding contacts;
       the first residue (index) is for N-H, the second for O.
       ResInd is a list of residue indices to consider."""
    CM = self.HBondContactMap(ResInd)
    Ind1, Ind2 = CM.nonzero()
    return zip(Ind1,Ind2)

  def SecondaryStructure(self):
    """Returns a quick approximation to secondary structure based
on hydrogen bonding patterns and dihedral angles."""
    def InRange(Dih, DihRange):
      Phi, Psi = Dih
      MinPhi, MaxPhi, MinPsi, MaxPsi = DihRange
      return (Phi >= MinPhi and Phi <= MaxPhi and Psi >= MinPsi and Psi <= MaxPsi)
    def ListInRange(DihList, DihRangeList):
      for i in range(0, len(DihList)):
        if not InRange(DihList[i], DihRangeList[i]): return False
      return True
    def SameSS(SS, Label):
      for i in range(len(Label)):
        if not SS[i] == " " and not Label[i] == " " and not SS[i] == Label[i]:
          return False
      return True
    def UpdateSS(OldSS, NewSS):
      SS = ""
      for i in range(len(OldSS)):
        if OldSS[i] == " ":
          SS += NewSS[i]
        else:
          SS += OldSS[i]
      return SS
    N = len(self)
    Dih = [self.PhiPsi(i) for i in range(N)]
    Done = [False for i in range(N)]
    SS = " "*N
    CM = self.HBondContactMap()
    for (Label, NLocal, DihRangeList, CMTmpl) in SSTypes:
      #sort through nonlocal contacts
      NSeg = len(Label)
      NNon = NSeg - NLocal
      NContact = (CMTmpl > 0).sum()
      CMSeg = zeros((NSeg, NSeg), int)
      for i in range(N - NLocal + 1):
        if Done[i]: continue
        #check dihedrals
        if not ListInRange(Dih[i:i+NLocal], DihRangeList[:NLocal]):
          continue
        #start contact map
        CMSeg[:NLocal, :NLocal] = CM[i:i+NLocal, i:i+NLocal]
        #loop over nonlocal segments
        if NNon == 0:
          LoopN = 1
        else:
          LoopN = N - NNon + 1
        for j in range(LoopN):
          if not SameSS(SS[i:i+NLocal] + SS[j:j+NNon], Label): continue
          #make sure ranges are nonoverlapping
          if not(j + NNon - 1 < i or i + NLocal - 1 < j): continue
          #check nonlocal dihedrals
          if not ListInRange(Dih[j:j+NNon], DihRangeList[NLocal:]):
            continue
          #finish contact map
          CMSeg[:NLocal, NLocal:] = CM[i:i+NLocal, j:j+NNon]
          CMSeg[NLocal:, :NLocal] = CM[j:j+NNon, i:i+NLocal]
          CMSeg[NLocal:, NLocal:] = CM[j:j+NNon, j:j+NNon]
          #check the contact map
          if (CMSeg * CMTmpl).sum() < NContact: continue
          #label the ss types
          SS = SS[:i] + UpdateSS(SS[i:i+NLocal], Label[:NLocal]) + SS[i+NLocal:]
          SS = SS[:j] + UpdateSS(SS[j:j+NNon], Label[NLocal:]) + SS[j+NNon:]
          Done[i] = True
          break
    return SS




#======== TESTING ========

def TestLib(N = 25):
  "Runs comparison tests between the lib and the slower routines."
  global USELIB
  if not USELIB:
    print "Library currently not loaded."
    return
  p = ProteinClass(Seq = "ACDEFGHIKLMNPQRSTVWY")
  p.UseBonds = False
  for i in range(len(p)): p.RotateToPhiPsi(i, -60,-45)
  for i in range(0,2):
    if i == 0:
      USELIB = False
      print "Not using library..."
    else:
      USELIB = True
      print "Using library..."
    #Rg test with CB atoms
    p.ResAtom = "CB"
    StartTime = time.time()
    for i in range(N):
      x = p.RadiusOfGyration()
    print "RGCB   ", x, "%.2f sec" % (time.time() - StartTime)
    #Rg test with centroids
    p.ResAtom = "*"
    StartTime = time.time()
    for i in range(N):
      x = p.RadiusOfGyration()
    print "RGCEN  ", x, "%.2f sec" % (time.time() - StartTime)
    #Rg test with hydrophobics
    ResInd = p.ResInd(Hydrophobic = True)
    StartTime = time.time()
    for i in range(N):
      x = p.RadiusOfGyration(ResInd)
    print "RGPHO  ", x, "%.2f sec" % (time.time() - StartTime)
    #hbond test
    ResInd = range(0,len(p),2)
    StartTime = time.time()
    for i in range(N):
      x,y,z,w = p.HBondLJScore(), p.HBondLJScore(ResInd), \
                p.HBondLJScore(CosPower=1), p.HBondLJScore(CosPower=2)
    print "HBOND  ", x, y, z, w, "%.2f sec" % (time.time() - StartTime)
    #hbond dipole test
    ResInd = range(0,len(p),2)
    StartTime = time.time()
    for i in range(N):
      x,y = p.HBondDipoleScore(), p.HBondDipoleScore(ResInd)
    print "HBOND-D", x, y, "%.2f sec" % (time.time() - StartTime)
    #hbond charge test
    ResInd = range(0,len(p),2)
    StartTime = time.time()
    for i in range(N):
      x,y,z = p.HBondChargeScore(), p.HBondChargeScore(ResInd), \
              p.HBondChargeScore(CosPower1=0,CosPower2=0)
    print "HBOND-C", x, y, z, "%.2f sec" % (time.time() - StartTime)
    #steric test
    p.UseBonds = False
    StartTime = time.time()
    for i in range(N):
      x,y = p.StericScore(), p.StericScore(ResNum=0)
    print "STERIC ", x, y, "%.2f sec" % (time.time() - StartTime)
    #steric test with bonds
    #prep for bonds by running score once
    p.UseBonds = True
    x = p.StericScore()
    StartTime = time.time()
    for i in range(N):
      x,y = p.StericScore(), p.StericScore(ResNum=0)
    print "STERICB", x, y, "%.2f sec" % (time.time() - StartTime)
    #LJ test without bonds
    p.UseBonds = False
    StartTime = time.time()
    for i in range(N):
      x,y = p.LJScore(), p.LJScore(ResNum = 0)
    print "LJ     ", x, y, "%.2f sec" % (time.time() - StartTime)
    #LJ test with bonds
    StartTime = time.time()
    p.UseBonds = True
    for i in range(N):
      x,y = p.LJScore(), p.LJScore(ResNum=0)
    print "LJB    ", x, y, "%.2f sec" % (time.time() - StartTime)
    #charge test
    p.UseBonds = False
    StartTime = time.time()
    for i in range(N):
      x,y = p.ChargeScore(), p.ChargeScore(ResNum=0)
    print "CHARGE ", x, y, "%.2f sec" % (time.time() - StartTime)
    #charge test
    p.UseBonds = True
    StartTime = time.time()
    for i in range(N):
      x,y = p.ChargeScore(), p.ChargeScore(ResNum=0)
    print "CHARGEB", x, y, "%.2f sec" % (time.time() - StartTime)
    #energy test
    p.UseBonds = False
    StartTime = time.time()
    for i in range(N):
      x,y = p.Energy(), p.Energy(ResNum=0)
    print "ENERGY ", x, y, "%.2f sec" % (time.time() - StartTime)
    #energy test
    p.UseBonds = True
    StartTime = time.time()
    for i in range(N):
      x,y = p.Energy(), p.Energy(ResNum=0)
    print "ENERGYB", x, y, "%.2f sec" % (time.time() - StartTime)
    #test of residue contacts
    p.UseBonds = False
    EpsMat = array([[0., 0.5], [0.5, 1.]], float)
    ResID = p.ResInd(Hydrophobic=True, ReturnMask=True).astype(int)
    StartTime = time.time()
    p.ResAtom = "*"
    for i in range(N):
      x = p.ResContactScore(ResID, EpsMat)
    print "RESCONT", x, "%.2f sec" % (time.time() - StartTime)
   

def Profile():
  import cProfile, __main__
  __main__.profilep = ProteinClass("ACDEFGHIKLMNPQRSTVWY")
  cProfile.run("profilep.RunMC(lambda q, rn: q.LJScore(rn), 10000, profilep.GetMCMoves(), 1.0, 0.2)")
  del __main__.profilep


#======== TEMPLATE LOADING =========

def AddTemplates(PdbFile):
  """Adds residues in PdbFile to the template library."""
  global Templates
  if os.path.isfile(PdbFile):
    p = ProteinClass(Pdb = PdbFile)
    #put residues into the template
    for i in range(0, len(p)):
      r = p[i]
      Templates[r.Res[0].Name] = r
  else:
    raise IOError, "File %s not found." % PdbFile

 

#make a global dictionary of residue templates
Templates = {}
#search the paths for the template file
for p in sys.path:
  fn = os.path.join(p, TemplateFile)
  if os.path.isfile(fn):
    AddTemplates(fn)
