#!/usr/bin/env python

#LAST MODIFIED: 05-02-07

import os, sys, math, copy
from numpy import *
from protein import *


#CONSTANTS

#minimization parameters
MCMinT = 1.e-10
MCMinStepsUpdateMax = 100

#tweaking parameters
MCTwkT = 0.1
MCTwkStepsUpdateTemp = 100

#refinement parameters
MCRefT = 1.0
MCRefStepsUpdateTemp = 100

#exploration parameters
MCExpT = 5.0
MCExpStepsUpdateTemp = 100

#overlap algorithm
OvrTi = 1.e12
OvrTf = 1.e-8
OvrSteps = 50
OvrFracPhiPsi = 0.5
OvrStepsUpdateMax = 100

#side chain optimizing algorithm
SmpSCSteps1 = 50
SmpSCSteps2 = 10
SmpSCSteps3 = 10
SmpSCT1i = 1.e12
SmpSCT1f = 1.e-8
SmpSCT2 = 0.01
SmpSCT3 = 1.e-10
SmpSCStepsUpdateTemp = 100
SmpSCStepsUpdateMax = 100
SmpSCFracPhiPsi = 0.1
SmpSCFracChain = 0.1

#bank algorithm
BankStepsExplore1 = 50
BankStepsExplore2 = 150
BankStepsRefine1 = 150
BankStepsRefine2 = 250
BankStepsTweak = 50
BankStepsMin = 30

#csa algorithm
CSASwapNRes = 5



def RMSD(Pos1, Pos2):
  #center
  p1 = Pos1 - mean(Pos1, axis=0)
  p2 = Pos2 - mean(Pos2, axis=0)
  d1, d2 = shape(p1)
  #calculate E0
  E0 = sum(p1*p1, axis=None) + sum(p2*p2, axis=None)
  #calculate correlation matrix
  C = dot(transpose(p2), p1)
  #get singular value decomp
  V, S, Wt = linalg.svd(C)
  #check determinants
  Vdet = linalg.det(V)
  Wtdet = linalg.det(Wt)
  #if it's a reflection, reflect along lowest eigenvalue
  if Vdet*Wtdet < 0.: S[-1] = -S[-1]
  #compute rmsd
  rmsdsq = (E0 - 2. * sum(S)) / float(d1)
  rmsdsq = max(rmsdsq, 0.)
  return sqrt(rmsdsq)



class BankItemClass:
  def __init__(self, p, Score):
    self.Score = Score
    self.Pos = p.ResPos()
    self.DihList = [p.PhiPsi(i) for i in range(len(p))]

class BankClass:
  def __init__(self, N, Prefix = "pep", DistFn = RMSD, DistTol = 2.,
               NDigit = 3, UseScore = True):
    "Initializes a new bank class."
    self.Items = []
    self.N = N
    self.DistTol = DistTol
    self.DistFn = DistFn
    self.Prefix = Prefix
    self.NDigit = max(NDigit, int(math.log(N, 10) + 1))
    self.MaxScore = -1.e200
    self.MaxInd = None
    self.UseScore = UseScore
  def __len__(self):
    return len(self.Items)
  def Pdb(self, ind):
    "Returns the name of the pdb file in the bank for a given index."
    return self.Prefix + "_%0*d.pdb" % (self.NDigit, ind)
  def __UpdateMax(self):
    self.MaxScore, self.MaxInd = max([(bij.Score, j) for (j, bij) in enumerate(self.Items)])
  def __GetInd(self, bi):
    """Returns the index for placing a new bankitemclass."""
    #find the closest
    ind = -1
    if not self.DistFn is None and not self.DistTol is None:
      for (j, bij) in enumerate(self.Items):
        if self.DistFn(bi.Pos, bij.Pos) < self.DistTol:
          ind = j
          break
    #check the results
    if ind == -1:
      if len(self.Items) < self.N:
        #add to the end of the list
        ind = len(self.Items)
      elif self.UseScore and bi.Score < self.MaxScore:
        #replace the max score, if sorting
        ind = self.MaxInd
    else:
      if not self.UseScore:
        #there is already one like this
        ind = -1
      elif bi.Score > self.Items[ind].Score:
        #replaced one is better score
        ind = -1
    return ind
  def Add(self, p, Score):
    """Adds a pdb or proteinclass conformation to the bank.
    Returns the index if added, or -1 if not."""
    if not isinstance(p, ProteinClass): p = ProteinClass(Pdb = p)
    bi = BankItemClass(p, Score)
    #get index
    ind = self.__GetInd(bi)
    #put the item where it goes
    if ind < 0:
      return -1
    elif ind == len(self.Items):
      self.Items.append(bi)
    else:
      self.Items[ind] = bi
    #save the pdb
    bi.Pdb = self.Pdb(ind)
    p.WritePdb(bi.Pdb)
    self.__UpdateMax()
    return ind
  def Sort(self):
    "Sorts the conformations (and files) in the bank by score."
    if len(self.Items) == 0: return
    self.Items.sort(cmp = lambda x, y: cmp(x.Score, y.Score))
    #make a dictionary of the old to new file names
    Map = {}
    for (j, bij) in enumerate(self.Items):
      Map[bij.Pdb] = self.Pdb(j)
    #now go and move the files around
    oldf, olddat = None, None
    while len(Map) > 0:
      if oldf is None or not oldf in Map:
        oldf = Map.keys()[0]
        olddat = file(oldf,"r").read()
      newf = Map.pop(oldf)
      if newf == oldf:
        oldf, olddat = None, None
      else:
        newdat = file(newf,"r").read()
        file(newf,"w").write(olddat)
        oldf, olddat = newf, newdat
    #now update the pdb filenames
    for (j, bij) in enumerate(self.Items):
      bij.Pdb = self.Pdb(j)
  def Trim(self, N):
    "Keeps only the first N bank items."
    for bi in self.Items[N:]:
      os.remove(bi.Pdb)
    self.Items = self.Items[:N]
    self.__UpdateMax()
  def WriteSummary(self):
    "Writes a text summary of conformations in the bank."
    m = max([len(os.path.basename(bi.Pdb)) for bi in self.Items] + [8])
    s = "%-*s  %s" % (m, 'Filename', 'Score')
    for bi in self.Items:
      s += "\n%-*s  %.3f" % (m, os.path.basename(bi.Pdb), bi.Score)
    file(self.Prefix + "-summary.txt", "w").write(s)
  def GetList(self):
    "Returns a list of (Score,Pdb)."
    return [(bi.Score, bi.Pdb) for bi in self.Items]

    
def MakeStepFn(SnapSteps, SnapFile):
  if SnapFile is None or SnapSteps == 0:
    StepFn = None
  else:
    def StepFn(p, i, T, MCMoves):
      if (i+1) % SnapSteps == 0:
        file(SnapFile,"a").write(p.GetPdb())
  return StepFn


def MCMinimize(p, ScoreFn, NSteps, MCMoves, Verbose = False,
               SnapSteps = 100, SnapFile = None):
  """Finds a local minimum using MC moves."""
  for m in MCMoves: m.SetTweak()
  StepFn = MakeStepFn(SnapSteps, SnapFile)
  p.RunMC(ScoreFn, NSteps, MCMoves, Ti = MCMinT,
          StepsUpdateMax = MCMinStepsUpdateMax,
          StepFn = StepFn, Verbose = Verbose)  

def MCTweak(p, ScoreFn, NSteps, MCMoves, KeepMin = True, Verbose = False,
             SnapSteps = 100, SnapFile = None):
  """Attempts to do some very small adjustments around the current local minimum
  to find a structure of lower score."""
  for m in MCMoves: m.SetTweak()
  StepFn = MakeStepFn(SnapSteps, SnapFile)
  p.RunMC(ScoreFn, NSteps, MCMoves, Ti = MCTwkT,
          StepsUpdateTemp = MCTwkStepsUpdateTemp, 
          StepFn = StepFn, KeepMin = KeepMin, Verbose = Verbose)

def MCRefine(p, ScoreFn, NSteps, MCMoves, KeepMin = True, Verbose = False,
             SnapSteps = 100, SnapFile = None):
  """Attempts to do some localized searching around the current local minimum
  to find a structure of lowest score."""
  for m in MCMoves: m.SetRefine()
  StepFn = MakeStepFn(SnapSteps, SnapFile)
  p.RunMC(ScoreFn, NSteps, MCMoves, Ti = MCRefT, 
          StepsUpdateTemp = MCRefStepsUpdateTemp, 
          StepFn = StepFn, KeepMin = KeepMin, Verbose = Verbose)

def MCExplore(p, ScoreFn, NSteps, MCMoves, KeepMin = True, Verbose = False,
              SnapSteps = 100, SnapFile = None):
  """Searches broad portions of configuration space to find a low score."""
  for m in MCMoves: m.SetExplore()
  StepFn = MakeStepFn(SnapSteps, SnapFile)
  p.RunMC(ScoreFn, NSteps, MCMoves, Ti = MCExpT,
          StepsUpdateTemp = MCExpStepsUpdateTemp,
          StepFn = StepFn, KeepMin = KeepMin, Verbose = Verbose)


def RemoveOverlaps(p, ResIndPhiPsi = None, ResIndChi = None, MCLen = 1.,
                   Verbose = False, SnapSteps = 100, SnapFile = None):
  """Runs MC to remove overlaps by backbone and sidechain perturbations."""
  if ResIndPhiPsi is None: ResIndPhiPsi = range(len(p))
  if ResIndChi is None: ResIndChi = range(len(p))
  StepFn = MakeStepFn(SnapSteps, SnapFile)
  #define a scoring function
  def ScoreFn(q, rn):
    return q.LJScore(rn)
  MCMoves = p.GetMCMoves(ResIndPhiPsi = ResIndPhiPsi, ResIndChi = ResIndChi,
                         WeightPhiPsi = OvrFracPhiPsi,
                         WeightChi = 1. - OvrFracPhiPsi, WeightChain = 0.)
  MCLen = MCLen * sum([m.NMove for m in MCMoves])
  #run the optimizing
  for m in MCMoves: m.SetExplore()
  p.RunMC(ScoreFn, int(MCLen * OvrSteps), MCMoves, 
          Ti = OvrTi, Tf = OvrTf, StepsUpdateMax = OvrStepsUpdateMax,
          StepFn = StepFn, KeepMin = False, Verbose = Verbose)

  
def ResampleSC(p, ScoreFn, ResIndSC = None, ResIndBB = None,
               MCLen = 1., RefWeight = 0.1,
               Verbose = False, SnapSteps = 100, SnapFile = None,
               LocalOpt = False, MoveChains = True):
  """Runs MC to optimize sidechains for non-overlaps; allows some backbone flex."""
  if ResIndSC is None: ResIndSC = range(len(p))
  if ResIndBB is None: ResIndBB = ResIndSC
  if MoveChains:
    FracChain = SmpSCFracChain
  else:
    FracChain = 0.
  StepFn = MakeStepFn(SnapSteps, SnapFile)
  #get CA locations and positions
  Ind = p.AtomInd(AtomName = "CA")
  RefPos = p.Pos.take(Ind, axis=0)
  #define a new scoring function
  def NewScoreFn(q, rn):
    if rn is None:
      RefScore = sum((q.Pos.take(Ind, axis=0) - RefPos)**2, axis=None)
      return ScoreFn(q, rn) + RefWeight * RefScore
    else:
      return ScoreFn(q, rn)
  p.WritePdb("0.pdb")
  #first minimize the sidechains
  MCMoves = p.GetMCMoves(ResIndChi = ResIndSC, WeightPhiPsi = 0.,
                         WeightChi = 1., WeightChain = 0.)
  MCLen = MCLen * sum([m.NMove for m in MCMoves])
  if LocalOpt:
    for m in MCMoves: m.SetRefine()
    p.RunMC(ScoreFn, int(MCLen * SmpSCSteps1), MCMoves, 
            Ti = SmpSCT1f, Tf = SmpSCT1f, StepsUpdateMax = SmpSCStepsUpdateMax,
            StepFn = StepFn, KeepMin = False, Verbose = Verbose)
  else:
    for m in MCMoves: m.SetExplore()
    p.RunMC(ScoreFn, int(MCLen * SmpSCSteps1), MCMoves,
            Ti = SmpSCT1i, Tf = SmpSCT1f, StepsUpdateMax = SmpSCStepsUpdateMax, 
            StepFn = StepFn, KeepMin = False, Verbose = Verbose)
  #now sample the backbone and sidechains
  MCMoves = p.GetMCMoves(ResIndPhiPsi = ResIndBB, ResIndChi = ResIndSC,
                         WeightPhiPsi = SmpSCFracPhiPsi,
                         WeightChi = 1. - SmpSCFracPhiPsi - FracChain,
                         WeightChain = FracChain)
  for m in MCMoves: m.SetTweak()
  p.RunMC(NewScoreFn, int(MCLen * SmpSCSteps2), MCMoves,
          Ti = SmpSCT2, StepsUpdateTemp = SmpSCStepsUpdateTemp, 
          StepFn = StepFn, KeepMin = False, Verbose = Verbose)
  #now minimize everything
  p.RunMC(NewScoreFn, int(MCLen * SmpSCSteps3), MCMoves,
          Ti = SmpSCT3, StepsUpdateMax = SmpSCStepsUpdateMax, 
          StepFn = StepFn, KeepMin = False, Verbose = Verbose)
  

def RunBank(p, ScoreFn, MCMoves, NGen, NKeep, Prefix,
            DistFn = RMSD, DistTol = 2., MCLen = 1.,
            VerboseMode = 0, Summary = True,
            SnapSteps = 0):
  """Runs multiple annealing simulations to generate a bank of conformations."""
  Bank = BankClass(NKeep, Prefix = Prefix, DistFn = DistFn, DistTol = DistTol)
  MCLen = MCLen * sum([m.NMove for m in MCMoves])
  for i in range(NGen):
    SnapFile = "%s_a%0*d.pdb" % (Prefix, Bank.NDigit, i)
    if VerboseMode > 0: print "CONF %d/%d: exploration stage" % (i+1, NGen)
    MCExplore(p, ScoreFn, int(BankStepsExplore1 * MCLen),
              MCMoves, KeepMin = False, Verbose = VerboseMode > 1,
              SnapSteps = SnapSteps, SnapFile = SnapFile)
    MCExplore(p, ScoreFn, int(BankStepsExplore2 * MCLen),
              MCMoves, Verbose = VerboseMode > 1,
              SnapSteps = SnapSteps, SnapFile = SnapFile)
    if VerboseMode > 0: print "CONF %d/%d: refinement stage" % (i+1, NGen)
    MCRefine(p, ScoreFn, int(BankStepsRefine2 * MCLen),
             MCMoves, Verbose = VerboseMode > 1,
             SnapSteps = SnapSteps, SnapFile = SnapFile)
    if VerboseMode > 0: print "CONF %d/%d: tweaking stage" % (i+1, NGen)
    MCTweak(p, ScoreFn, int(BankStepsTweak * MCLen),
            MCMoves, Verbose = VerboseMode > 1,
            SnapSteps = SnapSteps, SnapFile = SnapFile)
    if VerboseMode > 0: print "CONF %d/%d: minimization stage" % (i+1, NGen)
    MCMinimize(p, ScoreFn, int(BankStepsMin * MCLen),
               MCMoves, Verbose = VerboseMode > 1,
               SnapSteps = SnapSteps, SnapFile = SnapFile)
    Bank.Add(p, ScoreFn(p,None))
  Bank.Sort()
  if Summary: Bank.WriteSummary()
  return Bank.GetList()


def RunCSA(p, ScoreFn, MCMoves, NGen, NKeep, NBank, Prefix, 
           DistFn = RMSD, DistTol = 2., MCLen = 1.,
           VerboseMode = 0, Summary = True,
           SnapSteps = 0): 
  if ResInd is None: ResInd = range(len(p))
  NKeep = min(NBank, NKeep)
  Bank = BankClass(NBank, Prefix = Prefix, DistFn = DistFn, DistTol = DistTol)
  MCLen = MCLen * sum([m.NMove for m in MCMoves])
  i = 0
  while i < NGen and len(Bank) < NBank:
    SnapFile = "%s_a%0*d.pdb" % (Prefix, Bank.NDigit, i)
    if VerboseMode > 0: print "CONF %d/%d: init exploration stage" % (i+1, NGen)
    MCExplore(p, ScoreFn, int(BankStepsExplore1 * MCLen),
              MCMoves, KeepMin = False, Verbose = VerboseMode > 1,
              SnapSteps = SnapSteps, SnapFile = SnapFile)
    MCExplore(p, ScoreFn, int(BankStepsExplore2 * MCLen),
              MCMoves, Verbose = VerboseMode > 1,
              SnapSteps = SnapSteps, SnapFile = SnapFile)
    if VerboseMode > 0: print "CONF %d/%d: init refinement stage" % (i+1, NGen)
    MCRefine(p, ScoreFn, int(BankStepsRefine2 * MCLen),
             MCMoves, Verbose = VerboseMode > 1,
             SnapSteps = SnapSteps, SnapFile = SnapFile)
    if VerboseMode > 0: print "CONF %d/%d: init tweaking stage" % (i+1, NGen)
    MCTweak(p, ScoreFn, int(BankStepsTweak * MCLen),
            MCMoves, Verbose = VerboseMode > 1,
            SnapSteps = SnapSteps, SnapFile = SnapFile)
    if VerboseMode > 0: print "CONF %d/%d: init minimization stage" % (i+1, NGen)
    MCMinimize(p, ScoreFn, int(BankStepsMin * MCLen),
              MCMoves, Verbose = VerboseMode > 1,
               SnapSteps = SnapSteps, SnapFile = SnapFile)
    Bank.Add(p, ScoreFn(p,None))
    i += 1
  #make a list of groups of residue numbers for swapping
  ResClust, Clust = [], []
  for j in range(len(p) + 1):
    if j in ResInd:
      Clust.append(j)
    elif len(Clust) > 0:
      ResClust.append(Clust)
      Clust = []
  ResGroups = []
  for Clust in ResClust:
    if len(Clust) < CSASwapNRes:
      ResGroups.append(Clust)
    else:
      for j in range(0, len(Clust) - CSASwapNRes + 1):
        ResGroups.append(Clust[j:j+CSASwapNRes])
  #run the bank hybridization
  while i < NGen:
    SnapFile = "%s_b%0*d.pdb" % (Prefix, Bank.NDigit, i)
    #pick two random bank structures
    n = len(Bank)
    x = int(random.random() * n)
    y = int(random.random() * n)
    while x==y:
      y = int(random.random() * n)
    bix, biy = Bank.Items[x], Bank.Items[y]
    #load the first positions
    p.Pos = ProteinClass(Pdb = bix.Pdb).Pos
    #pick a random residue stretch
    Group = ResGroups[int(random.random() * len(ResGroups))]
    #now swap the conformation
    for j in Group:
      Phi, Psi = biy.DihList[j]
      p.RotateToPhiPsi(j, Phi, Psi)
    #now run the MC
    if VerboseMode > 0: print "CONF %d/%d: making from bank %d and %d" % (i+1, NGen,x,y)
    if VerboseMode > 0: print "CONF %d/%d: bank refinement stage" % (i+1, NGen)
    MCRefine(p, ScoreFn, int(BankStepsRefine1 * MCLen),
             MCMoves, KeepMin = False, Verbose = VerboseMode > 1,
             SnapSteps = SnapSteps, SnapFile = SnapFile)
    MCRefine(p, ScoreFn, int(BankStepsRefine2 * MCLen),
             MCMoves, Verbose = VerboseMode > 1,
             SnapSteps = SnapSteps, SnapFile = SnapFile)
    if VerboseMode > 0: print "CONF %d/%d: bank tweaking stage" % (i+1, NGen)
    MCTweak(p, ScoreFn, int(BankStepsTweak * MCLen),
            MCMoves, Verbose = VerboseMode > 1,
            SnapSteps = SnapSteps, SnapFile = SnapFile)
    if VerboseMode > 0: print "CONF %d/%d: bank minimization stage" % (i+1, NGen)
    MCMinimize(p, ScoreFn, int(BankStepsMin * MCLen),
               MCMoves, Verbose = VerboseMode > 1,
               SnapSteps = SnapSteps, SnapFile = SnapFile)
    ret = Bank.Add(p, ScoreFn(p,None))
    if VerboseMode > 0:
      if ret < 0:
        if VerboseMode > 0: print "CONF %d/%d: BANK REJECTED" % (i+1, NGen)
      else:
        if VerboseMode > 0: print "CONF %d/%d: BANK ACCEPTED TO POSITION %d" % (i+1, NGen, ret)
    i += 1      
  Bank.Sort()
  Bank.Trim(NKeep)
  if Summary: Bank.WriteSummary()
  return Bank.GetList()






