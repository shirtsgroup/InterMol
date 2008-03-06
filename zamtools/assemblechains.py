#!/usr/bin/env python

#LAST MODIFIED: 05-23-07

Usage = """Assembles two pdb chains using a Monte Carlo algorithm.

Usage     : assemblechains.py PDB1 PDB2 PREFIX [OPTIONS] [WEIGHTS]

*PDB      : input pdb files or lists of pdb files 
PREFIX    : prefix of output pdb files; will append 001, 002, etc
OPTIONS   : "--nkeep=N" to set maximum number of configs (dflt 50)
            "--ngen=N" to set number of configs to generate before sorting
            "--flexres1=1,4,5-10" to allow backbone flex in residues 1, 4, and 5-10 in pdb1
            "--flexres2=1,4,5-10" to allow backbone flex in residues 1, 4, and 5-10 in pdb2
            "--hires" to sample using all atoms
            "--len=x" to scale the # of MC trials by x (default 10)
            "--rmsd=x" to make the cutoff rmsd for diff. structures x (default 0.0)
            "--restr=12-34" requires residues 12 and 34 to be within 8.0 A one from another
WEIGHTS   : "--phob=x" to give hydrophobic residue pairs energy x (default 4.0)
            "--salt=x" to give ion pairing residues energy x (default 4.0)
            "--else=x" to give any other residue pairs energy x (default 2.0)
            "--rloc=x" to set the weight for a preferred residue-residue interaction (default 8.0)
"""

#check for instructions
import sys
if __name__ == "__main__" and len(sys.argv) < 3:
  print Usage
  sys.exit()

from numpy import *
from geometry import *
import os, sys, math, shutil, glob, rmsd
import sequence, pdbtools, protein, proteinalg, geometry
import scripttools


#energy function weights
RestraintWeight = 2.0
PhobicContactWeight = 4.0
SaltContactWeight = 4.0
OtherContactWeight = 2.0
HBondWeight = 0.
HBondCosPower = 2
HBondCoef = 1.
LJWeight = 1.0
RefWeight = 1.0
ResWeight = 8.0
SpaceBuffer = 15.

#aa restraint cutoff
ResCutoff = 8.0

#conformational space annealing?
UseCSA = False
CSANBank = 50


def MakeScoreFn(p, RestList = None):
  #calculate variables
  ResID = 3 * ones(len(p), int)
  ResID[p.ResInd(Hydrophobic = True)] = 0
  ResID[p.ResInd(ChargeSign = 1)] = 1
  ResID[p.ResInd(ChargeSign = -1)] = 2
  EpsMat = -ones((4,4), float) * OtherContactWeight
  EpsMat[0,0] = -PhobicContactWeight
  EpsMat[1,2] = -SaltContactWeight
  EpsMat[2,1] = EpsMat[1,2]
  ResInd1 = range(p.ChainResNums[0], p.ChainResNums[1])
  ResInd2 = range(p.ChainResNums[1], p.ChainResNums[2])
  Dist = p.RadiusOfGyration(ResInd = ResInd1) + p.RadiusOfGyration(ResInd = ResInd2) + SpaceBuffer
  if RestList is None: RestList = []
  def ScoreFn(q,rn):
    if rn is None:
      #apply any res-res restraints
      ResScore = 0.
      if len(RestList) > 0:
        Pos = q.ResPos()
        for (a,b) in RestList:
          ResScore += max(0.0, geometry.Length(Pos[a] - Pos[b]) - ResCutoff) 
      #restraint interaction (to keep the peptide from floating away)
      l = geometry.Length(q.Centroid(ChainNum=0) - q.Centroid(ChainNum=1))
      RestraintScore = max(0., l - Dist)**2
      #get hydrophobic interaction
      ContactScore = q.ResContactScore(ResID, EpsMat)
      #get the Sterics
      LJScore = q.LJScore()
      #get the hbond energy
      HBondScore = q.HBondChargeScore(Coef = HBondCoef,
                                      CosPower1=HBondCosPower,
                                      CosPower2=HBondCosPower)
      return   RestraintWeight * RestraintScore \
             + ContactScore \
             + LJWeight * LJScore \
             + HBondWeight * HBondScore \
             + ResWeight * ResScore

    else:
      return LJWeight * q.LJScore(rn)
  return ScoreFn


def MakeProtein(p1, p2):
  """Makes a combined proteinclass."""
  n1, n2 = len(p1), len(p2)
  p = p1.Copy().Concat(p2)
  p.NewChain(ResInd=range(n1,n1+n2))
  Shift = p1.Pos.max(axis=0) - p2.Pos.min(axis=0) + 5.
  p.Translate(Shift, ChainNum = 1)
  p.Center(ChainNum = 0)
  return p, n1, n2


def RunAssembly(p1, p2, Prefix, NKeep,
  MCLen = 1., NGen = None, RMSDTol = 4.0,
  FlexRes1 = [], FlexRes2 = [], LowRes = True, RestList = None,
  Verbose = False, Summary = True):

  if NGen is None: NGen = NKeep
  if RestList is None: RestList = []

  #make the new proteinclass with both proteins added
  p, n1, n2 = MakeProtein(p1, p2)
  p1Res = array(range(n1), int)
  p2Res = array(range(n2), int) + n1

  #get which residues to sample
  SampResSC = range(len(p))
  SampResBB = FlexRes1 + [i + n1 for i in FlexRes2]

  #comment on restraints
  if Verbose:
    for (a,b) in RestList:
      print "Restraining residue %s%d to residue %s%d with weight %f" % \
            (p.Res[a].Name, a+1, p.Res[b].Name, b+1, ResWeight)
  
  #set protein class options and prep for score
  if LowRes:
    psamp = p.CoarseGrained()
    psamp.ResAtom = "cen"
  else:
    psamp = p.Copy()
    psamp.ResAtom = "CB"
  ScoreFn = MakeScoreFn(psamp, RestList = RestList)

  #define a distance function for comparison
  def DistFn(Pos1, Pos2):
    return sqrt(sum((Pos1[p2Res] - Pos2[p2Res])**2) / n2)
  
  #run MC algorithm
  MCMoves = psamp.GetMCMoves(ResIndPhiPsi = SampResBB, ResIndChi = SampResSC,
                             WeightPhiPsi = 0.25, WeightChi = 0.25, WeightChain = 0.5)
##  s = sum([m.NMove for m in MCMoves]) * MCLen
##  proteinalg.MCExplore(psamp, ScoreFn, int(proteinalg.BankStepsExplore1 * s), MCMoves,
##                       KeepMin = False, SnapSteps = 25, SnapFile = "c:/stage1.pdb")
##  proteinalg.MCExplore(psamp, ScoreFn, int(proteinalg.BankStepsExplore2 * s), MCMoves,
##                       KeepMin = True, SnapSteps = 25, SnapFile = "c:/stage2.pdb")
##  proteinalg.MCRefine(psamp, ScoreFn, int(proteinalg.BankStepsRefine2 * s), MCMoves,
##                      KeepMin = True, SnapSteps = 25, SnapFile = "c:/stage3.pdb")
##  proteinalg.MCTweak(psamp, ScoreFn, int(proteinalg.BankStepsTweak * s), MCMoves,
##                     KeepMin = True, SnapSteps = 25, SnapFile = "c:/stage4.pdb")
##  proteinalg.MCMinimize(psamp, ScoreFn, int(proteinalg.BankStepsMin * s), MCMoves,
##                        SnapSteps = 25, SnapFile = "c:/stage5.pdb")
##  f = file("c:/stages.pdb","w")
##  for i in range(1,6):
##    fn = "c:/stage%d.pdb" % i
##    f.write(file(fn).read())
##    os.remove(fn)
##  f.close()
##  sys.exit()

  if UseCSA:
    Gen = proteinalg.RunCSA(psamp, ScoreFn, MCMoves,
                            NGen, NKeep, max(NKeep, CSANBank), Prefix,
                            DistFn = DistFn, DistTol = RMSDTol,
                            MCLen = MCLen,
                            VerboseMode = int(Verbose), Summary = Summary)
  else:
    Gen = proteinalg.RunBank(psamp, ScoreFn, MCMoves,
                             NGen, NKeep, Prefix, 
                             DistFn = DistFn, DistTol = RMSDTol,
                             MCLen = MCLen,
                             VerboseMode = int(Verbose), Summary = Summary)

  #optimize sidechains
  for (Score, Pdb) in Gen:
    if Verbose: print "Optimizing side chains for %s" % Pdb
    psamp = protein.ProteinClass(Pdb = Pdb)
    psamp.Template(TemplateProtein = p)
    psamp.ResAtom = "CB"
    ScoreFn = MakeScoreFn(psamp, RestList = RestList)
    proteinalg.ResampleSC(psamp, ScoreFn, ResIndSC = SampResSC,
                          ResIndBB = p2Res, MCLen = MCLen,
                          RefWeight = RefWeight, Verbose = False,
                          LocalOpt = not LowRes)
    psamp.WritePdb(Pdb)
    
  #return scores and filenames
  return Gen


#Argument parsing
def GetResList(Arg):
  if Arg is None or Arg == "": return []
  ResList = []
  for s in """'"[]""":
    Arg = Arg.replace(s,"")
  for l in Arg.split(","):
    if "-" in l:
      a, b = [int(x) for x in l.split("-")]
      ResList.extend(range(a-1, b))
    else:
      a = int(l)
      ResList.append(a - 1)
  ResList.sort()
  ResList = [x for (i,x) in enumerate(ResList) if not x in ResList[i+1:]]
  return ResList

#command line running
if __name__ == "__main__":
  Args = scripttools.ParseArgs(sys.argv[1:], {"nkeep":50, "len":10., "ngen":0,
                                              "phob":PhobicContactWeight,
                                              "salt":SaltContactWeight,
                                              "else":OtherContactWeight,
                                              "rmsd":0.0,
                                              "rloc":ResWeight,
                                              "flexres1":"", "flexres2":""})
  PdbFile1, PdbFile2 = Args[0], Args[1]
  Prefix = Args[2]
  MCLen = Args["len"]
  NKeep, NGen = Args["nkeep"], Args["ngen"]
  LowRes = not "hires" in Args["FLAGS"]
  FlexRes1 = GetResList(Args["flexres1"])
  FlexRes2 = GetResList(Args["flexres2"])
  if NGen == 0: NGen = None
  PhobicContactWeight = Args["phob"]
  SaltContactWeight = Args["salt"]
  OtherContactWeight = Args["else"]
  RMSDTol = Args["rmsd"]
  ResWeight = Args["rloc"]
  if "restr" in Args:
    a, b = [int(x) for x in Args["restr"].split("-")]
    RestList = [(a-1, b-1)]
  else:
    RestList = []
  p1 = protein.ProteinClass(Pdb = PdbFile1)
  p2 = protein.ProteinClass(Pdb = PdbFile2)
  p1.Template()
  p2.Template()
  RunAssembly(p1, p2, Prefix, NKeep, MCLen = MCLen, NGen = NGen, RMSDTol = RMSDTol,
    FlexRes1 = FlexRes1, FlexRes2 = FlexRes2, LowRes = LowRes, RestList = RestList,
    Verbose = True, Summary = True)
