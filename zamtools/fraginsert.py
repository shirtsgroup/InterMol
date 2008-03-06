#!/usr/bin/env python

#LAST MODIFIED: 05-01-07

Usage = """Uses a Monte Carlo algorithm with a library of fragments.

Usage     : fraginsert.py makelib PREFIX PDB1 PDB2 ... [OPTIONS]
   OR     : fraginsert.py run SEQ PREFIX FRAGLIB1 FRAGLIB2 ... [OPTIONS]
   OR     : fraginsert.py run SEQ PREFIX PDB1 PDB2 ... [OPTIONS]

SEQ       : amino acid sequence
PREFIX    : prefix of output pdb files; will append 001, 002, etc
FRAGLIB   : precompiled library file; ends in .dat
PDB*      : input pdb files for fragment 
OPTIONS   : "--nkeep=N" to set maximum number of configs (dflt 50)
            "--ngen=N" to set number of configs to generate before sorting
            "--len=x" to scale the no of MC trials/conf by x (default 1)
            "--ntrim=x" to trim frag pdbs by x residues on each end (default 0)
            "--minlen=x" to specify minimum fragment length (default 3)
            "--maxlen=x" to specify maximum fragment length (default unlimited)
            "--fragfrac=x" to specify fraction fragment moves (default 0.8)
            "--tol=x" to specify rmsd tolerance (default 2.0)
            "--weight=x" to force all fragment weights to x
"""

#check for instructions
import sys
if __name__ == "__main__" and len(sys.argv) < 3:
  print Usage
  sys.exit()


from numpy import *
import os, sys, shutil, glob, cPickle
import sequence, protein, proteinalg
import scripttools

#GLOBALS:
DEBUG = False

#energy function weights
RgAWeight = 0.
RgHWeight = 0.
StericWeight = 0.01
PhobicContactWeight = 10.0
HBondWeight = protein.HBondChargeRatioToSteric1 * StericWeight
HBondCosPower = 2
HBondCoef = 1.
LJWeight = 1. / protein.HBondChargeRatioToLJ1
RefWeight = 1.0

#mc variables
LowRes = True
VerboseMC = False
CSANBank = 20
UseCSA = False




def MakeScoreFn1(p):
  #calculate variables
  PhobicResInd = p.ResInd(Hydrophobic = True)
  ResID = -ones(len(p), int)
  ResID[PhobicResInd] = 0
  EpsMat = -ones((1,1), float)
  def ScoreFn(q,rn):
    if rn is None:
      #get the Rg terms
      if RgAWeight > 0 or RgHWeight > 0:
        ResPos = q.ResPos()
        RgA = q.RadiusOfGyration(ResPos = ResPos)
        RgH = q.RadiusOfGyration(ResInd = PhobicResInd, ResPos = ResPos)
      else:
        RgA, RgH = 0., 0.
      #get hydrophobic interaction
      if PhobicContactWeight > 0:
        PhobicScore = q.ResContactScore(ResID, EpsMat)
      else:
        PhobicScore = 0.
      #get the Sterics
      StericScore = q.StericScore()
      #get the hbond energy
      HBondScore = q.HBondChargeScore(Coef = HBondCoef,
                                      CosPower1=HBondCosPower, CosPower2=HBondCosPower)
      return   PhobicContactWeight * PhobicScore \
             + StericWeight * StericScore \
             + HBondWeight * HBondScore \
             + RgAWeight * RgA + RgHWeight * RgH
    else:
      return StericWeight * q.StericScore(rn)
  return ScoreFn

def MakeScoreFn2(p):
  #calculate variables
  PhobicResInd = p.ResInd(Hydrophobic = True)
  ResID = -ones(len(p), int)
  ResID[PhobicResInd] = 0
  EpsMat = -ones((1,1), float)
  def ScoreFn(q,rn):
    if rn is None:
      #get the Rg terms
      if RgAWeight > 0 or RgHWeight > 0:
        ResPos = q.ResPos()
        RgA = q.RadiusOfGyration(ResPos = ResPos)
        RgH = q.RadiusOfGyration(ResInd = PhobicResInd, ResPos = ResPos)
      else:
        RgA, RgH = 0., 0.
      #get hydrophobic interaction
      if PhobicContactWeight > 0:
        PhobicScore = q.ResContactScore(ResID, EpsMat)
      else:
        PhobicScore = 0.
      #get the Sterics
      LJScore = q.LJScore()
      #get the hbond energy
      HBondScore = q.HBondChargeScore(Coef = HBondCoef,
                                      CosPower1=HBondCosPower, CosPower2=HBondCosPower)
      return   PhobicContactWeight * PhobicScore \
             + LJWeight * LJScore \
             + HBondWeight * HBondScore \
             + RgAWeight * RgA + RgHWeight * RgH
    else:
      return LJWeight * q.LJScore(rn)
  return ScoreFn


def MakeLibFromPdbs(PdbList, MinLen = 3, MaxLen = 1000000, NResTrim = 0,
                    PctStr = r"_(\d+)pc", Weight = None, Verbose = False):
  """Extracts fragments from Lib relevant to p;
  Lib is a list of (Weight, Seq, DihList) tuples."""
  import re
  Lib = []
  for Pdb in PdbList:
    p = protein.ProteinClass(Pdb = Pdb)
    #use first chain only
    if len(p.Chains) > 1: p = p[:p.ChainResNums[1]]
    #trim if necessary
    if NResTrim > 0: p = p[NResTrim:-NResTrim]
    #check length
    if len(p) < MinLen or len(p) > MaxLen: continue
    #get the probabilities
    if Weight is None:
      ThisWeight = re.search(PctStr, Pdb)
      try:
        ThisWeight = max(1., float(ThisWeight.group(1))) / 100.
      except (AttributeError, IndexError, ValueError):
        ThisWeight = 1.
    else:
      ThisWeight = Weight
    if Verbose:
      print "Adding %s with weight %.2f" % (Pdb, ThisWeight)
    #add lib entry
    DihList = [p.PhiPsi(i) for i in range(len(p))]
    Lib.append((ThisWeight, p.Seq, DihList))
  return Lib

def MakeFragsFromLib(p, Lib, MinLen = 3, MaxLen = 1000000, Weight = None,
                     ResProb = None):
  """Extracts fragments from Lib relevant to p;
  Lib is a list of (Weight, Seq, DihList) tuples."""
  if ResProb is None: ResProb = ones(len(p), float)
  Frags = []
  for (ThisWeight, Seq, DihList) in Lib:
    if not Weight is None: ThisWeight = Weight
    Maps = sequence.GetAllMaps(p.Seq, Seq, MinLen = MinLen)
    for Map in Maps:
      ThisLen = min(MaxLen, len(Map))
      for i in range(0, len(Map) - ThisLen + 1):
        a, b = Map.a + i, Map.a + i + ThisLen
        c, d = Map.c + i, Map.c + i + ThisLen
        #check chain
        if not p.ResChain(a) == p.ResChain(b): continue
        #check probabilities
        if any(ResProb[a:b] <= 0.): continue
        ThisWeight *= average(ResProb[a:b])
        Frags.append((ThisWeight, a, DihList[c:d]))
  return Frags

def MakeFragsFromProt(p, pFrag, MinLen = 3, MaxLen = 1000000, Weight = 1.,
                      ResProb = None):
  """Extracts fragments from pFrag relevant to p."""
  if ResProb is None: ResProb = ones(len(p), float)
  Maps = GetAllMaps(p.Seq, pFrag.Seq, MinLen = MinLen)
  DihList = [pFrag.PhiPsi(i) for i in range(len(pFrag))]
  Frags = []
  for Map in Maps:
    ThisLen = min(MaxLen, len(Map))
    for i in range(0, len(Map) - ThisLen + 1):
      a, b = Map.a + i, Map.a + i + ThisLen
      c, d = Map.c + i, Map.c + i + ThisLen
      #check probabilities
      if any(ResProb[a:b] <= 0.): continue
      ThisWeight *= average(ResProb[a:b])
      Frags.append((Weight, a, DihList[c:d]))
  return Frags


def GetFragMoves(p, Frags, Weight = 1.0, CenterChain = 0):
  """Returns a list of MCMoveClass's for fragment insertion.
Frags is a list of (Weight, StartRes, [[Phi0,Psi0], [Phi1,Psi1], ...])
Weight is the probability weight of this move.
CenterChain is the chain number to center after each move; None for none."""
  #make the move classe
  MCMove = protein.pfunc.MCFragInsertClass(p, Frags, CenterChain)
  MCMove.W = Weight
  if MCMove.Active():
    return [MCMove]
  else:
    return []

  
def OptimizeMC(p, Frags, Prefix, NKeep,
  ResInd = None, FracFragMove = 0.8,
  MCLen = 1., NGen = None, RMSDTol = 2.0,
  Verbose = False, Summary = True):

  if NGen is None: NGen = NKeep
  if ResInd is None: ResInd = range(len(p))
  
  #set protein class options and prep for score
  if LowRes:
    psamp = p.CoarseGrained()
    psamp.ResAtom = "cen"
  else:
    psamp = p.Copy()
    psamp.ResAtom = "CB"
  ScoreFn = MakeScoreFn1(psamp)
  
  #run MC algorithm
  MCMoves = psamp.GetMCMoves(ResIndPhiPsi = ResInd, ResIndChi = ResInd,
                             WeightPhiPsi = (1. - FracFragMove) * 0.5,
                             WeightChi = (1. - FracFragMove) * 0.5)
  MCMoves.extend(GetFragMoves(p, Frags, Weight = FracFragMove))

  if UseCSA:
    Gen = proteinalg.RunCSA(psamp, ScoreFn, MCMoves,
                            NGen, NKeep, max(NKeep, CSANBank), Prefix,
                            DistTol = RMSDTol, MCLen = MCLen,
                            VerboseMode = int(Verbose), Summary = Summary)
  else:
    Gen = proteinalg.RunBank(psamp, ScoreFn, MCMoves,
                             NGen, NKeep, Prefix,
                             DistTol = RMSDTol, MCLen = MCLen,
                             VerboseMode = int(Verbose), Summary = Summary)

#ISSUES:
# 1) how to control mc scaling when fracfragsteps = 1. (since phi psi go to zero)
# 2) how to fix moves to allow rearranging from bad conformations

  #optimize sidechains
  for (Score, Pdb) in Gen:
    if Verbose: print "Optimizing side chains for %s" % Pdb
    psamp = protein.ProteinClass(Pdb = Pdb)
    if LowRes: psamp.Template()
    psamp.ResAtom = "CB"
    ScoreFn = MakeScoreFn2(psamp, LoopRes)
    proteinalg.ResampleSC(psamp, ScoreFn, MCLen = MCLen,
                          RefWeight = RefWeight, Verbose = False,
                          LocalOpt = not LowRes)
    psamp.WritePdb(Pdb)
    
  #return scores and filenames
  return Gen




#command line running
if __name__ == "__main__":
  Args = scripttools.ParseArgs(sys.argv, {"nkeep":50, "len":1.0, "ngen":0,
                                          "ntrim":2, "minlen":3, "maxlen":1000000,
                                          "tol":2.0, "fragfrac":0.8, "weight":-1.})
  MCLen = Args["len"]
  NKeep, NGen, NResTrim = Args["nkeep"], Args["ngen"], Args["ntrim"]
  if NGen == 0: NGen = NKeep
  MinLen, MaxLen = Args["minlen"], Args["maxlen"]
  RMSDTol, FracFragMove = Args["tol"], Args["fragfrac"]
  Weight = Args["weight"]
  if Weight < 0: Weight = None
  Cmd = Args[1].lower()
  if Cmd == "makelib":
    Prefix = Args[2]
    PdbList = [Args[x] for x in range(3, Args["NARG"])]
    print "Making fragment library"
    Lib = MakeLibFromPdbs(PdbList, NResTrim = NResTrim, Weight = Weight, Verbose = True)
    cPickle.dump(Lib, file(Prefix + ".dat", "w"))
    print "Saved library %s.dat" % Prefix
  elif Cmd == "run":
    Seq, Prefix = Args[2], Args[3]
    fl = [Args[x] for x in range(4, Args["NARG"])]
    PdbList = [x for x in fl if ".pdb" in x.lower()]
    LibList = [x for x in fl if ".dat" in x.lower()]
    Lib = []
    for LibName in LibList:
      print "Loading library %s" % LibName
      Lib.extend(cPickle.load(file(LibName, "r")))
    if len(PdbList) > 0:
      print "Adding PDBs to fragment library"
      l = MakeLibFromPdbs(PdbList, NResTrim = NResTrim, Weight = Weight, Verbose = True)
      Lib.extend(l)
    print "Making protein"
    p = protein.ProteinClass(Seq = Seq)
    p.OptimizeSC()
    print "Finding fragment alignments to protein sequence"
    Frags = MakeFragsFromLib(p, Lib, MinLen = MinLen, MaxLen = MaxLen)
    print "Running optimization"
    OptimizeMC(p, Frags, Prefix, NKeep, FracFragMove = FracFragMove,
               MCLen = MCLen, NGen = NGen, RMSDTol = RMSDTol,
               Verbose = True, Summary = True)




