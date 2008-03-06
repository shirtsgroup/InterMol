#!/usr/bin/env python


#LAST MODIFIED: 03-01-07


Usage = """Assembles two pdb structures using a Monte Carlo algorithm.

Usage     : assemble.py PDB1 PDB2 LOOPSEQ PREFIX [OPTIONS]
   OR     : assemble.py -s SEQ SS PREFIX [OPTIONS]

PDB*      : input pdb files or lists of pdb files 
LOOPSEQ   : one-letter sequence codes
PREFIX    : prefix of output pdb files; will append 001, 002, etc
SEQ       : amino acid sequence
SS        : secondary structure sequence (can include spaces or 'X')
OPTIONS   : "--nkeep=N" to set maximum number of configs (dflt 50)
            "--ngen=N" to set number of configs to generate before sorting
            "--len=x" to scale the no of MC trials/conf by x (default 1)
"""

#check for instructions
import sys
if __name__ == "__main__" and len(sys.argv) < 3:
  print Usage
  sys.exit()


from numpy import *
import os, sys, math, shutil, glob, rmsd
import sequence, pdbtools, protein, proteinalg, geometry
import scripttools

#GLOBALS:
DEBUG = False

#energy function weights
RgAWeight = 0.
RgHWeight = 0.
StericWeight = 0.01
PhobicContactWeight = 1.0
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

#dictionary of PhiPsi angles
SSDict = {"E":(-135., 135), "H":(-58., -47.)}



def MakeScoreFn1(p, LoopRes):
  #calculate variables
  ResID = -ones(len(p), int)
  ResID[p.ResInd(Hydrophobic = True)] = 0
  EpsMat = -ones((1,1), float)
  HBondResInd = array([x for x in range(len(p)) if not x in LoopRes], int)
  RgHResInd = array([x for x in p.ResInd(Hydrophobic = True) if not x in LoopRes], int)
  RgAResInd = array([x for x in range(len(p)) if not x in LoopRes], int)
  def ScoreFn(q,rn):
    if rn is None:
      #get the Rg terms
      if RgAWeight > 0 or RgHWeight > 0:
        ResPos = q.ResPos()
        RgA = q.RadiusOfGyration(ResInd = RgAResInd, ResPos = ResPos)
        RgH = q.RadiusOfGyration(ResInd = RgHResInd, ResPos = ResPos)
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
      HBondScore = q.HBondChargeScore(ResInd = HBondResInd, Coef = HBondCoef,
                                      CosPower1=HBondCosPower, CosPower2=HBondCosPower)
      return   PhobicContactWeight * PhobicScore \
             + StericWeight * StericScore \
             + HBondWeight * HBondScore \
             + RgAWeight * RgA + RgHWeight * RgH
    else:
      return StericWeight * q.StericScore(rn)
  return ScoreFn

def MakeScoreFn2(p, LoopRes):
  #calculate variables
  ResID = -ones(len(p), int)
  ResID[p.ResInd(Hydrophobic = True)] = 0
  EpsMat = -ones((1,1), float)
  HBondResInd = array([x for x in range(len(p)) if not x in LoopRes], int)
  RgHResInd = array([x for x in p.ResInd(Hydrophobic = True) if not x in LoopRes], int)
  RgAResInd = array([x for x in range(len(p)) if not x in LoopRes], int)
  def ScoreFn(q,rn):
    if rn is None:
      #get the Rg terms
      if RgAWeight > 0 or RgHWeight > 0:
        ResPos = q.ResPos()
        RgA = q.RadiusOfGyration(ResInd = RgAResInd, ResPos = ResPos)
        RgH = q.RadiusOfGyration(ResInd = RgHResInd, ResPos = ResPos)
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
      HBondScore = q.HBondChargeScore(ResInd = HBondResInd, Coef = HBondCoef,
                                      CosPower1=HBondCosPower, CosPower2=HBondCosPower)
      return   PhobicContactWeight * PhobicScore \
             + LJWeight * LJScore \
             + HBondWeight * HBondScore \
             + RgAWeight * RgA + RgHWeight * RgH
    else:
      return LJWeight * q.LJScore(rn)
  return ScoreFn


def OptimizeMC(p, LoopRes, Prefix, NKeep,
  MCLen = 1., NGen = None, RMSDTol = 2.0,
  Verbose = False, Summary = True):

  if NGen is None: NGen = NKeep
  
  #set protein class options and prep for score
  if LowRes:
    psamp = p.CoarseGrained()
    psamp.ResAtom = "cen"
  else:
    psamp = p.Copy()
    psamp.ResAtom = "CB"
  ScoreFn = MakeScoreFn1(psamp, LoopRes)
  
  #run MC algorithm
  MCMoves = psamp.GetMCMoves(ResIndPhiPsi = LoopRes, ResIndChi = LoopRes,
                             WeightPhiPsi = 0.5, WeightChi = 0.5, WeightChain = 0.)
  if len(LoopRes) == 0:
    Pdb = "%s_0.pdb" % Prefix
    p.WritePdb(Pdb)
    Gen = [(ScoreFn(p, None), Pdb)]
  elif UseCSA:
    Gen = proteinalg.RunCSA(psamp, ScoreFn, MCMoves,
                            NGen, NKeep, max(NKeep, CSANBank), Prefix,
                            DistTol = RMSDTol, MCLen = MCLen,
                            VerboseMode = int(Verbose), Summary = Summary)
  else:
    Gen = proteinalg.RunBank(psamp, ScoreFn, MCMoves,
                             NGen, NKeep, Prefix,
                             DistTol = RMSDTol, MCLen = MCLen,
                             VerboseMode = int(Verbose), Summary = Summary)

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


def MultiProtList(Pdbs, LoopSeqs, SeqBefore = [], SeqAfter = []):
  """Generates combinations of pdb files as proteinclassess."""
  
  #make a list of lists of pdb filenames along the chain  
  PdbList = []
  for p in Pdbs:
    if type(p) is list:
      PdbList.append(p)
    else:
      PdbList.append([p])
  N = len(PdbList)
  NPdb = [len(p) for p in PdbList]
  CurInd = [0]*N

  #loop over all combinations of pdb files
  pList, LoopResList = [], []
  while CurInd[-1] < NPdb[-1]:
    ThisPdbs = [PdbList[i][CurInd[i]] for i in range(N)]
    
    #prep the pdb files and make the protein class
    p = protein.ProteinClass(Seq = SeqBefore)
    LoopRes = range(len(p))
    for i, Pdb in enumerate(ThisPdbs):
      p = p - protein.ProteinClass(Pdb = Pdb)
      if i < len(LoopSeqs):
        Seq = LoopSeqs[i]
        LoopRes.extend(range(len(p), len(p) + len(Seq)))
        p = p - protein.ProteinClass(Seq = Seq)
    LoopRes.extend(range(len(p), len(p) + len(SeqAfter)))
    p = p - protein.ProteinClass(Seq = SeqAfter)
    pList.append(p)
    LoopResList.append(LoopRes)
    
    #increment the indices
    CurInd[0] += 1
    for i in range(0, N-1):
      if CurInd[i] >= NPdb[i]:
        CurInd[i] = 0
        CurInd[i+1] += 1
        
  return pList, LoopResList


def ZipGen(Gen, NKeep):
  x = 1. / float(max(NKeep,1))
  return [(1. - x*(i + 0.5), y) for (i,y) in enumerate(Gen)]

def OptimizeMCMulti(pList, LoopResList, Prefix, NKeep,
  MCLen = 1, Weights = None, NGenEach = None, RMSDTol = 2.0,
  Verbose = False, Summary = True):
  """Optimizes multiple protein classes with different loops."""
  
  if NGenEach is None: NGenEach = NKeep
  if Weights is None: Weights = ones(len(pList), float)
  if Weights.sum() == 0: Weights = ones(len(pList), float)
  Weights = Weights / Weights.sum()

  #loop over all combinations of pdb files  
  Gen1, Gen2 = [], []
  for i in range(0, len(pList)):
    p = pList[i].Copy()
    LoopRes = LoopResList[i]
    ThisPrefix = Prefix + "-%d" % i
    #run the assembly
    ThisGen = OptimizeMC(p, LoopRes, ThisPrefix, NKeep,
                         MCLen = MCLen, NGen = NGenEach, RMSDTol = RMSDTol,
                         Verbose = Verbose, Summary = False)
    #estimate the number of structures to take from this one
    NKeepThis = int(Weights[i] * NKeep + 0.5)
    Gen1.extend(ZipGen(ThisGen[:NKeepThis], NKeepThis))
    Gen2.extend(ZipGen(ThisGen[NKeepThis:], NKeepThis))

  #now sort
  Gen1.sort(reverse = True)
  Gen2.sort(reverse = True)
  Gen1 = [l for (w,l) in Gen1]
  Gen2 = [l for (w,l) in Gen2]
  Gen = Gen1 + Gen2     

  #make a new bank
  Bank = proteinalg.BankClass(NKeep, Prefix = Prefix, DistTol = RMSDTol, UseScore = False)
  while len(Bank) < NKeep and len(Gen) > 0:
    (Score, Pdb) = Gen.pop(0)
    Bank.Add(Pdb, Score)
  if Summary: Bank.WriteSummary()

  #delete old files
  for (Score, Pdb) in Gen1 + Gen2:
    try:
      os.remove(Pdb)
    except IOError:
      print "Could not remove temporary file %s" % Pdb
  
  return Bank.GetList()


def MakeProt(Seq, SS):
  Seq = sequence.SeqToList(Seq)
  if not len(Seq) == len(SS):
    raise IndexError, "Seq len is %d but SS len is %d." % (len(Seq), len(SS))
  p = protein.ProteinClass(Seq = Seq)
  SS = SS.upper()
  LoopRes = []
  for i, s in enumerate(SS):
    if s in SSDict:
      Phi, Psi = SSDict[s]
      p.RotateToPhiPsi(i, Phi, Psi)
    else:
      LoopRes.append(i)
  #get rid of overlaps
  proteinalg.RemoveOverlaps(p, ResIndPhiPsi = LoopRes)
  return p, LoopRes

def OptimizeMCFromSS(Seq, SS, Prefix, NKeep,
  MCLen = 1., NGen = None, RMSDTol = 2.0,
  Verbose = False, Summary = True):
  Seq = sequence.SeqToList(Seq)
  if Verbose:
    print "Making initial protein and pruning for overlaps."
  p, LoopRes = MakeProt(Seq, SS)
  if Verbose:
    print "Loop residues: " + ", ".join(["%s%d" % (Seq[i], i+1) for i in LoopRes])
  OptimizeMC(p, LoopRes, Prefix, NKeep,
             MCLen = MCLen, NGen = NGen, RMSDTol = RMSDTol,
             Verbose = Verbose, Summary = Summary)


#command line running
if __name__ == "__main__":
  Args = scripttools.ParseArgs(sys.argv, {"nkeep":50, "len":1., "ngen":0})
  MCLen = Args["len"]
  NKeep, NGenEach = Args["nkeep"], Args["ngen"]
  if NGenEach == 0: NGenEach = None
  if "s" in Args["FLAGS"]:
    Seq, SS, Prefix = Args[1], Args[2], Args[3]
    OptimizeMCFromSS(Seq, SS, Prefix, NKeep, MCLen = MCLen, NGen = NGenEach,
                     Verbose = True, Summary = True)
  else:
    PdbFile1, PdbFile2 = Args[1], Args[2]
    if "," in PdbFile1: PdbFile1 = PdbFile1.replace("[","").replace("]","").split(",")
    if "," in PdbFile2: PdbFile2 = PdbFile2.replace("[","").replace("]","").split(",")
    LoopSeq = sequence.SeqToList(Args[3])
    Prefix = Args[4]
    pList, LoopResList = MultiProtList([PdbFile1, PdbFile2], [LoopSeq])
    OptimizeMCMulti(pList, LoopResList, Prefix, NKeep,
      MCLen = MCLen, NGenEach = NGenEach, Verbose = True, Summary = True)



