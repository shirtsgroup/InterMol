#!/usr/bin/env python

#LAST MODIFIED: 05-03-07

Usage = """Splices together old pdb files to make a new sequence.

Usage     : makeconf.py SEQ PREFIX [OPTIONS] PDBFILE1 PDBFILE2 PDBFILE3...

SEQ       : one-letter sequence string
PREFIX    : new pdbfiles will be named PREFIX01.pdb, PREFIX02.pdb, etc.
OPTIONS   : "--cap" to cap files
            "--dehydrogen" to remove hydrogens 
            "--overlapok" to skip overlap detection
            "--max=N" to set the max number of configs generated to N
            "--MC" to use MC assembly optimization (time consuming)
            "--len=x" to scale the MC sim length by x (default 1.)
PDBFILE*  : old pdb files (alignment automatically performed)
"""

#check for instructions
import sys
if __name__ == "__main__" and len(sys.argv) == 1:
  print Usage
  sys.exit()


from numpy import *
import os, sys, copy
import pdbtools, sequence, scripttools, protein, assemble

#GLOBALS
MinResDflt = 3
MinResMap = 5
OverlapDist = 0.5
MCNGen = 50
MinLoopRes = 6

       

def RecurseAlign(StartRes, StopRes, SeqMaps, MinRes = MinResDflt):
  """Produces a list of lists of mappings of sequences in SeqMaps;
   StartRes is inclusive, StopRes is not."""
  #check that target sequence is zero length
  if StopRes <= StartRes:
    return [[]]
  #if target is too small, return non matched sequence
  if StopRes - StartRes < MinRes:
    return [[sequence.SeqMapClass(a = StartRes, b = StopRes)]]  
  #start with the base alignment of non matched
  AlignList = [[sequence.SeqMapClass(a = StartRes, b = StopRes)]]
  #check all sequence possibilities for alignments
  for i in range(0, len(SeqMaps)):
    Map = SeqMaps[i].SubMap(a = StartRes, b = StopRes)
    #if the alignment took...
    if len(Map) >= MinRes:
      #search for more
      LeftList = RecurseAlign(StartRes, Map.a, SeqMaps, MinRes)
      RightList = RecurseAlign(Map.b, StopRes, SeqMaps, MinRes)
      for l in LeftList:
        for r in RightList:
          AlignList.append(l + [Map] + r)
  #remove duplicates
  i = 0
  StrAlignList = [str(x) for x in AlignList]
  while i < len(AlignList) - 1:
    StrAlign = str(AlignList[i])
    if StrAlign in StrAlignList[i+1:]:
      ind = StrAlignList.index(StrAlign, i+1)
      del AlignList[ind]
      del StrAlignList[ind]
    else:
      i += 1
  return AlignList
  

def Score(Align, LogWeights, SeqLen, SeqsLen):
  MinP = min(LogWeights) - 10
  s = 0.
  for Map in Align:
    if Map.ID is None:
      s += MinP * len(Map) / float(SeqLen)
    else:
      s += LogWeights[Map.ID] * len(Map) / float(SeqsLen[Map.ID])
  return s

def LogWeight(Weight):
  if Weight <= 0.:
    return -1.e100
  else:
    return min(log(Weight), 0.)

def GetAlign(TarSeq, Seqs, SeqWeights = None, MultiMap = False):
  if SeqWeights is None or sum(SeqWeights) == 0:
    SeqWeights = ones(len(Seqs), float)
  else:
    SeqWeights = array(SeqWeights, float)
  SeqWeights = SeqWeights / sum(SeqWeights)
  LogWeights = [LogWeight(x) for x in SeqWeights]
  #make maps
  SeqMaps = []
  for i in range(len(Seqs)):
    if MultiMap:
      Maps = sequence.GetAllMaps(TargetSeq, Seqs[i], MinLen = MinResMap, ID = i)
      SeqMaps.extend(Maps)
    else:
      Map = sequence.SeqMapClass(Seq1 = TarSeq, Seq2 = Seqs[i], ID = i)
      if len(Map) > MinResMap: SeqMaps.append(Map)
  #perform alignments 
  AlignList = RecurseAlign(0, len(TarSeq), SeqMaps)
  #remove the fully unmatched alignments
  AlignList = AlignList[1:]
  #sort
  SeqsLen = [len(Seq) for Seq in Seqs]
  SeqLen = len(TarSeq)
  AlignList.sort(cmp = lambda x, y: cmp(Score(y, LogWeights, SeqLen, SeqsLen),
                                        Score(x, LogWeights, SeqLen, SeqsLen)))
  #get weights based on average priorities
  Weights = array([Score(Align, LogWeights, SeqLen, SeqsLen)
                   for Align in AlignList], float)
  Weights = exp(Weights - Weights.max())
  Weights = Weights / Weights.sum()
  return AlignList, Weights


def GetSeq(Pdb):
  """Returns the sequence for a PDB."""
  return sequence.Standardize(protein.ProteinClass(Pdb = Pdb).Seq)

def GetProt(Pdb):
  """Returns a proteinclass for a Pdb."""
  return protein.ProteinClass(Pdb = Pdb)

def MakePdbFromAlign(Align, TarSeq, Pdbs):
  p = protein.ProteinClass()
  for Map in Align:
    if Map.ID is None:
      p = p + protein.ProteinClass(Seq = Map.TrimSeq1(TarSeq))
    else:
      p = p + GetProt(Pdbs[Map.ID])[Map.c:Map.d]
  return p


def MakePdbAlignments(Seq, Pdbs, Prefix, MaxAlign = 15, MaxConf = None,
    Cap = True, CheckOverlap = True, Dehydrogen = False,
    PdbWeights = None, LogFile = None, MC = False,
    CapN = None, CapC = None, MCLen = 1., Verbose = False):
  """Makes pdb files for Seq using any overlaps with existing pdbs in Pdbs.
Generates new residues as needed and prioritizes for greatest overlap.
Seq: target sequence
Pdbs: pdb filenames or strings
Prefix: output pdb prefix
MaxAlign: maximum number of alignments to use
MaxConf: maximum number of pdb files to produce
Cap: True to add caps
Overlap: True to detect overlaps
Dehydrogen: true to dehydrogen
PdbWeights: positive weight numbers for each pdb file
LogFile: filename for log file
MC: True to use time-consuming MC optimization
MCLen: increase to do longer MC simulations
"""
  #check caps
  if CapN is None: CapN = Cap
  if CapC is None: CapC = Cap
  #check max conf
  if MaxConf is None: MaxConf = MaxAlign
  #make labels
  Labels = []
  for p in Pdbs:
    if "\n" in p:
      Labels.append("Pdb #%d" % Pdbs.index(p))
    else:
      Labels.append(p)
  #get sequences and standardize
  Seqs = [GetSeq(p) for p in Pdbs]
  TarSeq = sequence.SeqToList(Seq)
  #get the alignments
  AlignList, Weights = GetAlign(TarSeq, Seqs, PdbWeights)
  #make pdbs
  if not LogFile is None: f = file(LogFile, "w")
  if MC:
    pList, LoopResList = [], []
    #shorten the list, and modify weights
    AlignList = AlignList[0:MaxAlign]
    NAlign = len(AlignList)
    Weights = Weights[0:MaxAlign]
    Weights = Weights / Weights.sum()
    #sort through alignments and make proteinclasses
    for i, Align in enumerate(AlignList):
      s = "Making MC alignment using " + ", ".join(
        ["%s[%d:%d]" % (Labels[m.ID], m.c+1, m.d+1)
        for m in Align if not m.ID is None])
      s += "\n  weight is %.3f" % Weights[i]
      NewProt = MakePdbFromAlign(Align, TarSeq, Pdbs)
      pList.append(NewProt)
      #init the loop residue numbers
      LoopRes = []
      #search through all non-terminal extended residues
      for i, m in enumerate(Align):
        #add extended res to the loop if not on the terminal
        if m.ID is None and i > 0 and i < len(Align) - 1:
          a = int((m.a + m.b - 1)*0.5 - 0.5*MinLoopRes)
          b = int((m.a + m.b - 1)*0.5 + 0.5*MinLoopRes) + 1
          a = max(min(a, m.a), 0)
          b = min(max(b, m.b), len(TarSeq))
          LoopRes.extend(range(a, b))
        #now make sure there are enough res between fragments
        m0 = Align[i-1]
        if i > 0 and not m.ID is None and not m0.ID is None:
          a = int((m.a + m0.b)*0.5 - 0.5*MinLoopRes)
          b = int((m.a + m0.b)*0.5 + 0.5*MinLoopRes) + 1
          a = max(a, 0)
          b = min(b, len(TarSeq))
          LoopRes.extend(range(a, b))
      #now sort and get rid of duplicates
      LoopRes.sort()
      LoopRes = [x for (i,x) in enumerate(LoopRes) if not x in LoopRes[:i]]
      LoopResList.append(LoopRes)
      s += "\n  loop residues are " + str([x+1 for x in LoopRes])
      s += "\n"
      print s
      if not LogFile is None: f.write(s + "\n")
    #close the logfile
    if not LogFile is None: f.close()
    #run the assembly
    GenPdbs = assemble.OptimizeMCMulti(pList, LoopResList, Prefix, MaxConf,
      Weights = Weights, NGenEach = MCNGen, RMSDTol = 2.0,
      Verbose = Verbose, Summary = True, MCLen = MCLen)
    #now process the files to add caps or dehydrogen
    for (ScoreVal, PdbFile) in GenPdbs:
      Prot = protein.ProteinClass(Pdb = PdbFile)
      #add caps
      Prot = Prot.Cap(CapN = CapN, CapC = CapC)
      #dehydrogen
      if Dehydrogen: Prot.Dehydrogen()
      Prot.WritePdb(PdbFile)
  else:
    i = 0
    while i < MaxAlign and len(AlignList) > 0:
      Align = AlignList.pop(0)
      s = "Making alignment using " + ", ".join(
        ["%s[%d:%d]" % (Labels[m.ID], m.c+1, m.d+1)
        for m in Align if not m.ID is None])
      print s
      if not LogFile is None: f.write(s + "\n")
      Prot = MakePdbFromAlign(Align, TarSeq, Pdbs)
      #add caps
      Prot = Prot.Cap(CapN = CapN, CapC = CapC)
      if Dehydrogen: Prot.Dehydrogen()
      if CheckOverlap and Prot.HasOverlap(OverlapDist):
        s = "  OVERLAP DETECTED, skipping."
        print s
        if not LogFile is None: f.write(s + "\n")
      else:
        i += 1
        fn = Prefix + "%02d.pdb" % i
        s = "  Saving %s" %fn
        print s
        if not LogFile is None: f.write(s + "\n")
        Prot.WritePdb(fn)
    if not LogFile is None: f.close()


if __name__ == "__main__":
  def FileCmp(x, y):
    return cmp(os.path.basename(x), os.path.basename(y))
  Args = scripttools.ParseArgs(sys.argv, {"cap":None, "dehydrogen":None,
                                          "max":10, "overlapok":None,
                                          "MC":None, "len":1.})
  Seq = sequence.SeqToList(Args[1])
  Prefix = Args[2]
  Cap = Args["cap"]
  Dehydrogen = Args["dehydrogen"]
  CheckOverlap = not Args["overlapok"]
  MC = Args["MC"]
  MaxAlign = Args["max"]
  MCLen = Args["len"]
  Pdbs = [Args[i] for i in range(3, Args["NARG"])]
  Pdbs.sort(cmp = FileCmp, reverse = True)
  MakePdbAlignments(Seq, Pdbs, Prefix, MaxAlign = MaxAlign,
                    Cap = Cap, CheckOverlap = CheckOverlap,
                    MC = MC, MCLen = MCLen, Verbose = False)
