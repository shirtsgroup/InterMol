#!/usr/bin/env python

#LAST MODIFIED: 05-01-07

#DESCRIPTION: Provides itemized scoring functions for ProteinClass in protein.py


Usage = """Runs scoring metrics on pdb files.

Usage     : scoreconf.py [OPTIONS] PDBFILE1 PDBFILE2 ...
   OR     : scoreconf.py [OPTIONS] TRJFILE PRMTOPFILE [NSKIP NREAD NSTRIDE]

OPTIONS   : "--file=FILENAME" to save results to FILENAME
            "--normalize" to scale pairs and CO by 1/N, Rg by 1/N^(1/3)
            "--omitpath" to use only the base name of files
            "--width=N" to format for N characters wide
            "--verbose" for verbose output
            "--refpdb=PDB" trim pdb files to be aligned to file PDB
            "--stats" to report averages and std devs of scores
PDBFILE*  : pdb files to score
TRJFILE   : trajectory CRD file to score (can be gzipped)
PRMTOPFILE: PARM7 file
NSKIP     : number of frames to skip (default is 0)
NREAD     : number of frames to read; -1 is all (default all)
NSTRIDE   : read every NSTRIDE frames (default is 1)
"""
import sys
if __name__ == "__main__" and len(sys.argv) == 1:
  print Usage
  sys.exit()


from numpy import *
import sys, os
import protein, sequence, geometry, coords, rmsd
import scripttools


#LIST OF ALL FUNCTIONS AND LABELS
NRefFunc = 2
FuncList = [
  #label         format    function    normalization
  ("BBRMSD",     "%-9.3f", lambda p, pRef: rmsd.RMSDProteinClass(p, pRef, Center = True, Backbone = True, AlignSeq = True)[0],
                           lambda p, pRef: 1. ),
  ("NResAlgn",   "%-d",    lambda p, pRef: len(sequence.SeqMapClass(pRef.Seq, p.Seq)),
                           lambda p, pRef: 1. ),
  ("NRes",       "%-d",    lambda p: len(p.Seq),
                           lambda p: 1 ),
  ("NHHPr",      "%-9.3f", lambda p: p.NPhobicContact(),
                           lambda p: len(p) ),
  ("NSaltPr",    "%-9.3f", lambda p: p.NSaltContact(),
                           lambda p: len(p) ),
  ("RgA",        "%-9.3f", lambda p: p.RadiusOfGyration(),
                           lambda p: float(len(p))**(1./3.) ),
  ("RgH",        "%-9.3f", lambda p: p.RadiusOfGyration(ResInd = p.ResInd(Hydrophobic=True)),
                           lambda p: float(len(p))**(1./3.) ),
  ("RgC",        "%-9.3f", lambda p: p.RadiusOfGyration(ResInd = p.ResInd(ResCharged=True)),
                           lambda p: float(len(p))**(1./3.) ),
  ("RgH/RgA",    "%-9.3f", lambda p: p.RadiusOfGyration(ResInd = p.ResInd(Hydrophobic=True)) /
                                     p.RadiusOfGyration(),
                           lambda p: 1. ),
  ("FracAInt",   "%-9.3f", lambda p: p.ResFracInterior(),
                           lambda p: 1. ),
  ("FracHInt",   "%-9.3f", lambda p: p.ResFracInterior(ResInd = p.ResInd(Hydrophobic=True)),
                           lambda p: 1. ),
  ("FracCInt",   "%-9.3f", lambda p: p.ResFracInterior(ResInd = p.ResInd(ResCharged=True)),
                           lambda p: 1. ),
  ("CoordA",     "%-9.3f", lambda p: p.MeanCoordination(),
                           lambda p: 1. ),
  ("CoordH",     "%-9.3f", lambda p: p.MeanCoordination(FromResInd = p.ResInd(Hydrophobic=True)),
                           lambda p: 1. ),
  ("CrH/CrA",    "%-9.3f", lambda p: p.MeanCoordination(FromResInd = p.ResInd(Hydrophobic=True)) /
                                     p.MeanCoordination(),
                           lambda p: 1. ),
  ("CO-A",       "%-9.3f", lambda p: p.MeanContactOrder(),
                           lambda p: len(p) ),
  ("CO-H",       "%-9.3f", lambda p: p.MeanContactOrder(ResInd = p.ResInd(Hydrophobic=True)),
                           lambda p: len(p) ),
  ("NBBHBond",   "%-9.3f", lambda p: sum(p.HBondContactMap()),
                           lambda p: len(p) ),
  ("HBond",      "%-9.3f", lambda p: p.HBondChargeScore(),
                           lambda p: len(p) ),
  ("PhobCont",   "%-9.3f", lambda p: p.ResContactScore(p.ResInd(Hydrophobic=True, ReturnMask=True).astype(int),
                                                       array([[0,0],[0,-1]], float)),
                           lambda p: len(p) )
  ]

ColWidth = 9
ScoreLabels = [Label for (Label, Fmt, Func, Norm) in FuncList]
NScores = len(ScoreLabels)


def Scores(p, Normalize = False, pRef = None):
  "Returns a list of scores for all metrics."
  if Normalize:
    if pRef is None:
      ret = [0. for (Label, Fmt, Func, Norm) in FuncList[:NRefFunc]]
    else:
      ret = [Func(p,pRef)/float(Norm(p,pRef)) for (Label, Fmt, Func, Norm) in FuncList[:NRefFunc]]
    ret = ret + [Func(p)/float(Norm(p)) for (Label, Fmt, Func, Norm) in FuncList[NRefFunc:]]
  else:
    if pRef is None:
      ret = [0. for (Label, Fmt, Func, Norm) in FuncList[:NRefFunc]]
    else:
      ret = [Func(p,pRef) for (Label, Fmt, Func, Norm) in FuncList[:NRefFunc]]
    ret = ret + [Func(p) for (Label, Fmt, Func, Norm) in FuncList[NRefFunc:]]
  return ret

def ScoreDict(Vals):
  "Converts a list of scores to a dictionary."
  d = {}
  for i in range(0, len(Vals)):
    (Label, Fmt, Func, Norm) = FuncList[i]
    d[Label] = Vals[i]
  return d

def ScoreStrList(Vals, Tag = None, TagWidth = ColWidth):
  "Returns a list of strings for a list of score values."
  if Tag is None:
    ValStr = []
  else:
    ValStr = [Tag.ljust(TagWidth)[:TagWidth]]
  for i in range(0, len(FuncList)):
    (Label, Fmt, Func, Norm) = FuncList[i]
    ValStr.append((Fmt % Vals[i]).ljust(ColWidth))
  return ValStr

def ScoreStr(Vals, Tag = None, TagWidth = ColWidth):
  "Returns a string for a list of score values."
  return " ".join(ScoreStrList(Vals, Tag, TagWidth))

def HeadStrList(Tag = None, TagWidth = ColWidth):
  "Returns list of header strings."
  Headers = [Label.ljust(ColWidth) for Label in ScoreLabels]
  if not Tag is None:
    Headers = [Tag.ljust(TagWidth)[:TagWidth]] + Headers
  return Headers

def HeadStr(Tag = None, TagWidth = ColWidth):
  "Returns a header string."
  return " ".join(HeadStrList(Tag, TagWidth))

def FormatTable(Table, Width = 0, KeepFirstCol = True):
  """Formats a table to a specified width."""
  if Width == 0:
    return "\n".join([" ".join(Row) for Row in Table])
  else:
    #find where each of the columns goes
    NSec, Pos = 0, 0
    NCol = len(Table[0])
    if KeepFirstCol:
      SecInd = [1]
    else:
      SecInd = [0]
    for i in range(0, NCol):
      HeadLen = len(Table[0][i])
      if Pos + HeadLen > Width:
        if KeepFirstCol:
          Pos = len(Table[0][0]) + 1
        else:
          Pos = 0
        SecInd.append(i)
      Pos += HeadLen + 1
    if not SecInd[-1] == NCol: SecInd.append(NCol)
    NSec = len(SecInd) - 1
    s = ""
    for i in range(NSec):
      a, b = SecInd[i], SecInd[i+1]
      if KeepFirstCol:
        s += "\n".join([" ".join(Row[0:1] + Row[a:b]) for Row in Table])
      else:
        s += "\n".join([" ".join(Row[a:b]) for Row in Table])
      s += "\n\n"
    return s
    
    

def ScorePdbs(PdbFiles, Normalize = False, BaseName = False,
              Verbose = False, Width = 0, RefPdb = None,
              ReportStats = False):
  """Returns a list of values for each pdb, a list of average
  values, and a string summary of PdbFile scores."""
  if RefPdb is None:
    pRef = None
  else:
    pRef = protein.ProteinClass(Pdb = RefPdb)
  #start a data table
  if BaseName:
    PdbNames = [os.path.basename(f) for f in PdbFiles]
  else:
    PdbNames = PdbFiles
  PdbWidth = max([len(Pdb) for Pdb in PdbNames]) + 1
  Table = [HeadStrList("Pdb", PdbWidth)]
  PdbScores = []
  AvgScores = zeros(NScores, float)
  StdScores = zeros(NScores, float)
  N = 0
  #sort through the pdb files
  for i in range(0, len(PdbFiles)):
    Pdb = PdbFiles[i]
    if not os.path.isfile(Pdb):
      Table.append(["%s not found." % Pdb])
      PdbScores.append([])
    else:
      p = protein.ProteinClass(Pdb = Pdb)
      if not pRef is None:
        Map = sequence.SeqMapClass(Seq1 = pRef.Seq, Seq2 = p.Seq)
        p = p[Map.c:Map.d]
      ThisScores = Scores(p, Normalize, pRef)
      Table.append(ScoreStrList(ThisScores, PdbNames[i], PdbWidth))
      PdbScores.append(ThisScores)
      AvgScores += ThisScores
      StdScores += array(ThisScores, float)**2
      N += 1
    if Verbose: print "Analyzed pdb %s" % Pdb
  StdScores = list(sqrt(StdScores / N - (AvgScores / N)**2))
  AvgScores = list(AvgScores / N)
  if ReportStats:
    Table.append(ScoreStrList(AvgScores, "AVG", PdbWidth))
    Table.append(ScoreStrList(StdScores, "STDEV", PdbWidth))
  Summary = FormatTable(Table, Width)
  return PdbScores, AvgScores, Summary

def ScoreTrj(Trj, Normalize = False, Verbose = False, Width = 0,
             RefPdb = None, ReportStats = False):
  """Returns a list of values for each frame, a list of average
  values, and a string summary for trajectory frames."""
  if RefPdb is None:
    pRef = None
  else:
    pRef = protein.ProteinClass(Pdb = RefPdb)
  #start a data table
  Table = [HeadStrList("Index", 8)]
  TrjScores = []
  AvgScores = zeros(NScores, float)
  StdScores = zeros(NScores, float)
  N = 0.000001
  #link proteinclass
  p = protein.ProteinClass()
  p.LinkTrj(Trj)
  #sort through the trajectory
  Trj.Reset()
  for Pos in Trj:
    ThisScores = Scores(p, Normalize, pRef)
    Ind = "%-8d" % Trj.Index
    Table.append(ScoreStrList(ThisScores, Ind, 8))
    TrjScores.append(ThisScores)
    AvgScores += ThisScores
    StdScores += array(ThisScores, float)**2
    N += 1
    if Verbose and N % 10 == 0: print "Scanned frame %d" % Trj.Index
  #unlink trj class
  p.UnlinkTrj()
  StdScores = list(sqrt(StdScores / N - (AvgScores / N)**2))
  AvgScores = list(AvgScores / N)
  if ReportStats:
    Table.append(ScoreStrList(AvgScores, "AVG", 8))
    Table.append(ScoreStrList(StdScores, "STDEV", 8))
  Summary = FormatTable(Table, Width)
  return TrjScores, AvgScores, Summary  

  
if __name__ == "__main__":
  Args = scripttools.ParseArgs(sys.argv, {"width":0})
  Norm = "normalize" in Args["FLAGS"]
  BaseName = "omitpath" in Args["FLAGS"]
  ReportStats = "stats" in Args["FLAGS"]
  Width = Args["width"]
  Verbose = Args.get("verbose", None)
  RefPdb = Args.get("refpdb", None)
  if "pdb" in Args[1]:
    #get pdb name column with
    PdbFiles = Args["ARGS"][1:]
    #run the scoring routines
    PdbScores, AvgScores, Summary = ScorePdbs(PdbFiles, Norm, BaseName,
                                              Verbose = Verbose, Width = Width,
                                              RefPdb = RefPdb, ReportStats = ReportStats)
  else:
    #get trajectory files and make class
    TrjFile = Args[1]
    PrmtopFile = Args[2]
    if Args["NARG"] > 3:
      NSkip, NRead, NStride = int(Args[3]), int(Args[4]), int(Args[5])
    else:
      NSkip, NRead, NStride = 0, -1, 1
    Trj = coords.TrjClass(TrjFile, PrmtopFile, NSkip=NSkip, NRead=NRead,
                          NStride=NStride)
    TrjScores, AvgScores, Summary = ScoreTrj(Trj, Norm, Verbose = Verbose,
                                             Width = Width, RefPdb = RefPdb,
                                             ReportStats = ReportStats)
    Trj.Close()
  #write to file
  if "file" in Args["FLAGS"]:
    OutFile = Args["file"]
    file(OutFile, "w").write(Summary)
  #print results
  print Summary

   
    







  
