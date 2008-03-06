#!/usr/bin/env python

#LAST MODIFIED: 04-18-07

Usage = """Runs mesostring analysis on an amber trajectory.
Produces mesoresults.txt."

Usage     : mesoanalysis.py trjfile prmtopfile [nskip nread nstride]

trjfile   : trajectory CRD file (can be gzipped)
prmtopfile: PARM7 file
nskip     : number of configs in trajectory to skip (default is 0)
nread     : number of configs in trajectory to read; -1 is all (default -1)
nstride   : read configs every nstride frames (default is 1)
"""

#check for instructions
import sys
if __name__ == "__main__" and len(sys.argv) == 1:
  print Usage
  sys.exit()
  

import sys, os, gzip
import coords, sequence, protein
import geometry
from numpy import *
import copy  #this must be after numpy

#GLOBALS
#phi, psi ranges
PhiRange = [(0.,120.), (-180.,0.), (120.,180.)]
PsiRange = [(-180.,180.), (-135.,45.)]
#for glycine
PhiRangeG = [(-180,0.), (0.,180.)]
PsiRangeG = [(-135.,45.), (-45.,135.)]


def InRange(Phi, Psi, PhiRange, PsiRange):
  return (Phi >= PhiRange[0] and Phi <= PhiRange[1]
          and Psi >= PsiRange[0] and Psi <= PsiRange[1])

def MesoState(Phi, Psi, Gly):
  "Returns the mesostate for a phi, psi pair."
  if Gly:
    if InRange(Phi, Psi, PhiRangeG[0], PsiRangeG[0]):
      return "a"
    elif InRange(Phi, Psi, PhiRangeG[1], PsiRangeG[1]):
      return "l"
    else:
      return "b"
  else:
    if InRange(Phi, Psi, PhiRange[0], PsiRange[0]):
      return "l"
    elif InRange(Phi, Psi, PhiRange[1], PsiRange[1]):
      return "a"
    elif InRange(Phi, Psi, PhiRange[2], PsiRange[1]):
      return "a"
    else:
      return "b"

def MesoString(Dihedrals, Seq):
  "Returns a string of mesostates for a list of dihedrals."
  ms = ""
  Seq = sequence.SeqToList(Seq)
  for i in range(0, len(Dihedrals)):
    if sequence.Cap(Seq[i]): continue
    (Phi, Psi) = Dihedrals[i]
    ms += MesoState(Phi, Psi, Seq[i] == "GLY")
  return ms

def MesoMask(p):
  "Makes a mask for mesostring analysis."
  #Mask=0 for skip, =1 for normal, =2 for GLY
  Mask = []
  for (i,r) in enumerate(p.Res):
    rn = r.Name.upper().strip()
    if rn == "GLY":
      Mask.append(2)
    elif sequence.Cap(rn) or rn in protein.NoDihedrals:
      Mask.append(0)
    else:
      Mask.append(1)
  return Mask

def MesoStringMask(Dihedrals, Mask):
  "Returns a string of mesostates give a mask."
  #Mask=0 for skip, =1 for normal, =2 for GLY
  ms = ""
  for i in range(0, len(Dihedrals)):
    if Mask[i] == 0:
      continue
    else:
      (Phi, Psi) = Dihedrals[i]
      ms += MesoState(Phi, Psi, Mask[i] == 2)
  return ms

def PdbMesoString(PdbFile):
  "Returns a string of mesostates for a pdb file."
  p = protein.ProteinClass(Pdb = PdbFile)
  Dih = [p.PhiPsi(i) for i in range(0, len(p))]
  return Mesostring(Dih, p.Seq)

def RunAnal(CoordsObj, OutputPath = None, Prefix = None):
  "Runs mesostring analysis of a coords object."
  #link proteinclass
  p = protein.ProteinClass()
  p.LinkCoordsObj(CoordsObj)
  #loop through configurations
  ConfMeso = []
  CoordsObj.Reset()
  Mask = MesoMask(p)
  while True:
    Pos = CoordsObj.GetNextCoords()
    if Pos is None: break
    #parse dihedrals
    Dih = []
    for i in range(0, len(p)):
      Dih.append(p.PhiPsi(i))
    ConfMeso.append(MesoStringMask(Dih, Mask))
  #unlink
  p.UnlinkCoordsObj()
  #compute the number of instances of each mesostring
  MesoDict = {}
  for s in ConfMeso:
    if s in MesoDict:
      MesoDict[s] += 1
    else:
      MesoDict[s] = 1
  #put into a list and sort by population
  MesoPop = [(p,s) for (s,p) in MesoDict.iteritems()]
  MesoPop.sort()
  MesoPop.reverse()
  MesoPop = [(s,p) for (p,s) in MesoPop]
  #calculate the mesostate entropy
  Tot = float(len(ConfMeso))
  MesoEntropy = log(Tot) - sum([p * log(p) for (s,p) in MesoPop]) / Tot
  #write the output
  ConfIndices = CoordsObj.GetPastIndices()
  if not OutputPath is None:
    WriteAnalysis(OutputPath, ConfMeso, MesoPop, MesoEntropy, Prefix, ConfIndices)
  return ConfMeso, MesoPop, MesoEntropy

def WriteAnalysis(OutputPath, ConfMeso, MesoPop, MesoEntropy,
  Prefix = None, ConfIndices = None):
  "Writes an analysis file from mesostring results."
  if OutputPath is None: return
  if ConfIndices is None:
    ConfIndices = range(0, len(ConfMeso))
  s = "MESOSTRING ENTROPY:\n%.4f\n" % MesoEntropy
  s += "\nMESOSTRING POPULATION:\nMesostring, Number of configs, Percent\n"
  Tot = float(len(ConfMeso))
  for (ms,p) in MesoPop:
    s += "%s    %-7d  %.2f\n" % (ms, p, float(p)*100./Tot)
  s += "\nCONFIGURATION MESOSTRINGS:\nConfig number, Mesostring\n"
  for i in range(0,len(ConfMeso)):
    s += "%-7d  %s\n" % (ConfIndices[i], ConfMeso[i])
  #write file
  fn = "mesoresults.txt"
  if not Prefix is None: fn = Prefix.strip() + "." + fn
  f = file(os.path.join(OutputPath, fn), "w")
  f.write(s)
  f.close()
    
def RunAnalTrj(TrjFile, PrmtopFile, OutputPath = None, Prefix = None,
  NSkip = 0, NRead = None, NStride = 1):
  "Runs a mesostring analysis of a trajectory."
  #make the coords object
  CoordsObj = coords.TrjClass(TrjFile, PrmtopFile,
      NSkip = NSkip, NRead = NRead, NStride = NStride)
  ConfMeso, MesoPop, MesoEntropy = RunAnal(CoordsObj, OutputPath, Prefix)
  return ConfMeso, MesoPop, MesoEntropy

def RunAnalPdbs(PdbFileList, OutputPath = None, Prefix = None):
  "Runs a mesostring analysis of a list of pdb files."
  #make the coords object
  CoordsObj = coords.PdbListClass(PdbFileList)
  ConfMeso, MesoPop, MesoEntropy = RunAnal(CoordsObj, OutputPath, Prefix)
  return ConfMeso, MesoPop, MesoEntropy


#command-line running
if __name__ == "__main__":
  TrjFile = sys.argv[1]
  PrmtopFile = sys.argv[2]
  if len(sys.argv) > 3:
    NSkip = int(sys.argv[3])
    NRead = int(sys.argv[4])
    NStride = int(sys.argv[5])
  else:
    NSkip = 0
    NRead = None
    NStride = 1
  print "Running mesostring analysis..."
  RunAnalTrj(TrjFile, PrmtopFile, OutputPath = ".",
    NSkip = NSkip, NRead = NRead, NStride = NStride)
