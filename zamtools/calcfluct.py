#!/usr/bin/env python

#LAST MODIFIED: 06-20-06

Usage = """Calculates residue-specific fluctuations between a trajectory
and a reference pdb file.  Produces a new pdb file with "-fluct" appended.

Usage     : calcfluct.py TRJFILE PRMTOPFILE REFPDB [NSKIP NREAD NSTRIDE]

trjfile   : trajectory CRD file (can be gzipped)
prmtopfile: PARM7 file
nskip     : number of configs in trajectory to skip (default is 0)
nread     : number of configs in trajectory to read; -1 is all (default -1)
nstride   : read configs every nstride frames (default is 1)
"""

#check for instructions
import sys
if len(sys.argv) == 1:
  print Usage
  sys.exit()

import coords, rmsd, protein, os
from numpy import *

TrjFile = sys.argv[1]
PrmtopFile = sys.argv[2]
PdbFile = sys.argv[3]
if len(sys.argv) > 4:
  NSkip = int(sys.argv[4])
  NRead = int(sys.argv[5])
  NStride = int(sys.argv[6])
else:
  NSkip = 0
  NRead = -1
  NStride = 1

RefPos = coords.GetPdbCoords(PdbFile)
t = coords.TrjClass(TrjFile, PrmtopFile, 
    NRead = NRead, NSkip = NSkip, NStride = NStride)

NAtom = len(RefPos)
print "Found %d atoms" % NAtom
dSq = zeros(NAtom, float)
N = 0

for Pos in t:
  N += 1
  r = rmsd.RMSD(RefPos, Pos, Align = True)
  dSq += ((RefPos - Pos)**2).sum(axis=1)
  if N % 100 == 0: print "Analyzed %d frames" % N
dSq = dSq/N
dSq = sqrt(dSq)
print "Analyzed %d frames" % N

#s = "FLUCTUATION RESULTS\n"
#s += "atom number, root-mean-square fluctuation\n"
#for i in range(0,NAtom):
#  s += "%-4d %-8.2f\n" % (i+1, dSq[i])
#file("fluctresults.txt", "w").write(s)

p = protein.ProteinClass(Pdb = PdbFile)
i = 0
for r in p.Res:
  for a in r.Atoms:
    a.BFactor = dSq[i]
    i += 1
p.WritePdb(os.path.basename(PdbFile).replace(".pdb","-fluct.pdb"))
