#!/usr/bin/env python

#LAST MODIFIED: 04-18-07

Usage = """Runs clustering on pdb files in a directory or on an amber trajectory.
Produces clustresults.txt, clust####.txt, and clust####.pdb files."

Usage     : cluster.py PDBFILES [OPTIONS]
   OR     : cluster.py TRJFILE PRMTOPFILE [OPTIONS]

PDBFILES  : list of pdb files
TRJFILE   : trajectory CRD file (can be gzipped)
PRMTOPFILE: PARM7 file
OPTIONS   : "--rmsd=X" to set the rmsd tolerance (default 2.0)
            "--nskip=X" number of configs in trajectory to skip (default is 0)
            "--nread=X" number of configs in trajectory to read; -1 is all (default -1)
            "--nstride=X" read configs every nstride frames (default is 1)
            "--allatom" to use all atom rmsd (default is backbone)
            "--alphacarbon" to use just alpha carbon (default is backbone)
            "--maxclust=X" maximum number of clusters; negative values will force all
                           configurations to a cluster (default 10)
            "--maxiter=X" maximum number of iterations (default 10)
            "--prefix=X" to change the output prefix (default "clust")
"""

#check for instructions
import sys
if len(sys.argv) == 1:
  print Usage
  sys.exit()


import sys, glob, os
import rmsd, coords, scripttools

#parse arguments
Args = scripttools.ParseArgs(sys.argv)
NSkip = int(Args.get("nskip", 0))
NRead = int(Args.get("nread", -1))
NStride = int(Args.get("nstride", 1))
MaxCluster = int(Args.get("maxclust", 10))
MaxIter = int(Args.get("maxiter", 10))
RmsdTol = float(Args.get("rmsd", 2.0))
Prefix = Args.get("prefix", "clust")

#decide mask
if "allatom" in Args["FLAGS"]:
  Mask = coords.NoMask
  print "Using all atoms"
elif "alphacarbon" in Args["FLAGS"]:
  Mask = coords.AlphaCarbonMask
  print "Using alpha carbons"
else:
  Mask = coords.BackboneMask
  print "Using backbone atoms"

#decide mode
if ".pdb" in Args[1]:
  Mode = 0
  PdbFiles = Args["ARGS"][1:]
  f = coords.PdbListClass(PdbFiles, Mask = Mask)
else:
  Mode = 1
  TrjFile, PrmtopFile = Args[1], Args[2]
  f = coords.TrjClass(TrjFile, PrmtopFile, Mask = Mask,
                      NRead = NRead, NSkip = NSkip, NStride = NStride)

#run the cluster algorithm
Pos, ClustNum, ClustPop, ConfRmsd, ClustRmsd = rmsd.ClusterMSS(f,
  RmsdTol, MaxIter = MaxIter, MaxCluster = MaxCluster, Method = 0)
Indices = f.GetPastIndices()

if Mode == 0:
  #save output pdbs  
  for i in range(len(Pos)):
    pc = 100. * float(ClustPop[i]) / float(sum(ClustPop))
    fn = Prefix + "%02d_%.0fpc.pdb" % (i+1, pc)
    if -i in ClustNum:
      j = ClustNum.index(-i)
    else:
      j = 0
    coords.SavePdbCoords(Pos[i], fn, PdbFiles[j])
  #save results
  rmsd.SaveClustResults(Pos, ClustNum, ClustPop, ConfRmsd, ClustRmsd,
    ConfIndices = Indices, Verbose = False, Prefix = Prefix)
  
elif Mode == 1:
  #save output configs
  for i in range(len(Pos)):
    pc = 100. * float(ClustPop[i]) / float(len(ClustNum))
    fn = Prefix + "%02d_%.0fpc.pdb" % (i+1, pc) 
    coords.SavePdbCoordsAmb(Pos[i], PrmtopFile, fn)
  #save results
  rmsd.SaveClustResults(Pos, ClustNum, ClustPop, ConfRmsd, ClustRmsd,
    ConfIndices = Indices, Verbose = False, Prefix = Prefix)

