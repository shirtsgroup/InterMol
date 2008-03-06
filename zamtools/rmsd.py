#!/usr/bin/env python

#LAST MODIFIED: 05-01-07

Usage = """Calculates rmsd between pdb files.

Usage     : rmsd.py [OPTIONS] PDBFILE1 PDBFILE2 ...

OPTIONS   : "--first" will only compare the first to all other pdb files
            "--align" to align the sequences before rmsd comparison
            "--compres='1,3-5'" to use residues 1 and 3-5 to minimize RMSD
            "--calcres='1-5,9'" to calculate rmsd only for residues 1-5 and 9
"""

#check for instructions
import sys
if __name__ == "__main__" and len(sys.argv) == 1:
  print Usage
  sys.exit()

from numpy import *  
import copy, os, coords, random
import geometry, sequence, protein, scripttools


def RMSD(Pos1, Pos2, Align = False, Center = True,
         CompInd = None, CalcInd = None, Verbose = False,
         RetAlignment = False):
  """Calculates the RMSD between two conformations.
* Pos1: array of dimensions [N,3] for conformation 1
* Pos2: array of dimensions [N,3] for conformation 2
* Align: True = modify Pos2 to be aligned to Pos1
  (default is False)
* CompInd: indices in [0,N) for positions in Pos to perform alignment
* CalcInd: indices in [0,N) for positions in Pos to compute RMSD
* Verbose: true to report shape/size mismatches
* RetAlignment: True to return the translation vectors and rotation matrix"""
  #clean indices
  AllInd = arange(len(Pos1), dtype = int)
  if all(CompInd == AllInd): CompInd = None
  if all(CalcInd == AllInd): CalcInd = None
  #get indices
  if CompInd is None:
    p1, p2, n = Pos1, Pos2, len(Pos1)
  else:
    p1, p2, n = Pos1.take(CompInd, axis=0), Pos2.take(CompInd, axis=0), len(CompInd)
  #check for correct shapes
  if not shape(p1) == shape(p2):
    if Verbose: print "Position vectors are not the same size."
    return
  elif not len(shape(p1)) == 2:
    if Verbose: print "Position vectors are not the correct rank."
    return
  #get alignment
  Pos1Vec, Pos2Vec, RotMat, Resid = geometry.AlignmentRMSD(p1, p2, Center = Center)
  #compute rmsd
  if not all(CompInd == CalcInd):
    if CalcInd is None:
      p1, p2, n = Pos1, Pos2, len(Pos1)
    else:
      p1, p2, n = Pos1.take(CalcInd, axis=0), Pos2.take(CalcInd, axis=0), len(CalcInd)
    p1 = p1 + Pos1Vec
    p2 = dot(p2 + Pos2Vec, RotMat)
    Resid = sum((p1 - p2)**2, axis=None)
  r = sqrt(Resid / float(n))
  #align Pos2 to Pos1
  if Align: Pos2[:,:] = dot(Pos2 + Pos2Vec, RotMat) - Pos1Vec
  if RetAlignment:
    return r, Pos1Vec, Pos2Vec, RotMat
  else:
    return r

def RMSDProteinClass(p1, p2, Center = True, Backbone = True, AlignSeq = False,
                     CompResInd = None, CalcResInd = None):
  "Returns (RMSD,NRes)"
  #align the sequences
  if AlignSeq:
    Map = sequence.SeqMapClass(p1.Seq, p2.Seq)
    p1 = p1[Map.a:Map.b]
    p2 = p2[Map.c:Map.d]
  #filter res indices
  n = min(len(p1), len(p2))
  if not CompResInd is None: CompResInd = [i for i in CompResInd if i < n]
  if not CalcResInd is None: CalcResInd = [i for i in CalcResInd if i < n]
  #get the right atoms
  BBAtoms = ["N", "CA", "C"]
  if Backbone:
    BBInd1 = p1.AtomInd(AtomName = BBAtoms)
    BBInd2 = p2.AtomInd(AtomName = BBAtoms)
    Pos1, Pos2 = p1.Pos.take(BBInd1, axis=0), p2.Pos.take(BBInd2, axis=0)
    CompInd = p1.AtomInd(AtomName = BBAtoms, ResNum = CompResInd)
    CalcInd = p1.AtomInd(AtomName = BBAtoms, ResNum = CalcResInd)
    CompInd = [i for (i,j) in enumerate(BBInd1) if j in CompInd]
    CalcInd = [i for (i,j) in enumerate(BBInd1) if j in CalcInd]
  else:
    CompInd = p1.AtomInd(ResNum = CompResInd)
    CalcInd = p1.AtomInd(ResNum = CalcResInd)
    Pos1, Pos2 = p1.Pos, p2.Pos
  #compute rmsd
  r = RMSD(Pos1, Pos2, Align = False, Center = Center,
           CompInd = CompInd, CalcInd = CalcInd)
  #determine num of residues used in rmsd calculation
  if CalcResInd is None:
    NRes = len(p1)
  else:
    NRes = len(CalcResInd)
  return r, NRes

def RMSDPdb(PdbFile1, PdbFile2, Center = True, Backbone = True,
            AlignSeq = False, CompResInd = None, CalcResInd = None):
  p1 = protein.ProteinClass(Pdb = PdbFile1)
  p2 = protein.ProteinClass(Pdb = PdbFile2)
  return RMSDProteinClass(p1, p2, Center = Center, Backbone = Backbone,
                          AlignSeq = AlignSeq, CompResInd = CompResInd,
                          CalcResInd = CalcResInd)

def ClusterMSS(CoordsObj, Cutoff, MaxIter = 3, MaxCluster = None,
  MaxClusterWork = None, Method = 0, AvgConf = False, ClusterTol = 0.0,
  Verbose = True):
  """Clusters conformations in a trajectory based on RMSD distance.
* CoordsObj: an object exposing the functions GetNextCoords() which
  returns an array object of the next set of coordinates (or None
  if at the end) and has optional list variable Mask, and
  Reset() which moves back to the first set of coordinates
* Cutoff: maximum RMSD distance of a configuration to a cluster
* MaxIter: maximum number of iterations to perform
* MaxCluster: maximum number of clusters; negative values will force
  all configs to be members of a cluster (default is none)
* MaxClusterWork: maximum numeber of working clusters
* Method: 0 to cluster such that each RMSD between a configuration
  and the average cluster configuration is below Cutoff; 1 to
  cluster such that the average such RMSD is below Cutoff; 2 to
  use both 0 and 1; 3 is similar to 0 except no alignment is
  performed
* AvgConf: True will return an average of configurations in each
  cluster; False returns the configuration with minimum RMSD
  from the average configuration (default is False)
* ClusterTol: float less than 1.0; cluster pairs whose mutual
  RMSD is less than ClusterTol * Cutoff will be combined.  Values
  close to 1.0 may lead to large numbers of iterations."""
  def ClusterRMSD(PosSum, PosSqSum, NAdd):
    "Calculates the average rmsd among configurations within a cluster."
    d = 1. / float(size(PosSum,0))
    b = 1. / float(NAdd)
    rmsdsq = (b*sum(sum(PosSqSum)) - b*b*sum(sum(PosSum**2))) * d
    rmsdsq = max([rmsdsq,0.])
    return sqrt(rmsdsq)
  def BasicRMSD(Pos1, Pos2):
    "Calculates the rmsd between two configurations without alignment."
    rmsdsq = sum(sum((Pos1-Pos2)**2)) / float(size(Pos1,0))
    rmsdsq = max([rmsdsq,0.])
    return sqrt(rmsdsq)
  Iteration = 0   #iteration number
  NAdd = []       #number of configs added to each cluster
  Pos = []        #list of cluster configuration arrays
  PosSq = []      #list of squared cluster config arrays
  FinalIters = 0  #number of iterations without additions/deletions of clusters
  while FinalIters < 2:
    Iteration += 1
    FinalIters += 1
    if Iteration > MaxIter:
      if Verbose: print "Did not converge within maximum number of iterations"
      break
    if Verbose: print "Cluster iteration %d" % Iteration
    if Verbose: print "Starting with %d clusters" % len(Pos)
    CoordsObj.Reset()
    CurInd = 0    #index of the current config
    ClustNum = [] #cluster number of each configuration, starting at 1
    CurPos = CoordsObj.GetNextCoords()
    NAddThis = [0]*len(NAdd)   #number of configs added to each cluster this iteration
    while not CurPos is None:
      ind = -1  #cluster number assigned to this config; -1 means none
      #calculate the rmsd between this configuration and each cluster config,
      #but stop when a rmsd is found which is below the cutoff
      for i in range(0,len(Pos)):
        if Method == 0:
          #rmsd between new config and average cluster config incl new config
          r = RMSD((Pos[i]+CurPos) / float(NAdd[i]+1), CurPos, Align = True, Center = True)
        elif Method == 1:
          #average rmsd of cluster with new config added
          r = ClusterRMSD(Pos[i] + CurPos, PosSq[i] + CurPos**2, NAdd[i]+1)
        elif Method == 2:
          #maximum of the two previous methods
          r1 = ClusterRMSD(Pos[i] + CurPos, PosSq[i] + CurPos**2, NAdd[i]+1)
          r2 = RMSD((Pos[i]+CurPos) / float(NAdd[i]+1), CurPos, Align = True, Center = True)
          r = max([r1,r2])
        else:
          #rmsd between new config and average cluster config incl new config, without alignment
          r = BasicRMSD((Pos[i]+CurPos) / float(NAdd[i]+1), CurPos)
        if r < Cutoff:
          #go with a cluster if rmsd is within the cutoff
          ind = i
          break
      if ind >= 0:
        #add the configuration to the cluster
        Pos[ind] = Pos[ind] + CurPos
        PosSq[ind] = PosSq[ind] + CurPos**2
        NAdd[ind] += 1
        NAddThis[ind] += 1
        ClustNum.append(ind+1)
      elif len(Pos) < MaxClusterWork or MaxClusterWork is None:
        #create a new cluster with this config, as long as it
        #doesn't exceed the maximum number of working clusters
        if Verbose: print "Adding new cluster using config %d (%d clusters total)" % (CurInd+1, len(Pos)+1)
        Pos.append(CurPos)
        PosSq.append(CurPos**2)
        NAdd.append(1)
        NAddThis.append(1)
        ClustNum.append(len(Pos))
        FinalIters = 0
      else:
        #cluster is nothing
        ClustNum.append(0)
        FinalIters = 0
      CurPos = CoordsObj.GetNextCoords()
      CurInd += 1
    i = 0
    #loop through clusters
    while i < len(Pos):
      #remove clusters that have no additions this iteration
      if NAddThis[i] == 0:
        if Verbose: print "Removing cluster %d" % (i+1,)
        del Pos[i]
        del PosSq[i]
        del NAdd[i]
        del NAddThis[i]
        for k in range(0,len(ClustNum)):
          if ClustNum[k] > i + 1:
            ClustNum[k] -= 1
          elif ClustNum[k] == i + 1:
            ClustNum[k] = -1
        FinalIters = 0
      else:
        i += 1
    i = 0
    #find clusters that are too close together and try to combine them;
    #should try to avoid doing this since it prevents convergence sometimes
    #by getting stuck into a combine-add new cluster loop
    while i < len(Pos):
      j = i + 1
      while j < len(Pos):
        if Method == 1:
          r = RMSD(Pos[i]/float(NAdd[i]), Pos[j]/float(NAdd[j]), Align = False, Center = True)
        elif Method == 2:
          r = ClusterRMSD(Pos[i] + Pos[j], PosSq[i] + PosSq[j], NAdd[i] + NAdd[j])
        elif Method == 2:
          r1 = RMSD(Pos[i]/float(NAdd[i]), Pos[j]/float(NAdd[j]), Align = False, Center = True)
          r2 = ClusterRMSD(Pos[i] + Pos[j], PosSq[i] + PosSq[j], NAdd[i] + NAdd[j])
          r = max([r1,r2])
        elif Method == 3:
          r = BasicRMSD(Pos[i]/float(NAdd[i]), Pos[j]/float(NAdd[j]))
        if r < Cutoff * ClusterTol:
          if Verbose: print "Combining cluster %d and %d" % (i+1,j+1)
          Pos[i] = Pos[i] + Pos[j]
          PosSq[i] = PosSq[i] + PosSq[j]
          NAdd[i] = NAdd[i] + NAdd[j]
          del Pos[j]
          del PosSq[j]
          del NAdd[j]
          for k in range(0,len(ClustNum)):
            if ClustNum[k] > j + 1:
              ClustNum[k] -= 1
            elif ClustNum[k] == j + 1:
              ClustNum[k] = i
          FinalIters = 0
        else:
          j += 1
      i += 1
  if Verbose: print "Calculating average structures"
  for i in range(0,len(Pos)):
    Pos[i] = Pos[i] / float(NAdd[i])
  #get cluster populations
  ClustPop, Pos, ClustNum = __CalcClustPop(Pos, ClustNum, MaxCluster, Verbose)
  #if there is a maximum cluster specification that's negative, force
  #everything to the closest cluster
  if not MaxCluster == None and MaxCluster < 0:
    ClustPop, Pos, ClustNum = __CalcForceClust(CoordsObj, ClustPop, Pos, ClustNum, Verbose)
  #calculate final rmsd values for configs and clusters
  ConfRmsd, ClustRmsd = __CalcRmsd(CoordsObj, Pos, ClustNum, AvgConf, Verbose)
  if Verbose: print "%d configurations sorted into %d clusters" % (CurInd,len(Pos))
  return Pos, ClustNum, ClustPop, ConfRmsd, ClustRmsd


def __CalcClustPop(Pos, ClustNum, MaxCluster = None, Verbose = True):
  if Verbose: print "Reordering clusters by population"
  #count the number of configs in each cluster and sort by most populous
  Counts = [(ClustNum.count(i+1)+ClustNum.count(-i-1), i) for i in range(0,len(Pos))]
  Counts.sort()
  Counts.reverse()
  #resort the position array
  Pos = [Pos[j] for (i,j) in Counts]
  #create a dictionary which will tell us the new cluster
  #number for a given old cluster number
  Trans = {0:0}
  for i in range(0,len(Pos)):
    ind = Counts[i][1] + 1
    Trans[ind] = i + 1
    Trans[-ind] = -i - 1
  #update ClustNum with the rearranged cluster numbers
  ClustNum = [Trans[i] for i in ClustNum]
  #update the cluster population
  ClustPop = [i for (i,j) in Counts]
  #crop off any extraneous clusters; clusterless configs
  #are assigned a cluster index of 0
  if len(Pos) > abs(MaxCluster) and not MaxCluster == None:
    del Pos[abs(MaxCluster):]
    del ClustPop[abs(MaxCluster):]
    for i in range(0,len(ClustNum)):
      if abs(ClustNum[i]) > abs(MaxCluster): ClustNum[i] = 0        
  return ClustPop, Pos, ClustNum


def __CalcForceClust(CoordsObj, ClustPop, Pos, ClustNum, Verbose = True):
  #count the number of clusterless configurations
  c = reduce(lambda x, y: x + int(y==0), ClustNum, 0)
  if Verbose: print "Forcing %d extraneous configurations to existing clusters" % c
  #find the nearest cluster to each clusterless config and assign it
  CoordsObj.Reset()
  CurPos = CoordsObj.GetNextCoords()
  for j in range(0,len(ClustNum)):
    if ClustNum[j] == 0:
      ind = -1
      minr = 0.
      for i in range(0,len(Pos)):
        r = RMSD(Pos[i], CurPos, Align = False, Center = True)
        if r < minr or ind < 0: 
          ind = i
          minr = r
      ClustNum[j] = ind + 1
      ClustPop[ind] += 1
    CurPos = CoordsObj.GetNextCoords()
  return ClustPop, Pos, ClustNum
  

def __CalcRmsd(CoordsObj, Pos, ClustNum, AvgConf = False, Verbose = True):
  if Verbose: print "Calculating cluster rmsd values"
  #calculate the pairwise cluster rmsd values
  ClustRmsd = zeros((len(Pos),len(Pos)),float)
  for i in range(0,len(Pos)):
    for j in range(i+1,len(Pos)):
      ClustRmsd[i,j] = RMSD(Pos[i], Pos[j], Align=False, Center=True)
      ClustRmsd[j,i] = ClustRmsd[i,j]
  if Verbose: print "Calculating final rmsd values"
  #loop through configs and find the one with the lowest
  #rmsd in each cluster
  CoordsObj.Reset()
  CurInd = 0
  CurPos = CoordsObj.GetNextCoords()
  ConfRmsd = -1. * ones((len(ClustNum),), float)
  if not AvgConf:
    MinRmsd = [-1]*len(Pos)
  while not CurPos is None:
    i = abs(ClustNum[CurInd]) - 1
    if i >= 0:
      ConfRmsd[CurInd] = RMSD(Pos[i], CurPos, Align=False, Center=True)
      if not AvgConf:
        if MinRmsd[i] < 0:
          MinRmsd[i] = CurInd
        elif ConfRmsd[MinRmsd[i]] > ConfRmsd[CurInd]:
          MinRmsd[i] = CurInd
    CurPos = CoordsObj.GetNextCoords()
    CurInd += 1
  #loop through the configs again and extract the
  #coords of the minimum-rmsd configs for each clust
  if not AvgConf:
    if Verbose: print "Finding nearest cluster structures"
    for ind in MinRmsd:
      ClustNum[ind] *= -1
    CoordsObj.Reset()
    CurInd = 0
    CurPos = CoordsObj.GetNextCoords(Mask = coords.NoMask)
    while not CurPos is None:
      if ClustNum[CurInd] < 0:
        Pos[MinRmsd.index(CurInd)] = CurPos
      CurInd += 1
      CurPos = CoordsObj.GetNextCoords(Mask = coords.NoMask)
  return ConfRmsd, ClustRmsd


def SaveClustResults(Pos, ClustNum, ClustPop, ConfRmsd, ClustRmsd,
  Prefix = "clust", ConfIndices = None, Verbose = False):
  #make the indices
  if ConfIndices is None:
    ConfIndices = range(0, len(ConfRmsd))
  #calculate the percent
  ClustPct = [100.*float(ClustPop[i])/float(len(ClustNum))
              for i in range(0,len(Pos))]
  #save results file
  s = "CLUSTER POPULATION:\nCluster number, Number of configs, Percent\n"
  s += "\n".join(["%-5d  %-7d  %.2f" % (i+1, ClustPop[i], ClustPct[i])
                  for i in range(0,len(Pos))])
  s += "\n\nCLUSTER-TO-CLUSTER RMSD:\nCluster number, Cluster number, RMSD\n"
  for i in range(0,len(Pos)):
    for j in range(i+1,len(Pos)):
      s += "%-5d  %-5d  %-8.2f\n" % (i+1, j+1, ClustRmsd[i,j])
  s += "\n\nCLUSTER MEMBERSHIP:\nConfig number, Cluster number, RMSD\n"
  s += "\n".join(["%-7d  %-5d  %-8.2f" % (ConfIndices[i], ClustNum[i],
    ConfRmsd[i]) for i in range(0,len(ConfRmsd))]) + "\n"
  fn = Prefix + "results.txt"
  file(fn,"w").write(s)
  #save cluster files
  if Verbose:
    for i in range(0,len(Pos)):
      data = [(ConfRmsd[j], ConfIndices[j]) for j in range(0,len(ConfRmsd))
              if abs(ClustNum[j]) == i + 1]
      data.sort()
      s = "\n".join(["%-7d  %-8.2f" % (j[1],j[0]) for j in data]) + "\n"
      fn = Prefix + "%04d-%dpc.txt" % (i+1, int(ClustPct[i]))
      file(fn, "w").write(s)



#======== COMMAND-LINE RUNNING ========

def GetResList(Arg):
  if Arg is None or Arg == "": return None
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

def RepRMSD(r):
  if r is None:
    return 'NA'
  else:
    return "%-8.2f" % r

if __name__ == "__main__":
  Args = scripttools.ParseArgs(sys.argv[1:])
  #get pdbs
  Pdbs = [f for f in Args["ARGS"] if os.path.isfile(f)]
  #get options
  Align = "align" in Args["FLAGS"]
  CompResInd = GetResList(Args.get("compres", None))
  CalcResInd = GetResList(Args.get("calcres", None))
  #check for files
  for f in [f for f in Args["ARGS"] if not os.path.isfile(f)]:
    print "Could not find %s." % f
  N = len(Pdbs)
  Filenames = [os.path.basename(x) for x in Pdbs]
  if N <= 1:
    print "Nothing to compare."
    sys.exit()
  if "first" in Args["FLAGS"]:
    MaxLen1 = len(Filenames[0])
    MaxLen2 = max([len(fn) for fn in Filenames[1:]])
    print "%-*s %-*s %-8s %-8s %-5s" % (MaxLen1, "Pdb1", MaxLen2, "Pdb2",
                                        "BB_RMSD", "All_RMSD", "NRes")
    p1 = protein.ProteinClass(Pdb = Pdbs[0])
    for j in range(1, N):
      p2 = protein.ProteinClass(Pdb = Pdbs[j])
      x1, NRes = RMSDProteinClass(p1, p2, Backbone = True, AlignSeq = Align,
                                  CompResInd = CompResInd, CalcResInd = CalcResInd)
      x2, NRes = RMSDProteinClass(p1, p2, Backbone = False, AlignSeq = Align,
                                  CompResInd = CompResInd, CalcResInd = CalcResInd)
      print "%-*s %-*s %-8s %-8s %-5d" % (MaxLen1, Filenames[0],
                                          MaxLen2, Filenames[j],
                                          RepRMSD(x1), RepRMSD(x2), NRes)
  else:
    MaxLen = max([len(fn) for fn in Filenames])
    print "%-*s %-*s %-8s %-8s %-5s" % (MaxLen, "Pdb1", MaxLen, "Pdb2",
                                        "BB_RMSD", "All_RMSD", "NRes")
    Prots = [protein.ProteinClass(Pdb = f) for f in Pdbs]
    for i in range(0, N):
      for j in range(i+1, N):
        x1, NRes = RMSDProteinClass(Prots[i], Prots[j], Backbone = True, AlignSeq = Align,
                                    CompResInd = CompResInd, CalcResInd = CalcResInd)
        x2, NRes = RMSDProteinClass(Prots[i], Prots[j], Backbone = False, AlignSeq = Align,
                                    CompResInd = CompResInd, CalcResInd = CalcResInd)
        print "%-*s %-*s %-8s %-8s %-5d" % (MaxLen, Filenames[i],
                                            MaxLen, Filenames[j],
                                            RepRMSD(x1), RepRMSD(x2), NRes)
      
      