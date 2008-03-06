#!/usr/bin/env python

#CONTAINS ROUTINES FOR ANALYZING TRAJECTORIES;
#ANALAGOUS TO PTRAJ

#LAST MODIFIED: 05-11-07


from numpy import *
import os, sys, glob, gzip
import coords, sequence

#GLOBALS
BlockLen = 13
DistFmt = "%.3f"
VerboseDflt = False
DistMethodDesc = ["alpha carbons", "beta carbons", "residue centroid"]



def MakeAtomGroups(PairAtoms):
  """Finds all atom groups from PairAtoms.
Returns GroupAtoms, PairGroups.
* PairAtoms: list of pairs of lists of atom numbers
* GroupAtoms: list of atoms in each group
* PairGroups: list of group numbers for each pair"""
  if len(PairAtoms) == 0:
    return [], []
  GroupAtoms = []
  PairGroups = []
  for Atoms1, Atoms2 in PairAtoms:
    #put single entries in a list
    if not type(Atoms1) is list:
      Atoms1 = [Atoms1]
    if not type(Atoms2) is list:
      Atoms2 = [Atoms2]  
    #check to see if groups are already there
    if not Atoms1 in GroupAtoms:
      GroupAtoms.append(Atoms1)
    if not Atoms2 in GroupAtoms:
      GroupAtoms.append(Atoms2)
    #update the group pair list
    Group1 = GroupAtoms.index(Atoms1)
    Group2 = GroupAtoms.index(Atoms2)
    PairGroups.append((Group1,Group2))
  return GroupAtoms, PairGroups  

    
def GetResPairAtoms(PrmtopFile, PairList, DistMethod):
  """Makes a list of all the atoms involved in residue pairs.
Returns PairAtoms.
* PrmtopFile: string, path of prmtop file
* PairList: list of pairs of residue numbers
* DistMethod: int, 0 = CAs, 1 = CBs, 2 = residue COM
* PairAtoms: list of pairs of lists of atom numbers"""
  AtomNames = coords.GetPrmtopAtomNames(PrmtopFile)
  AtomRes = coords.GetPrmtopAtomRes(PrmtopFile)
  NAtom = len(AtomNames)
  #make a list of (residue num, atom name)
  AtomDat = [(AtomRes[i], AtomNames[i].strip()) for i in range(0,NAtom)]
  PairAtoms = []
  def GetInd(Dat, ResInd, NameList):
    #first search through name list
    for Name in NameList:
      if (ResInd, Name) in Dat: return Dat.index((ResInd, Name))
    #just get any atom with the residue
    for (i, (ResInd2, Name)) in enumerate(Dat):
      if ResInd2 == ResInd: return i
    raise IndexError, "Could not find atoms in res %d." % ResInd
  for a,b in PairList:
    if DistMethod == 0:
      Atoms1 = GetInd(AtomDat, a, ["CA", "CH3", "N"])
      Atoms2 = GetInd(AtomDat, b, ["CA", "CH3", "N"])
    elif DistMethod == 1:
      Atoms1 = GetInd(AtomDat, a, ["CB", "CA", "CH3", "N"])
      Atoms2 = GetInd(AtomDat, b, ["CB", "CA", "CH3", "N"])
    else:
      Atoms1 = [i for i in range(0,NAtom) if AtomRes[i]==a]
      Atoms2 = [i for i in range(0,NAtom) if AtomRes[i]==b]
    PairAtoms.append((Atoms1, Atoms2))
  return PairAtoms


def GetSaltPairAtoms(PrmtopFile, MinCO = 1):
  """Makes a list of all the atoms involved in salt-forming atom pairs.
Returns PairAtoms, PairList.
* PrmtopFile: string, path of prmtop file
* PairAtoms: list of pairs of lists of atom numbers
* PairList: list of residue number pairs
* MinCO: minimum contact order to consider"""
  AtomNames = coords.GetPrmtopAtomNames(PrmtopFile)
  AtomRes = coords.GetPrmtopAtomRes(PrmtopFile)
  Seq = coords.GetPrmtopSeq(PrmtopFile)
  NAtom, NRes = len(AtomNames), len(Seq)
  #make a list of (residue num, atom name)
  AtomDat = [(AtomRes[i], AtomNames[i].strip()) for i in range(0,NAtom)]
  PairAtoms = []
  PairList = []
  for i in range(0, NRes):
    if not Seq[i] in SaltAtoms: continue
    Charge1, Names1 = SaltAtoms[Seq[i]]
    for j in range(i+MinCO, NRes):
      if not Seq[j] in SaltAtoms: continue
      Charge2, Names2 = SaltAtoms[Seq[j]]
      #check opp charge
      if not Charge1*Charge2 < 0: continue
      for Name1 in Names1:
        Atom1 = AtomDat.index((i, Name1))
        for Name2 in Names2:
          Atom2 = AtomDat.index((j, Name2))
          PairAtoms.append((Atom1, Atom2))
          PairList.append((i,j))
  return PairAtoms, PairList


def GetPhobicPairAtoms(PrmtopFile, DistMethod, MinCO = 3):
  """Makes a list of all the atoms involved in hydrophobic residue pairs.
Returns PairAtoms, PairList.
* PrmtopFile: string, path of prmtop file
* DistMethod: int, 0 = CAs, 1 = CBs, 2 = residue COM
* PairAtoms: list of pairs of lists of atom numbers
* PairList: list of residue number pairs
* MinCO: minimum contact order to consider"""
  Seq = coords.GetPrmtopSeq(PrmtopFile)
  NRes = len(Seq)
  PairList = []
  for i in range(0, NRes):
    for j in range(i+MinCO, NRes):
      if sequence.HydrophobicPair(Seq[i], Seq[j]):
        PairList.append((i,j))
  PairAtoms = GetResPairAtoms(PrmtopFile, PairList, DistMethod)
  return PairAtoms, PairList
  

def SaveTrjDists(TrjFile, PrmtopFile, OutputPath, PairAtoms,
  PairLabels = None, NSkip = 0, NRead = None, NStride = 1,
  Prefix = "", Verbose = VerboseDflt):
  """Saves distances from a trajectory to a gzipped file.
* TrjFile: string, path of trajectory file
* PrmtopFile: string, path of prmtop file
* OutputPath: string, path to save distance file(s)
* PairAtoms: list of pairs of lists of atom numbers
* PairLabels: list of labels for each distance pair
* NSkip: number of configurations to skip
* NRead: maximum number of configurations to read (default is all)
* NStride: stride between configuration frames (default is 1)
* Prefix: string, prefix to add to the name of the distance file
* Verbose: boolean, display verbose messages?"""
  #make sure any pairs are specified
  if len(PairAtoms) == 0:
    return
  #convert to groups
  GroupAtoms, PairGroups = MakeAtomGroups(PairAtoms)
  NGroup = len(GroupAtoms)
  #make trj object
  Trj = coords.TrjClass(TrjFile, PrmtopFile, NSkip = NSkip, NRead = NRead, NStride = NStride)
  #check labels
  if PairLabels is None:
    PairLabels = ["pair%d" % i for i in range(0,len(PairAtoms))]
  #open the output file output the header
  ColHead = ["index".ljust(BlockLen)] + [s.ljust(BlockLen) for s in PairLabels]
  fn = os.path.join(OutputPath, Prefix + "dist.txt.gz")
  f = gzip.GzipFile(fn, "w")
  f.write("".join(ColHead) + "\n")
  #initialize variables
  GroupPos = zeros((NGroup,3), float)
  #now parse the trajectory
  for Pos in Trj:
    #calculate the centroid of group coordinates
    GroupPos.fill(0.)
    for i in range(0,NGroup):
      Group = GroupAtoms[i]
      for AtomNum in Group:
        GroupPos[i,:] += Pos[AtomNum,:]
      GroupPos[i,:] /= float(len(Group))
    #calculate the distances
    DistList = []
    for Group1, Group2 in PairGroups:
      Dist = sqrt( sum((GroupPos[Group1,:] - GroupPos[Group2,:])**2) )
      DistList.append(Dist)
    #output to file
    s = ("%d" % Trj.Index).ljust(BlockLen-1) + " "
    s += " ".join([(DistFmt % d).ljust(BlockLen-1) for d in DistList])
    s += "\n"
    f.write(s)
  f.close()
  Trj.Close()

def GetTrjList(DataPath, ReplicaInd = None):
  """Returns the mdtrj files for a replica exchange simulation."""
  if ReplicaInd is None:
    TrjList = glob.glob(os.path.join(DataPath, '*mdtrj.crd')) + \
              glob.glob(os.path.join(DataPath, '*mdtrj.crd.gz'))
    TrjList.sort()
    ReplicaInd = range(len(TrjList))
    if len(TrjList) == 0:
      raise IndexError, "No trajectories found in %s" % DataPath
  else:
    TrjList = []
    for i in ReplicaInd:
      fn = os.path.join(DataPath, '%d.mdtrj.crd' % i)
      if os.path.isfile(fn):
        TrjList.append(fn)
        continue
      fn += ".gz"
      if os.path.isfile(fn):
        TrjList.append(fn)
      else:
        raise IOError, "Trajectory %d in %s not found" % (i, DataPath)
  return TrjList, ReplicaInd

def SaveAllTrjDists(DataPath, OutputPath, PairAtoms,
  PairLabels = None, ReplicaInd = None, NSkip = 0, NRead = None, NStride = 1,
  Prefix = "", Verbose = VerboseDflt):
  """Saves distances from all trajectories in a path to gzipped files.
* DataPath: path with trajectories and prmtop files
* OutputPath: string, path to save distance files
* PairAtoms: list of pairs of lists of atom numbers
* PairLabels: list of labels for each distance pair
* ReplicaInd: indices of states to be considered (default is all)
* NSkip: number of configurations to skip (default 0)
* NRead: maximum number of configurations to read (default is all)
* NStride: stride between configuration frames (default is 1)
* Prefix: will be added after the number, eg "0.prefixdist.txt.gz"
* Verbose: boolean, display verbose messages?"""
  #check path
  if not os.path.isdir(OutputPath): os.mkdir(OutputPath)
  #find all trajectories
  TrjList, ReplicaInd = GetTrjList(DataPath, ReplicaInd)  
  #for each trajectory
  for TrjFile in TrjList:
    #get trj prefix
    if Verbose: print "Processing trajectory %s" % TrjFile
    TrjPrefix = os.path.basename(TrjFile).replace("mdtrj.crd", "").replace(".gz", "").strip()
    PrmtopFile =  os.path.join(DataPath, TrjPrefix + "prmtop.parm7")
    SaveTrjDists(TrjFile, PrmtopFile, OutputPath, PairAtoms,
      PairLabels = PairLabels, NSkip = NSkip, NRead = NRead, NStride = NStride,
      Prefix = TrjPrefix + Prefix, Verbose = Verbose)   
  if Verbose: print "Done processing trajectory distances"
  

def SaveTrjResDists(TrjFile, PrmtopFile, OutputPath, PairList,
  NSkip = 0, NRead = None, NStride = 1,
  Prefix = "", DistMethod = 0, StartRes = 0, Verbose = VerboseDflt):
  """Saves distances from a trajectory to a gzipped file.
* TrjFile: string, path of trajectory file
* PrmtopFile: string, path of prmtop file
* OutputPath: string, path to save distance file(s)
* PairList: list of (a,b) tuples, residue pairs to monitor
* NSkip: number of configurations to skip
* NRead: maximum number of configurations to read (default is all)
* NStride: stride between configuration frames (default is 1)
* Prefix: string, prefix to add to the name of the distance file
* DistMethod: int, 0 = CAs, 1 = CBs, 2 = residue COM
* StartRes: int, the named index of the first residue in the coord
            files, such that a PairList entry of (3,4) with StartRes
            equal to 2 means a pair between the 2nd and 3rd residues
            in the trajectory
* Verbose: boolean, display verbose messages?"""
  #make sure any pairs are specified
  if len(PairList) == 0:
    return
  #find numbers using a different starting residue
  NewPairList = [(a-StartRes, b-StartRes) for (a,b) in PairList]
  #make atom list
  PairAtoms = GetResPairAtoms(PrmtopFile, NewPairList, DistMethod)
  #make labels
  PairLabels = ["%d,%d " % (a+1,b+1) for (a,b) in PairList]
  #save the distances
  SaveTrjDists(TrjFile, PrmtopFile, OutputPath, PairAtoms,
    PairLabels = PairLabels, NSkip = NSkip, NRead = NRead, NStride = NStride,
    Prefix = Prefix, Verbose = Verbose)

    
def SaveAllTrjResDists(DataPath, OutputPath, PairList,
  ReplicaInd = None, NSkip = 0, NRead = None, NStride = 1,
  DistMethod = 0, StartRes = 0, Prefix = "", Verbose = VerboseDflt):
  """Saves distances from all trajectories in a path to gzipped files.
* DataPath: path with trajectories and prmtop files
* OutputPath: string, path to save distance files
* PairList: list of (a,b) tuples, residue pairs to monitor
* ReplicaInd: indices of replicas to be considered (default is all)
* NSkip: number of configurations to skip (default 0)
* NRead: maximum number of configurations to read (default is all)
* NStride: stride between configuration frames (default is 1)
* DistMethod: int, 0 = CAs, 1 = CBs, 2 = residue COM
* StartRes: int, the named index of the first residue in the coord
            files, such that a PairList entry of (3,4) with StartRes
            equal to 2 means a pair between the 2nd and 3rd residues
            in the trajectory
* Prefix: will be added after the number, eg "0.prefixdist.txt.gz"
* Verbose: boolean, display verbose messages?"""
  #check path
  if not os.path.isdir(OutputPath): os.mkdir(OutputPath)
  #find all trajectories
  TrjList, ReplicaInd = GetTrjList(DataPath, ReplicaInd)       
  #for each trajectory
  for TrjFile in TrjList:
    #get trj prefix
    if Verbose: print "Processing trajectory %s" % TrjFile
    TrjPrefix = os.path.basename(TrjFile).replace("mdtrj.crd", "").replace(".gz", "").strip()
    PrmtopFile =  os.path.join(DataPath, TrjPrefix + "prmtop.parm7")
    SaveTrjResDists(TrjFile, PrmtopFile, OutputPath, PairList,
      NSkip = NSkip, NRead = NRead, NStride = NStride,
      Prefix = TrjPrefix + Prefix, DistMethod = DistMethod,
      StartRes = StartRes, Verbose = Verbose)   
  if Verbose: print "Done processing trajectory distances"
  

def DeleteAllTrjDists(Path):
  "Deletes all distance files in a given Path."
  fl = glob.glob(os.path.join(Path, "*dist.txt.gz"))
  for f in fl:
    os.remove(f)



