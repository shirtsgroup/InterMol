#!/usr/bin/env python

#LAST MODIFIED: 05-11-07

Usage = """Runs analysis on a zam fragment.

Usage     : zamanal.py FRAGPATH [OPTIONS] [STARTRES STOPRES]

ANALTIME  : time in ps at end of trajectories to use for analysis
FRAGPATH  : path of fragment (e.g., "r1-8a")
FRAGPAIRS : either "all" or "hydrophobic"
OPTIONS   : "--nocleanup" to keep the distance files
            "--nopunch" to skip punch analysis
            "--targettemp=X" to change target temperature
            "--phobicpairs" to just do hydrophobic residue pairs
            "--analtime=X" to change analysis time (default 1ns)
            "--skiptime=X" to change time skipped (default all but analtime)
            "--keepcaps" to include caps in sequence and residue numbers
STARTRES  : number of the starting residue (not including cap)
STOPRES   : number of the stopping residue (not including cap)
"""
import sys
if __name__ == "__main__" and len(sys.argv) == 1:
  print Usage
  sys.exit()



from numpy import *
import os, math, cPickle, time, sys, gzip, glob, shutil, copy
import rexsock, rexanalysis, mdsim, mdtrj, coords, rmsd, protein, sequence
import punchanalysis, mesoanalysis, scoreconf
import zamdata
import scripttools


#analysis variables
DistMethod = 1      #contact distances calc between 0 = CA, 1 = CB, 2 = residue centroid
PmfNBin = 11        #num of pmf bins
PmfBinMin = 4.      #minimum pmf distance
PmfBinMax = 15.     #maximum pmf distance
ContactCut = 8.0    #distance defining a contact cutoff
TargetTemp = 270.   #target temperature for calculations
DeltaTemp = 200.    #maximum temperature range to consider in WHAM calculations
DeltaTempUmbr = 200.#maximum temperature range to consider in WHAM when umbrella sampling
RemoveRest = True   #remove restraints using wham?
ChainPmfNBin = 46   #number of bins in chain pmf
ChainPmfMin = 4.    #minimum chain pmf distance
ChainPmfMax = 50    #maximum chain pmf distance
UmbrNBin = 50       #num of pmf bins for umbrella pmfs
UmbrWindow = 10.    #num of angstroms on each side of umbrella range to do min, max for pmf
UmbrDistMethod = 1  #umbrella distances calc between 0 = CA, 1 = CB, 2 = residue centroid

#clustering parameters
MaxCluster = 10     #maximum number of cluster configurations
ClustRmsdTol = 2.0  #tolerance rmsd between clusters
ClustMaxIter = 10   #maximum number of iterations

#punch analysis
PunchCritProb = 0.5 #probability threshold for punch calculation
PunchCritCoop = 0.4 #cooperativity threshold for punch calculation

#score parameters
PmfCutDist = 9.0     #distance in angstroms for taking pmf difference
ScoreSlope = -14.0   #slope for critical scores in (Cut,PMFDiff) plane
ScoreInt = 11.0      #intercept for critical scores in (Cut,PMFDiff) plane

#bestn analysis
BestnMinECO = 3      #minimum ECO for generating bestn pairs
BestnN = 3           #number of contact pairs to choose

#energy analysis
EneTerms = ["TEMP","ETOT","EPOT","EREST","EVDW","EEL","ESURF"]
                     #terms to average from the energy function

#score analysis
ScoreNStride = 20    #stride in frames for calculating various scores

#timeouts
AnalTimeout = None #timeout in seconds for the analysis routines
SnapTimeout = None #timeout in seconds for the snapshot routines



def GetAnalFrames(Frag):
  """Returns NFrameSkip, NFrameRead for analysis frames."""
  NFrameTot = int(Frag.RexTime / Frag.FrameTime + 0.5)
  NFrameRead = min(int(Frag.AnalTime / Frag.FrameTime + 0.5), NFrameTot)
  if Frag.SkipTime < 0:
    NFrameSkip = NFrameTot - NFrameRead
  else:
    NFrameSkip = int(Frag.SkipTime / Frag.FrameTime + 0.5)
    NFrameSkip = min(NFrameSkip, NFrameTot - NFrameRead)
  return NFrameSkip, NFrameRead


def GetTargetReplica(Frag):
  """Returns the index of the replica whose temperature is
  closest to but less than the target analysis temperature."""
  DataPath = os.path.join(Frag.BasePath, "data/")
  Temps = file(os.path.join(DataPath, "temps.txt")).read().split()
  Temps = array([float(x) for x in Temps], float)
  Temps = where(Temps <= TargetTemp, Temps, min(Temps))
  TarRep = argmax(Temps)
  return TarRep, Temps[TarRep]


def RunAnal(Frag):
  "Calculates residue contact distances and runs WHAM analysis."
  #setup vars
  DataPath = os.path.join(Frag.BasePath, "data/")
  AnalPath = os.path.join(Frag.BasePath, "anal/")
  if Frag.CapN:
    ResShift = 1
  else:
    ResShift = 0

  #make the analysis path
  if not os.path.isdir(AnalPath):
    os.mkdir(AnalPath)

  #get number of frames to analyze
  NFrameSkip, NFrameRead = GetAnalFrames(Frag)

  #get pairlist  
  Frag.PairList = Frag.GetPairList()

  #find the replica closest to the target temperature
  TarRep, TarTemp = GetTargetReplica(Frag)
  print "Using replica %d with temperature %.2f for analysis." % (TarRep, TarTemp)

  #run wham on all replicas
  ReplicaInd, TempK, NFrame, NFrameSkip, NFrameRead, FK, LogwKN = \
    rexanalysis.RunWham(DataPath, NFrameSkip, NFrameRead,
    TargetTemp = TargetTemp, RemoveRest = RemoveRest, Verbose = True)
  #output results
  s = "WHAM FREE ENERGIES\nindex, T, F/kT\n" + \
      "\n".join(["%-3d %-8.2f %-11.4f" % (ind, TempK[i], FK[i])
                 for (i, ind) in enumerate(ReplicaInd)])
  file(os.path.join(AnalPath, "whamresults.txt"), "w").write(s)

  print "Skipping %d frames and reading %d frames" % (NFrameSkip, NFrameRead)  

  #find the maximum temperature difference
  if len(Frag.UmbrRest) > 0:
    dT = DeltaTempUmbr
  else:
    dT = DeltaTemp
  #make sure target temperature is there at least
  dT = max(abs(TargetTemp - TarTemp)*1.001, dT) 
  #get rid of other temperatures
  MinTemp, MaxTemp = TargetTemp - dT, TargetTemp + dT
  Ind = [i for (i,j) in enumerate(ReplicaInd) if TempK[i] < MaxTemp and TempK[i] >= MinTemp]
  ReplicaInd, TempK, FK = ReplicaInd[Ind], TempK[Ind], FK[Ind]
  LogwKN = LogwKN.take(Ind, axis=0)
  
  #calculate pmfs and cutoff probs for res-res pairs
  Frag.PmfDist, Frag.Pmf, Frag.CutProb  = rexanalysis.CalcDistProbs(DataPath, AnalPath,
    NFrameSkip, NFrameRead, Frag.PairList, DistMethod,
    PmfNBin, PmfBinMin, PmfBinMax, ContactCut,
    TargetTemp = TargetTemp, ReplicaInd = ReplicaInd, NFrame = NFrame, LogwKN = LogwKN,
    StartRes = Frag.StartRes-ResShift, Verbose = True, KeepDist = True)
  
  #calculate chain pmfs
  p = protein.ProteinClass(Pdb = os.path.join(DataPath, "0.current.pdb"))
  NChain = len(p.Chains)
  if NChain > 1:
    AtomGroups = [range(p.ChainAtomNums[i], p.ChainAtomNums[i+1]) for i in range(NChain)]
    PairAtoms = [[AtomGroups[i], AtomGroups[j]] for i in range(NChain) for j in range(i+1,NChain)]
    PairLabels = ["%d-%d" % (i+1, j+1) for i in range(NChain) for j in range(i+1,NChain)]
    Pmf, PmfDist = rexanalysis.CalcAtomDistProbs(DataPath, AnalPath, NFrameSkip, NFrameRead,
      PairAtoms, ChainPmfNBin, ChainPmfMin, ChainPmfMax, TargetTemp = TargetTemp,
      ReplicaInd = ReplicaInd, NFrame = NFrame, LogwKN = LogwKN, PairLabels = PairLabels,
      Prefix = "chain", Verbose = True, KeepDist = True)

  #calculate umbrella pmfs
  if len(Frag.UmbrRest) > 0:
    rMin = max(min(Frag.UmbrDist) - UmbrWindow, 0.)
    rMax = max(Frag.UmbrDist) + UmbrWindow
    rCut = average(Frag.UmbrDist)
    PmfDist, Pmf, CutProb  = rexanalysis.CalcDistProbs(DataPath, AnalPath,
      NFrameSkip, NFrameRead, Frag.UmbrRest, UmbrDistMethod, UmbrNBin,
      rMin, rMax, rCut, TargetTemp = TargetTemp, ReplicaInd = ReplicaInd,
      NFrame = NFrame, LogwKN = LogwKN, StartRes = Frag.StartRes-ResShift,
      Prefix = "umbr", Verbose = True, KeepDist = True)  

  #calculate the score
  CalcPmfScores(Frag)

  #run the energy and pair analyses
  RunEneAnal(Frag)
  RunScoreAnal(Frag)

  #run the swapped restraint analysis
  RunSwapRestAnal(Frag)

  #run the mesostring analysis for the lowest temperature
  TrjFile = os.path.join(DataPath, "%d.mdtrj.crd.gz" % TarRep)
  PrmtopFile = os.path.join(DataPath, "%d.prmtop.parm7" % TarRep)
  ConfMeso, MesoPop, MesoEntropy = mesoanalysis.RunAnalTrj(TrjFile, PrmtopFile,
    OutputPath = AnalPath, NSkip = NFrameSkip, NRead = NFrameRead)
  #set fragment variables
  NTot = float(len(ConfMeso))
  Frag.MesoEntropy = MesoEntropy
  Frag.MesoGroundPop = float(MesoPop[0][1]) / NTot
  #set the structured score;
  #this is just the negative, rezeroed mesostring entropy
  Frag.StructScore = math.log(NTot) - MesoEntropy

  #run the punch analysis
  if Frag.RunPunch:
    PunchPath = os.path.join(AnalPath, "punch/")
    Frag.PunchStabTot, Frag.PunchCoopTot, Frag.PunchStab, Frag.PunchCoop = \
      punchanalysis.RunAnal(AnalPath, PunchPath, Frag.PairList,
      ContactCut, PunchCritProb, PunchCritCoop, NFrameRead,
      TargetTemp, ReplicaInd, LogwKN)

    #run the bestn analysis using punch stab + cutprob/10
    Scores = [Frag.PunchStab[i] * 0.5 + Frag.CutProb[i] * 0.5
              for i in range(0, len(Frag.PairList))]
    RunBestn(Frag, Scores)
  else:
    #run the bestn analysis using cutoff probabilities
    RunBestn(Frag, Frag.CutProb)    

  #perform a clustering analysis
  #first delete old clusters
  for fn in glob.glob(os.path.join(AnalPath, "clust*.pdb")): os.remove(fn)
  TrjFn = os.path.join(Frag.BasePath, "data/%d.mdtrj.crd.gz" % TarRep)
  PrmtopFn = os.path.join(Frag.BasePath, "data/%d.prmtop.parm7" % TarRep)
  Trj = coords.TrjClass(TrjFn, PrmtopFn, Mask = coords.CrdBackboneMask,
                        NSkip = NFrameSkip, NRead = NFrameRead)
  Pos, ClustNum, ClustPop, ConfRmsd, ClustRmsd = rmsd.ClusterMSS(Trj,
    ClustRmsdTol, MaxIter = ClustMaxIter, MaxCluster = MaxCluster, Method = 0)
  #save the clusterinfo
  Prefix = os.path.join(Frag.BasePath, "anal/clust")
  rmsd.SaveClustResults(Pos, ClustNum, ClustPop, ConfRmsd, ClustRmsd,
    Prefix = Prefix, ConfIndices = Trj.GetPastIndices(), Verbose = False)
  Frag.ClustPop = ClustPop
  #save the pdb files
  Pdbs = []
  for i in range(0,len(Pos)):
    pc = 100. * float(ClustPop[i]) / float(len(ClustNum))
    fn = os.path.join(AnalPath, "clust%02d_%.0fpc.pdb" % (i+1, pc))
    Pdbs.append(fn)
    coords.SavePdbCoordsAmb(Pos[i], PrmtopFn, fn)
  #score the clusters
  PdbScores, AvgScores, Summary = scoreconf.ScorePdbs(Pdbs, BaseName=True)
  file(os.path.join(AnalPath, "clustmetrics.txt"), "w").write(Summary)

  #save all the analysis data
  Frag.SaveAnalData()
  

def CalcPmfScores(Frag):
  "Calculates scores from pmfs and cutoff probs and saves to file."
  #record scores in a score file
  c = 10
  ScoreFile = file(os.path.join(Frag.BasePath, "anal/scores.txt"), "w")
  ScoreFile.write("contact".ljust(c) + " " + "cont_prob".ljust(c) + " "
                  + "delta_pmf".ljust(c) + " " + "score".ljust(c) + "\n")
  Frag.PmfScore = []
  for i in range(0, len(Frag.PairList)):
    (a,b) = Frag.PairList[i]
    #get a score
    DeltaPmf, Score = ContScore(Frag.CutProb[i], Frag.Pmf[i], Frag.PmfDist)
    #write score to a file
    ContName = "%d,%d" % (a+1, b+1)
    ContProbStr = "%.2f" % Frag.CutProb[i]
    DeltaPmfStr = "%.2f" % DeltaPmf
    ScoreStr = "%.2f" % Score
    ScoreFile.write(ContName.ljust(c) + " " + ContProbStr.ljust(c) + " "
                    + DeltaPmfStr.ljust(c) + " " + ScoreStr.ljust(c) + "\n")
    Frag.PmfScore.append(Score)
  #close the score file
  ScoreFile.close()
  
def ContScore(CutProb, Pmf, PmfDist):
  "Assigns a score to a contact based on analysis data."
  #find the cutoff index
  dr = PmfDist[1] - PmfDist[0]
  Ind = int((PmfCutDist - PmfDist[0])/dr + 0.5) + 1
  #disallow the last bin, which contains extra entries from
  #outside of the calculated pmf range
  LastInd = len(Pmf) - 1
  Ind = max(min(Ind, LastInd-1), 0)
  #find the values of the pmf minima before and after the cutoff,
  #take their difference, and return it as the score
  PmfDiff = min(Pmf[Ind:LastInd]) - min(Pmf[:Ind])
  #calculate the relative critical amount
  CutCrit = (PmfDiff - ScoreInt) / ScoreSlope
  return PmfDiff, CutProb - CutCrit


def RunSwapRestAnal(Frag):
  "Runs analysis on sorted restraints."
  if len(Frag.SwapRest) == 0 or Frag.RexTimeSwapRest == 0.: return
  print "Running swap restraint analysis."
  #get the restraint number in each replica
  RepRest = [i % len(Frag.SwapRest) for i in range(Frag.NRep)]
  #get cycles
  NCycleTot = int(Frag.RexTime / Frag.CycleTime + 0.5)
  NCycleRead = min(int(Frag.AnalTime / Frag.CycleTime + 0.5), NCycleTot)
  if Frag.SkipTime < 0:
    NCycleSkip = NCycleTot - NCycleRead
  else:
    NCycleSkip = int(Frag.SkipTime / Frag.CycleTime + 0.5)
    NCycleSkip = min(NCycleSkip, NCycleTot - NCycleRead)
  NCycleSwap = min(int(Frag.RexTimeSwapRest / Frag.CycleTime + 0.5), NCycleTot)
  #sort through the replica numbers and calculate the
  #average state number
  NRest = len(Frag.SwapRest)
  RestStateNum1 = zeros(NRest, float)
  Counts1 = ones(NRest, float) * 0.0001
  RestStateNum2 = zeros(NRest, float)
  Counts2 = ones(NRest, float) * 0.0001
  f = gzip.GzipFile(os.path.join(Frag.BasePath, "data/statenum.txt.gz"), "r")
  i = 0
  while True:
    s = f.readline()
    if len(s) == 0: break
    i += 1
    StateNum = [int(x) for x in s.split()]
    if i <= NCycleSwap:
      for j in range(Frag.NRep):
        RestStateNum1[RepRest[j]] += StateNum[j]
        Counts1[RepRest[j]] += 1.
    if i > NCycleSkip:
      for j in range(Frag.NRep):
        RestStateNum2[RepRest[j]] += StateNum[j]
        Counts2[RepRest[j]] += 1.
  f.close()
  #average
  RestStateNum1 /= Counts1 * Frag.NRep
  RestStateNum2 /= Counts2 * Frag.NRep
  #write results
  s = "SWAPPED RESTRAINTS RESULTS\n"
  s += "restraint, frac position during swaps, frac position during anal\n"
  for (i, (a,b)) in enumerate(Frag.SwapRest):
    #correct for a, b are numpy types -- this should be fixed in other parts of the code at some point
    a, b = int(a), int(b)
    s += "%s%d-%s%d %-01.4f %-01.4f\n" % (Frag.Seq[a - Frag.StartRes], a+1,
                                          Frag.Seq[b - Frag.StartRes], b+1,
                                          RestStateNum1[i],
                                          RestStateNum2[i])
  file(os.path.join(Frag.BasePath, "anal/swapresults.txt"), "w").write(s)  
  

def BestnQuick(Frag, ScorePairList, N, MaxKeepLev = 5, MaxKeepGlob = 100,
               CurDepth = 0, CurScore = 0, CurList = []):
  """"Runs an approximate best-n algorithm by just keeping the top MaxKeepLev
scores at each depth."""
  def PairECO(a1, b1, a2, b2):
    "Returns a simple metric of ECO."
    a1, b1 = min(a1, b1), max(a1, b1)
    a2, b2 = min(a2, b2), max(a2, b2)
    return abs(a1-a2) + abs(b1-b2) + 1
  if CurDepth == N: return [(CurScore, CurList)]
  #first prune out any low eco contacts
  if CurDepth == 0:
    ScorePairList = [(s,(a,b)) for (s,(a,b)) in ScorePairList
                      if not Frag.ResWithinECO(a, b, OtherRest = CurList,
                      MaxECO = BestnMinECO - 1, MaxDepth = 1)]
    ScorePairList.sort(reverse = True)
  if len(ScorePairList) == 0: return [(CurScore, CurList)]
  #now recurse the top ones
  RetList = []
  i = 0
  for (s, (a,b)) in ScorePairList[:MaxKeepLev]:
    i += 1
    #remove any low-eco contacts
    NewScorePairList = [(s2, (a2,b2)) for (s2, (a2,b2)) in ScorePairList[i:]
                        if PairECO(a, b, a2, b2) >= BestnMinECO]
    RetList.extend(BestnQuick(Frag, NewScorePairList, N, MaxKeepLev, MaxKeepGlob,
                              CurDepth + 1, CurScore + s, CurList + [(a,b)]))
  RetList.sort(reverse = True)
  RetList = RetList[:MaxKeepGlob]
  return RetList

def RunBestn(Frag, Scores):
  "Runs the bestn analysis for a fragment."
  print "Running bestn analysis"
  
  def WriteResults(Frag, FileName, BestnList):
    def Rep(Group):
      return ", ".join(["%s%d-%s%d" % (Frag.Seq[a - Frag.StartRes], a+1,
        Frag.Seq[b - Frag.StartRes], b+1) for (a,b) in Group])
    #make the output string and write it
    if len(BestnList) == 0:
      s = "no contacts found\n"
    else:
      s = "BEST SCORE AND GROUP"
      s += "\n%.4f  %s" % (BestnList[0][0], Rep(BestnList[0][1]))
      s += "\n\nALL SCORES AND GROUPS"
      for (GroupScore, GroupList) in BestnList:
        s += "\n%.4f  %s" % (GroupScore, Rep(GroupList))
      s += "\n"
    file(FileName, "w").write(s)

  #get the path
  AnalPath = os.path.join(Frag.BasePath, "anal/")
  
  NPair = len(Frag.PairList)
  
  #first do all pairs
  if Frag.AllPairs:
    ScorePairList = [(Scores[i], Frag.PairList[i]) for i in range(NPair)]
    BestnList = BestnQuick(Frag, ScorePairList, BestnN)
    WriteResults(Frag, os.path.join(AnalPath, "bestn-all.txt"), BestnList)
    
  #do hydrophobics
  ScorePairList = [(Scores[i], Frag.PairList[i]) for i in range(NPair)]
  ScorePairList = [(s, (a,b)) for (s, (a,b)) in ScorePairList
                   if sequence.HydrophobicPair(Frag.Seq[a-Frag.StartRes],
                                               Frag.Seq[b-Frag.StartRes])]
  BestnList = BestnQuick(Frag, ScorePairList, BestnN)
  WriteResults(Frag, os.path.join(AnalPath, "bestn-phobics.txt"), BestnList)

  #do all non-charged pairs
  if Frag.AllPairs:
    ScorePairList = [(Scores[i], Frag.PairList[i]) for i in range(NPair)]
    ScorePairList = [(s, (a,b)) for (s, (a,b)) in ScorePairList
                     if not sequence.Charged(Frag.Seq[a-Frag.StartRes])
                     and not sequence.Charged(Frag.Seq[b-Frag.StartRes])]
    BestnList = BestnQuick(Frag, ScorePairList, BestnN)
    WriteResults(Frag, os.path.join(AnalPath, "bestn-noncharged.txt"), BestnList)
    
  #set the bestn contacts using the non-salt pairs as a result
  ind = [Frag.PairList.index((a,b)) for (a,b) in BestnList[0][1]]
  Frag.BestnCont = []
  for i in ind:
    (a,b) = Frag.PairList[i]
    Frag.BestnCont.append((Scores[i], a, b))


def RunEneAnal(Frag):
  "Calculates average energies."
  print "Running energy analysis"
  #count the replicas
  DataPath = os.path.join(Frag.BasePath, "data/")
  AnalPath = os.path.join(Frag.BasePath, "anal/")
  Temps = file(os.path.join(DataPath, "temps.txt"), "r").read().split()
  Temps = [float(x) for x in Temps]
  NRep = len(Temps)
  NFrameSkip, NFrameRead = GetAnalFrames(Frag)
  #average energies
  s = "AVERAGE ENERGIES FROM ANALYSIS FRAMES"
  s += "\n%-3s %-8s " % ("Rep", "TempSet")+ " ".join(["%-10s" % t for t in EneTerms])
  for i in range(0, NRep):
    print "Getting trajectory %d energies" % i
    Ret = mdsim.GetConcatData(str(i), DataPath, Vars = EneTerms)
    #convert to a string for output
    s += "\n%-3d %-8.2f " % (i, Temps[i])
    s += " ".join(["%-10.2f" % mean(Ret[t][NFrameSkip:]) for t in EneTerms])
  s += "\n"
  file(os.path.join(AnalPath, "eneresults.txt"), "w").write(s)
  

def RunScoreAnal(Frag):  
  "Calculates all sorts of metrics."
  print "Running hydrophobic, salt pair, and other metric analysis"
  #count the replicas
  DataPath = os.path.join(Frag.BasePath, "data/")
  AnalPath = os.path.join(Frag.BasePath, "anal/")
  Temps = file(os.path.join(DataPath, "temps.txt"), "r").read().split()
  Temps = [float(x) for x in Temps]
  NRep = len(Temps)
  NFrameSkip, NFrameRead = GetAnalFrames(Frag)
  s = "AVERAGE METRICS FOR ANALYSIS FRAMES\n"
  s += "%-3s %-8s %s\n" % ("Rep", "TempSet", scoreconf.HeadStr())
  #get the average scores for each run
  for i in range(0, NRep):
    print "Examining trajectory %d metrics" % i
    TrjFile = os.path.join(DataPath, "%d.mdtrj.crd.gz" % i)
    PrmtopFile = os.path.join(DataPath, "%d.prmtop.parm7" % i)
    #make trj class
    Trj = coords.TrjClass(TrjFile, PrmtopFile, NSkip = NFrameSkip,
                          NRead = -1, NStride = ScoreNStride)
    #run the scoring
    TrjScores, AvgScores, Summary = scoreconf.ScoreTrj(Trj)
    s += "%-3d %-8.2f %s\n" % (i, Temps[i], scoreconf.ScoreStr(AvgScores))
    #cleanup
    Trj.Close()
  file(os.path.join(AnalPath, "metricresults.txt"), "w").write(s)



def ExtractSnaps(Frag):
  "Gets representative snapshots from a REX simulation."
  #first make the path
  AnalPath = os.path.join(Frag.BasePath, "anal/")
  SnapPath = os.path.join(Frag.BasePath, "conf/")
  if not os.path.isdir(SnapPath): os.mkdir(SnapPath)
  #get the cluster pdbs
  Pdbs = glob.glob(os.path.join(AnalPath, "clust*.pdb"))
  #make the protein classes and get their ss
  ps = [protein.ProteinClass(Pdb = x).Decap() for x in Pdbs]
  ps = [(p, p.SecondaryStructure()) for p in ps]
  #make a weights list
  Weights = [Frag.ClustPop[i] for (i,p) in enumerate(ps)] 
  #remove redundant secondary structure
  def SSInSS(SS1, SS2):
    for i in range(len(SS1)):
      if SS1[i] != SS2[i] and SS1[i] != " ":
        return False
    return True
  if Frag.SnapTrimSS:
    ps2, Weights2 = [], []
    for (i, (p1, ss1)) in enumerate(ps):
      Test = [SSInSS(ss1, ss2) for (p2, ss2) in ps if not p1 is p2]
      if Test.count(True) == 0:
        ps2.append((p1, ss1))
        Weights2.append(Weights[i])
    ps = ps2
    Weights = Weights2
  #trim by loose ends
  ps = [(p, ss, Frag.StartRes, Frag.StopRes) for (p, ss) in ps]
  if Frag.SnapTrimEnds:
    ps2 = []
    for (p, ss, StartRes, StopRes) in ps:
      a, b = 0, len(p)
      while a < len(p) and ss[a] == " ":
        a += 1
      while b > 0 and ss[b-1] == " ":
        b -= 1
      if b > a: ps2.append((p[a:b], ss, StartRes+a, StartRes+b-1))
    ps = ps2
  #save the pdbs
  Frag.SnapPdbs, Frag.SnapResRange = [], []
  ps = [p for (p, ss, StartRes, StopRes) in ps]        
  for i in range(0, len(ps)):
    fn = "snap%02d.pdb" % i
    ps[i].WritePdb(os.path.join(SnapPath, fn))
    Frag.SnapPdbs.append(fn)
    Frag.SnapResRange.append((StartRes, StopRes))
  #save the weights
  Tot = float(sum(Weights))
  Frag.SnapWeight = [float(w) / Tot for w in Weights]


def CleanUp(Frag, ConserveDisk):
  """Removes analysis data from the fragment class after it's used
in order to save memory and disk space."""
  Frag.ClearAnalData()
  #remove distance files
  AnalPath = os.path.join(Frag.BasePath, "anal/")
  mdtrj.DeleteAllTrjDists(AnalPath)


def RunAnalTask(sr, Frag, Task):
  """Runs all the analysis routines as a socketring task."""
  #check to see if we're running
  if not Task is None:
    if Task.Done:
      #replace everything with the output
      Frag.LoadsData(Task.Result)
      sr.RemoveTask(Task)
      return True, None
    else:
      return False, Task
  else:
    #make the task
    s = "import zamanal, zamdata"
    s += "\nFrag = zamdata.FragDataClass()"
    s += "\nFrag.LoadsData(SockData)"
    s += "\nzamanal.RunAnal(Frag)"
    s += "\nSockResult = Frag.DumpsData()"
    Task = sr.AddTask(s, Data = Frag.DumpsData(), Timeout = AnalTimeout)
    return False, Task
    
def RunSnapTask(sr, Frag, Task):
  """Runs snapshot getting routine as a socketring task."""
  #check to see if we're running
  if not Task is None:
    if Task.Done:
      #replace everything with the output
      Frag.LoadsData(Task.Result)
      sr.RemoveTask(Task)
      return True, None
    else:
      return False, Task
  else:
    #make the task
    s = "import zamanal, zamdata"
    s += "\nFrag = zamdata.FragDataClass()"
    s += "\nFrag.LoadsData(SockData)"
    s += "\nzamanal.ExtractSnaps(Frag)"
    s += "\nSockResult = Frag.DumpsData()"
    Task = sr.AddTask(s, Data = Frag.DumpsData(), Timeout = SnapTimeout)
    return False, Task


if __name__ == "__main__":
  #make a dummy fragment so we can run from the command line
  def GetNum(s):
    return int("".join([x for x in s if x in "0123456789"]))
  Args = scripttools.ParseArgs(sys.argv)
  p = Args[1]
  AnalTime = float(Args.get("analtime", 1000))
  SkipTime = float(Args.get("skiptime", -1))
  AllPairs = not "phobicpairs" in Args["FLAGS"]
  RunCleanUp = not "nocleanup" in Args["FLAGS"]
  RunPunch = not "nopunch" in Args["FLAGS"]
  KeepCaps = "keepcaps" in Args["FLAGS"]
  TargetTemp = float(Args.get("targettemp", TargetTemp))
  if Args["NARG"] > 3:
    StartRes = int(sys.argv[2]) - 1
    StopRes = int(sys.argv[3]) - 1
  else:
    LastPath = p.split("/")[-1]
    StartRes = GetNum(LastPath.split("-")[0]) - 1
    StopRes = GetNum(LastPath.split("-")[1]) - 1
  f = os.path.join(p, "data/seq.txt")
  if not os.path.isfile(f):
    print "Could not find sequence for fragment."
    sys.exit()
  #get the sequence and get rid of caps
  Seq = file(f, "r").read().split()
  #make a dummy fragment
  if KeepCaps:
    CapN, CapC = False, False
  else:
    CapN = Seq[0] in ["ACE"]
    CapC = Seq[-1] in ["NME", "NHE"]
    if CapN: Seq = Seq[1:]
    if CapC: Seq = Seq[:-1]
  Frag = zamdata.FragDataClass(StartRes = StartRes, StopRes = StopRes, Seq = Seq)
  Frag.CapN, Frag.CapC = CapN, CapC
  Frag.BasePath = p
  Frag.AnalTime = AnalTime
  Frag.SkipTime = SkipTime
  Frag.AllPairs = AllPairs
  Frag.RunPunch = RunPunch
  #get the frame time
  f = os.path.join(p, "data/rexstats.txt")
  if not os.path.isfile(f):
    print "Could not find rexstats.txt for fragment"
    sys.exit()
  for l in file(f, "r").readlines():
    if l.startswith("Time in ps per frame"):
      Frag.FrameTime = float(l.split(":")[-1])
    elif l.startswith("Elapsed time in ps"):
      Frag.RexTime = float(l.split(":")[-1])
  #run everything
  Frag.DelAnalysis()
  RunAnal(Frag)
  ExtractSnaps(Frag)
  if RunCleanUp: CleanUp(Frag, False)

    