#!/usr/bin/env python

#LAST MODIFIED: 05-19-07

#TODO:

from numpy import *
import sys, os, gzip, time, math, copy, cPickle
import mdsim, mdtrj, coords
try:
  import whamlib
except ImportError:
  print "Could not import whamlib."

#globals
kB= mdsim.kB
TargetTempDflt = 270.
DeltaTempDflt = 1.e10
PmfMax = 1.0e4
VerboseDflt = False


#======= GET INFORMATION ABOUT REPLICAS ========

def GetTrjFiles(DataPath, Ind):
  "Returns the names of the trajectory files."
  fn1 = os.path.join(DataPath, "%d.mdtrj.crd" % Ind)
  fn2 = os.path.join(DataPath, "%d.prmtop.parm7" % Ind)
  if os.path.isfile(fn1):
    return fn1, fn2
  fn1 += ".gz"
  if os.path.isfile(fn1):
    return fn1, fn2
  else:
    raise IOError, "Could not find trajectory %d in %s" % (Ind, DataPath)

def GetStats(DataPath, MinTemp = -1., MaxTemp = 1.e100):
  "Reads temperatures and gets data length from a REX simulation."
  # read temperatures
  TempFile = os.path.join(DataPath, "temps.txt")
  TempList = file(TempFile,'r').read().split()
  TempK = array([float(x) for x in TempList], float)
  ReplicaInd = nonzero(logical_and(TempK < MaxTemp, TempK >= MinTemp))[0]
  TempK = TempK[ReplicaInd]
  BetaK = 1.0 / TempK / kB 
  #get data length of each trajectory
  NFrameList = [mdsim.GetConcatNFrames(str(i), DataPath) for i in ReplicaInd]
  #find the minimum
  NFrame = min(NFrameList)
  #check to see if everything is the same
  if all(NFrameList == NFrame):
    print "Warning: not all trajectories have same number of frames:", NFrameList
  #return data
  return ReplicaInd, TempK, BetaK, NFrame

def CheckFrames(NFrame, NFrameSkip, NFrameRead):
  """Returns new NFrameSkip and NFrameRead calculated to fit
  into the total number of frames."""
  #check number of frames to read
  if NFrameRead < 0:
    #read to the end
    NFrameSkip = min(NFrameSkip, NFrame)
    NFrameRead = NFrame - NFrameSkip
  else:
    #truncate skipped frames if too long
    NFrameRead = min(NFrameRead, NFrame)
    NFrameSkip = min(NFrame, NFrameRead + NFrameSkip) - NFrameRead
  return NFrameSkip, NFrameRead

def GetNReplica(DataPath, ReplicaInd = None):
  "Returns the nuber of replicas to use."
  if ReplicaInd is None:
    #count temperatures
    TempList = file(os.path.join(DataPath, "temps.txt"),'r').read().split()
    NRep = len(TempList)
    ReplicaInd = arange(NRep, dtype=int)
  else:
    NRep = len(ReplicaInd)


#======== GET ENERGIES FROM REPLICAS ========
    
def GetEnergies(DataPath, ReplicaInd, Vars, NFrameSkip, NFrameRead):
  "Returns arrays of energies corresponding to the terms specified."
  NComp, NRep = len(Vars), len(ReplicaInd)
  Energies = zeros((NComp, NRep, NFrameRead), float)
  for i in ReplicaInd:
    dat = mdsim.GetConcatData(str(i), DataPath, Vars = Vars)
    for (j, v) in enumerate(Vars):
      Energies[j,i,:] = dat[v][NFrameSkip:NFrameSkip+NFrameRead]
  return Energies

def GetRestEnergies(DataPath, ReplicaInd, RestGroups, NFrameSkip, NFrameRead):
  "Returns computed components of the restraint energy function."
  #RestGroups is a list of RestList lists
  NComp, NRep = len(RestGroups), len(ReplicaInd)
  Energies = zeros((NComp, NRep, NFrameRead), float)
  if NComp == 0: return Energies
  #sort through replicas
  for i in ReplicaInd:
    Prefix = os.path.join(DataPath, "%d")
    TrjFn, PrmtopFn = GetTrjFiles(DataPath, i)
    Trj = coords.TrjClass(TrjFn, PrmtopFn, NSkip = NFrameSkip, NRead = NFrameRead)
    #cycle through coordinate
    j = 0
    Trj.Reset()
    for Pos in Trj:
      #cycle through energy components
      for (k, RestList) in enumerate(RestGroups):
        #skip if there are no restraints; else calculate energies
        if len(RestList) > 0:
          Energies[k,i,j] = mdsim.RestEnergy(Pos, RestList)
      j += 1
  return Energies


#======== PARSE RESTRAINTS IN REPLICAS ========

def GetRestWeights(AllRest):
  """Groups restraints into common parts with scale factors for replicas."""
  RestWeights = []
  N = len(AllRest)
  for (i, RestList) in enumerate(AllRest):
    for NewRest in RestList:
      #compare to existing restraints
      j = 0
      while j < len(RestWeights):
        Rest, Scale = RestWeights[j]
        #see if this is a factor of the other
        r = mdsim.RestCompare(Rest, NewRest)
        if r > 1:
          #replace main with this restraint
          RestWeights[j][0] = NewRest
          RestWeights[j][1] /= r
          RestWeights[j][1][i] = 1
          break
        elif r > 0:
          #add to scale for this group
          RestWeights[j][1][i] = r
          break
        else:
          j += 1
      #create new group if wasn't found
      if j >= len(RestWeights):
        Scale = zeros(N, float)
        Scale[i] = 1.
        RestWeights.append([NewRest, Scale])
  return RestWeights

def GetRestGroups(RestWeights):
  "Puts restraint scales into common groups."
  #now group these into separate groups
  RestGroups, Weights = [], []
  for (Rest, Scale1) in RestWeights:
    #see if it's already in there
    j = 0
    while j < len(RestGroups):
      RestList, Scale2 = RestGroups[j], Weights[j]
      if allclose(Scale1, Scale2):
        RestList.append(Rest)
        break
      else:
        j += 1
    #create new restgroup if not found
    if j >= len(RestGroups):
      RestGroups.append([Rest])
      Weights.append(Scale1)
  #convert to a 2-d array
  Weights = array(Weights, float)
  return RestGroups, Weights

def ParseRestTerms(DataPath, ReplicaInd):
  """Returns RestGroups and Weights from replica exchange data."""
  #first get restraint lists from data
  fn = os.path.join(DataPath, "restraints.dat.gz")
  if os.path.isfile(fn):
    AllRest = cPickle.load(gzip.GzipFile(fn, "r"))
  else:
    return [], []
  #filter for relevant indices
  AllRest = [x for (i,x) in enumerate(AllRest) if i in ReplicaInd]
  #put into groups
  RestWeights = GetRestWeights(AllRest)
  RestGroups, Weights = GetRestGroups(RestWeights)
  return RestGroups, Weights


#======== PERFORM WHAM CALCULATIONS ========

def RunWham(DataPath, NFrameSkip, NFrameRead,
  TargetTemp = TargetTempDflt, DeltaTemp = DeltaTempDflt,
  RemoveRest = True, Verbose = VerboseDflt):
  """Returns FK, LogwKN.
DataPath: where the replica exchange files are
NFrameSkip: number of frames to skip over
NFrameRead: number of frames to read
TargetTemp: target temperature for reweighting
DeltaTemp: maximum temperature range to use from replicas around TargetTemp
RemoveRest: True to remove the restraints in the weighting
Verbose: True to print messages"""
#HERE
  #check defaults
  TargetBeta = 1.0 / (TargetTemp * kB) 
  # get stats
  if Verbose: print "Getting REX stats"
  MinTemp, MaxTemp = TargetTemp - DeltaTemp, TargetTemp + DeltaTemp
  t1 = time.time()
  ReplicaInd, TempK, BetaK, NFrame = GetStats(DataPath, MinTemp, MaxTemp)
  NRep = len(ReplicaInd)
  t2 = time.time()
  if Verbose: print "Found %d frames, in %.0f sec." % (NFrame, t2 - t1)
  #check frame indicators
  NFrameSkip, NFrameRead = CheckFrames(NFrame, NFrameSkip, NFrameRead)
  if NFrameSkip == NFrame:
    if Verbose: print "Number of frames to skip equals total number of frames."
    return
  #get energies
  EIKN = GetEnergies(DataPath, ReplicaInd, ["EPOT", "EREST"],
                     NFrameSkip, NFrameRead)
  t3 = time.time()
  if Verbose: print "Got energy histories in %.0f sec" % (t3 - t2)
  #subtract off restraint energies from total potential energy
  EIKN[0] = EIKN[0] - EIKN[1]
  #get restraint groups
  RestGroups, WeightIK = ParseRestTerms(DataPath, ReplicaInd)
  if len(RestGroups) == 0:
    #just use current array, but make weights
    WeightIK = ones((2, NRep), float)
  elif len(RestGroups) > 0:
    #go through and recalculate energy terms
    ERest = GetRestEnergies(DataPath, ReplicaInd, RestGroups, NFrameSkip, NFrameRead)
    #concatenate to add the potential energy terms
    EIKN = concatenate((EIKN[:1], ERest))
    WeightIK = concatenate((ones((1, NRep), float), WeightIK))
    t4 = time.time()
    if Verbose: print "Got detailed restraint energy histories in %.0f sec" % (t4 - t3)
  #target weights
  NComp = EIKN.shape[0]
  if Verbose: print "%d terms found in energy function." % NComp
  TargetWeightI = ones(NComp, float)
  if RemoveRest: TargetWeightI[1:] = 0.
  #do wham
  if Verbose: print "Running WHAM"
  #run wham from outside module
  t5 = time.time()
  FK, LogwKN = LibWham(TempK, EIKN, TargetTemp, WeightIK, TargetWeightI, Verbose = Verbose)
  t6 = time.time()
  if Verbose: print "Ran WHAM in %.0f sec" % (t6 - t5)
  #return data
  return ReplicaInd, TempK, NFrame, NFrameSkip, NFrameRead, FK, LogwKN

def LibWham(TempK, EIKN, TargetTemp, WeightIK = None, TargetWeightI = None, kB = kB,
            MaxIter = 5000, RelTol = 1.e-8, Verbose = False):
  """Computes WHAM free energies using compiled Fortran lib. Returns FK, logwKN.
TempK: array of NReplica temperatures
EIKN: energy array for different components of the energy func;
  either dimen (N_replicas, N_frames) or (N_components, N_replicas, N_frames)
TargetTemp: target temperature
Weights: list of weight arrays for corresponding terms in energy function,
  each with dimension (NReplica). Note that it is assumed Energies gives the
  full energies, not the effective energies (multiplied by weights) in each
  of the replica Hamiltonians.
TargetWeights: list of target weights for the returned distribution
FK: dimensionless free energies (-F/kT) for each replica
logwKN: log weights of each configuration at target temp and weights
kN: different boltzmann constant, if desired
MaxIter: maximum number of wham iterations
AbsTol: stop wham iterations when difference between successive free energies
  is less than this value"""
  TargetBeta = 1.0 / (TargetTemp * kB)
  BetaK = 1.0 / (TempK * kB)
  #check defaults
  if len(EIKN.shape) == 2: EIKN = EIKN[newaxis,...]
  NComp, NRep, NFrame = EIKN.shape
  if WeightIK is None: WeightIK = ones((NComp, NRep), float)
  if TargetWeightI is None: TargetWeightI = ones(NComp, float)
  #do wham
  FK, NIter = whamlib.free_energy_multi(EIKN, WeightIK, BetaK, MaxIter, RelTol)
  LogwKN = whamlib.log_weight_multi(EIKN, WeightIK, BetaK, TargetWeightI, TargetBeta, FK)
  if Verbose:
    print "WHAM required %d iterations to converge to %8.3e relative tol." % (NIter, RelTol)
  return FK, LogwKN

def PyWham(TempK, EKN, ERestKN = None, TargetTemp = TargetTempDflt,
           NBin = 1000, NIter = 50, Verbose = VerboseDflt):
  "Computes WHAM free energies using a binning method, completely within Python."
  #This is a cheap, slow wham based on bins.  it's old code; don't use it!
  #get ranges
  EMin = EKN.min()
  EMax = EKN.max()
  EMax += (EMax - EMin) * 0.00001
  dE = (EMax - EMin) / float(NBin)
  invdE = 1./dE
  NRep, NFrame = EKN.shape
  BetaK = 1.0 / (TempK * kB)
  TargetBeta = 1.0 / (TargetTemp * kB)
  if ERestKN is None: ERestKN = zeros_like(EKN)
  #create histograms
  LogHist = zeros((NBin), float) + 1.e-8
  #bin energies
  ind = floor((EKN - EMin) * invdE).astype(int)
  for k in range(0,NRep):
    for n in range(0,NFrame):
      LogHist[ind[k,n]] += 1.
  #take logs
  LogHist = log(LogHist)
  LogNFrame = log(float(NFrame))
  HistIndVal = array([EMin + dE * (float(i)+0.5) for i in range(0,NBin)], float)
  #find average energies
  EAvgK = sum(EKN, 1) / float(NFrame)
  #create initial prob and free energy est
  fK = zeros((NRep), float)
  #initialize arrays
  LogP = zeros((NBin), float)
  LogwKN = zeros((NRep,NFrame), float)
  #initial iteration using bins
  for i in range(0, NIter):
    fKOld = fK.copy()
    for k in range(0, NRep):
      for j in range(0, NBin):
        LogNum = LogHist[j] - BetaK[k] * HistIndVal[j]
        LogDem = fK - BetaK * HistIndVal[j]
        LogDemMax = LogDem.max()
        LogP[j] = LogNum - log(sum(exp(LogDem-LogDemMax))) - LogDemMax - LogNFrame
      LogPMax = LogP.max()
      fK[k] = -log(sum(exp(LogP - LogPMax))) - LogPMax
    fK -= fK[0]
    if Verbose: print "Iteration %d:\n" % i
    if Verbose: print "\n".join(["%.4f -> %.4f" % (fKOld[k],fK[k]) for k in range(0,NRep)]) + "\n\n"
  #calculate individual probabilities
  for k in range(0, NRep):
    for n in range(0, NFrame):
      LogNum = -TargetBeta * (EKN[k,n] - ERestKN[k,n])
      LogDem = fK - BetaK * EKN[k,n]
      LogDemMax = LogDem.max()
      LogwKN[k,n] = LogNum - log(sum(exp(LogDem-LogDemMax))) - LogDemMax - LogNFrame
  #shift probabilities
  LogwKN = LogwKN - LogwKN.max()
  #return values
  return fK, LogwKN


#========= REWEIGHT PROBABILITIES ========

def ReweightProbs(DataPath, FileName, NBin,
  VarMin = None, VarMax = None, IndFunc = None, NVar = None, 
  NFrameSkip = 0, NFrameRead = None, ReplicaInd = None, LogwKN = None, 
  HeadRow = True, HeadCol = True, VarNames = None, 
  OutFile = None, AsLog = False, MultFactor = 1., Verbose = VerboseDflt):
  """Reweights a distribution of variable values taken from files at each temp.
  FileName will be added a prefix of "0.", "1.", ... for each temperature.
  ReplicaInd = None means that the file in FileName will be analyzed alone.
  IndFunc is an optional function of a string read from the variable files
  that returns the indices corresponding to the variables."""

  #check that all files are here
  if ReplicaInd is None:
    Files = [os.path.join(DataPath, FileName)]
  else:
    Files = [os.path.join(DataPath, "%d.%s" % (i, FileName)) for i in ReplicaInd]
  for fn in Files:
    if not os.path.isfile(fn):
      print "Cound not find %s" % fn
      return

  if FileName.strip().endswith(".gz"):
    FileMthd = gzip.GzipFile
  else:
    FileMthd = file

  #count the number of variables and get the names
  f = FileMthd(Files[0], "r")
  dat = f.readline()
  f.close()
  if NVar is None:
    NVar = len(dat.split())
    if HeadCol: NVar -= 1
  if VarNames is None:
    if HeadRow and IndFunc is None:
      VarNames = dat.split()[1:]
    else:
      VarNames = ["Var%d" % i for i in range(NVar)]

  if IndFunc is None:
    #initialize variable arrays
    VarDelta = (VarMax - VarMin) / float(NBin)
    ValOfInd = (arange(NBin) + 0.5)*VarDelta + VarMin
    ValFloat = True
    #make the indexing function using NBin, Varmin, Varmax
    def IndFunc(s):
      WordList = s.split()
      #get rid of head col
      if HeadCol: WordList = WordList[1:]
      #convert to float
      Vals = array([float(w) for w in WordList], float)
      #get bin number for values
      #put outliers in bin extrema
      Ind = clip( floor((Vals - VarMin) / VarDelta).astype(int), 0, NBin-1 )
      return Ind
  else:
    #just use consecutive numbering
    ValOfInd = array(range(0,NBin), int)
    ValFloat = False
    
  #initialize arrays    
  Prob = zeros((NVar, NBin), float)

  #init max free energy weight array
  if LogwKN is None:
    LogwKN = zeros((NVar,NBin), float)
  Minw = LogwKN.min()
  LogwMax = ones((NVar,NBin), float) * Minw
  
  #loop over files
  StageDesc = ["Finding maximum weights", "Binning observables"]
  for Stage in range(0,2):
    if Verbose: print StageDesc[Stage]
    
     #loop over temps
    for (k, fn) in enumerate(Files):
      #open data file
      f = FileMthd(fn, 'r')
      #skip over headers
      if HeadRow: f.readline()
      #skip over frames
      for n in range(NFrameSkip):
        f.readline()
      #start reading frames
      n = -1
      while True:
        #get values for variables for this frame
        s = f.readline()
        if s == "": break
        n += 1
        if not NFrameRead is None and n >= NFrameRead: break
        #find the indices
        Ind = IndFunc(s)
        if Stage == 0:
          #find the maximum
          for j in range(NVar):
            LogwMax[j, Ind[j]] = max(LogwMax[j, Ind[j]], LogwKN[k,n])
        else:
          #average the probabilities
          for j in range(NVar):
            Prob[j, Ind[j]] += exp( LogwKN[k,n] - LogwMax[j, Ind[j]] )
      f.close()
    
  if Verbose: print "Normalizing reweighted probabilities"
  #check if we want to put in log format rather than raw probabilities
  if AsLog:
    ClipMin = -PmfMax / abs(MultFactor)
    ProbSum = zeros(NVar, float)
    ProbMax = LogwMax.max(1)
    for j in range(NVar):
      ProbSum[j] = sum( Prob[j,:] * exp(LogwMax[j,:] - ProbMax[j]) ) 
    LogProbSum = log(ProbSum) 
    for j in range(NVar):
      for i in range(NBin):
       if Prob[j,i] == 0.0:
         Prob[j,i] = ClipMin
       else:
         Prob[j,i] = log(Prob[j,i]) + LogwMax[j,i] - LogProbSum[j] - ProbMax[j]
    Prob *= MultFactor
    Prob.clip(-PmfMax, PmfMax)
  else:
    #build back in the maximum values
    ProbMax = LogwMax.max(1)
    for j in range(NVar):
      for i in range(NBin):
        Prob[j,i] *= exp(LogwMax[j,i] - ProbMax[j])
    ProbSum = Prob.sum(1)
    for j in range(NVar):
      Prob[j,:] /= ProbSum[j]
    Prob *= MultFactor

  if not OutFile is None:
    WriteProbs(OutFile, VarNames, ValOfInd, Prob, ValFloat)  

  return ValOfInd, Prob


def WriteProbs(OutFile, VarNames, ValOfInd, Prob, ValFloat = True):
  "Writes reweighted probabilities to files."
  NVar, NBin = Prob.shape
  #check the type
  if ValFloat:
    ValFmt = "%-11.4f"
  else:
    ValFmt = "%-11d"
  #write file
  b = max([len(x) for x in VarNames])
  f = file(OutFile, 'w')
  f.write("var".ljust(b) + " " + " ".join([ValFmt % x for x in ValOfInd]) + '\n')
  for j in range(NVar):
    s = VarNames[j].ljust(b) + " "
    s += " ".join(["%-11.4f" % x for x in Prob[j,:]]) + "\n"
    f.write(s)
  f.close()


def WriteProbsTranspose(OutFile, VarNames, ValOfInd, Prob, ValFloat = True):
  "Writes reweighted probabilities to files, transposed to old-style format."
  NVar, NBin = Prob.shape
  #check the type
  if ValFloat:
    ValFmt = "%.4f"
  else:
    ValFmt = "%d"
  #write file
  b = max([13] + [len(x) for x in VarNames])
  f = file(OutFile, 'w')
  f.write("dist".ljust(b) + " " + " ".join([s.ljust(b) for s in VarNames]) + "\n")
  for i in range(0, NBin):
    f.write((ValFmt % ValOfInd[i]).ljust(b))
    for j in range(0, NVar):
      f.write(' ' + ("%.4f" % Prob[j,i]).ljust(b))
    f.write('\n')
  f.close()


def CalcDistProbs(DataPath, OutputPath, NFrameSkip, NFrameRead,
  PairList, DistMethod, NBin, MinR, MaxR, Cut,
  TargetTemp = TargetTempDflt, 
  TrjFile = None, PrmtopFile = None, 
  ReplicaInd = None, NFrame = None, LogwKN = None ,
  StartRes = 0, Prefix = "", Verbose = VerboseDflt, KeepDist = False):
  """Calculates PMFs and contact probabilities from either a single file or
  from REX/WHAM results, for a list of residue pairs."""

  #check for pairs
  if len(PairList) == 0:
    print "No contact pairs specified."
    return [], [], []

  #check for just one file
  if ReplicaInd is None:
    NFrame = coords.GetTrjLen(TrjFile, PrmtopFile)
  else:
    if not os.path.isdir(OutputPath):
      os.makedirs(OutputPath)

  #check frame bounds
  NFrameSkip, NFrameRead = CheckFrames(NFrame, NFrameSkip, NFrameRead)
  if NFrameSkip >= NFrame:
    print "Number of frames to skip equals total number of frames."
    return

  #output distance files
  if Verbose: print "Collecting distances"
  if ReplicaInd is None:
    mdtrj.SaveTrjResDists(TrjFile, PrmtopFile, OutputPath, PairList,
      NSkip = NFrameSkip, NRead = NFrameRead, DistMethod = DistMethod,
      Prefix = Prefix, StartRes = StartRes, Verbose = Verbose)
  else:
    mdtrj.SaveAllTrjResDists(DataPath, OutputPath, PairList,
      ReplicaInd = ReplicaInd, NSkip = NFrameSkip, NRead = NFrameRead,
      DistMethod = DistMethod, Prefix = Prefix, StartRes = StartRes, Verbose = Verbose)

  #calculate pmf
  if Verbose: print "Calculating PMFs"
  MultFactor = -(kB * TargetTemp)
  PmfDist, Pmf = ReweightProbs(OutputPath, Prefix + "dist.txt.gz", NBin, MinR, MaxR, 
    NFrameRead = NFrameRead, ReplicaInd = ReplicaInd, LogwKN = LogwKN,
    OutFile = os.path.join(OutputPath, Prefix + "pmf.txt"), AsLog = True,
    MultFactor = MultFactor, Verbose = True)

  #calculate contact probs  
  if Verbose: print "Calculating contact probabilities"
  CutDist, CutProb = ReweightProbs(OutputPath, Prefix + "dist.txt.gz", 2, 0., 2.*Cut, 
    NFrameRead = NFrameRead, ReplicaInd = ReplicaInd, LogwKN = LogwKN,
    OutFile = os.path.join(OutputPath, Prefix + "contactprob.txt"), AsLog = False, Verbose = True)
  #reformat CutProb
  CutProb = [x[0] for x in CutProb]
    
  #write stats file
  f = file(os.path.join(OutputPath, Prefix + "analstats.txt"), "w")
  f.write("Number of frames skipped: %d\n" % NFrameSkip)
  f.write("Number of frames read   : %d\n" % NFrameRead)
  f.write("Distance method         : %s\n" % mdtrj.DistMethodDesc[DistMethod])
  f.write("PMF number of bins      : %d\n" % NBin)
  f.write("PMF bin minimum         : %.2f\n" % MinR)
  f.write("PMF bin maximum         : %.2f\n" % MaxR)
  f.write("Cutoff distance         : %.2f\n" % Cut)
  f.write("Target temperature      : %.2f\n" % TargetTemp)
  f.write("Number of temperatures  : %d\n" % len(ReplicaInd))
  f.write("Replica indices         : %s\n" % repr(list(ReplicaInd)))
  f.close()

  if not KeepDist:
    mdtrj.DeleteAllTrjDists(OutputPath)

  #return probs
  return PmfDist, Pmf, CutProb


def CalcAtomDistProbs(DataPath, OutputPath, NFrameSkip, NFrameRead,
  PairAtoms, NBin, MinR, MaxR, 
  TargetTemp = TargetTempDflt, 
  TrjFile = None, PrmtopFile = None, 
  ReplicaInd = None, NFrame = None, LogwKN = None ,
  PairLabels = None, Prefix = "", Verbose = VerboseDflt, KeepDist = False):
  """Calculates PMFs and contact probabilities from either a single file or
  from REX/WHAM results, for a list of residue pairs."""

  #check for pairs
  if len(PairAtoms) == 0:
    print "No contact pairs specified."
    return [], [], []

  #check for just one file
  if ReplicaInd is None:
    NFrame = coords.GetTrjLen(TrjFile, PrmtopFile)
  else:
    if not os.path.isdir(OutputPath):
      os.makedirs(OutputPath)

  #check frame bounds
  NFrameSkip, NFrameRead = CheckFrames(NFrame, NFrameSkip, NFrameRead)
  if NFrameSkip >= NFrame:
    print "Number of frames to skip equals total number of frames."
    return [], [], []

  #output distance files
  if Verbose: print "Collecting distances"
  if ReplicaInd is None:
    mdtrj.SaveTrjDists(TrjFile, PrmtopFile, OutputPath, PairAtoms,
      PairLabels = PairLabels, NSkip = NFrameSkip, NRead = NFrameRead,
      Prefix = Prefix, Verbose = Verbose)
  else:
    mdtrj.SaveAllTrjDists(DataPath, OutputPath, PairAtoms,
      PairLabels = PairLabels, ReplicaInd = ReplicaInd, NSkip = NFrameSkip,
      NRead = NFrameRead, Prefix = Prefix, Verbose = Verbose)

  #calculate pmf
  if Verbose: print "Calculating PMFs"
  MultFactor = -(kB * TargetTemp)
  PmfDist, Pmf = ReweightProbs(OutputPath, Prefix + "dist.txt.gz", NBin, MinR, MaxR, 
    NFrameRead = NFrameRead, ReplicaInd = ReplicaInd, LogwKN = LogwKN,
    OutFile = os.path.join(OutputPath, Prefix + "pmf.txt"), AsLog = True,
    MultFactor = MultFactor, Verbose = True)
 
  #write stats file
  f = file(os.path.join(OutputPath, Prefix + "analstats.txt"), "w")
  f.write("Number of frames skipped: %d\n" % NFrameSkip)
  f.write("Number of frames read   : %d\n" % NFrameRead)
  f.write("PMF number of bins      : %d\n" % NBin)
  f.write("PMF bin minimum         : %.2f\n" % MinR)
  f.write("PMF bin maximum         : %.2f\n" % MaxR)
  f.write("Target temperature      : %.2f\n" % TargetTemp)
  f.write("Number of temperatures  : %d\n" % len(ReplicaInd))
  f.write("Replica indices         : %s\n" % repr(list(ReplicaInd)))
  f.close()

  if not KeepDist:
    mdtrj.DeleteAllTrjDists(OutputPath)

  #return probs
  return PmfDist, Pmf


def TestWham():
  "Runs tests for wham routines."
  #MAKE VARIABLES
  NRep = 10
  TempK = arange(270., 269. + 10.*NRep, 10.).astype(float)
  BetaK = 1./(TempK * kB)
  NFrame = 1000
  EIKN = 100.*random.rand(3*NRep*NFrame).astype(float).reshape((3, NRep, NFrame))
  for i in range(NRep):
    EIKN[0,i,:] += i*5.
    EIKN[1,i,:] += i*10.
    EIKN[2,i,:] += i*20.
  WeightIK = random.rand(3*NRep).astype(float).reshape((3, NRep))
  TargetWeightI = random.rand(3).astype(float)
  TargetTemp = TargetTempDflt
  TargetBeta = 1.0 / (TargetTempDflt * kB)
  #FIRST SIMPLE; NO RESTRAINTS
  print "==NO RESTRAINTS=="
  print "Running whamlib wham."
  FK1 = whamlib.free_energy(EIKN[0,:,:], BetaK, 100, 50, 5)
  LogwKN1 = whamlib.log_weight(EIKN[0,:,:], BetaK, TargetBeta, FK1)
  print "Running wrapped whamlib wham."
  FK2, LogwKN2 = LibWham(TempK, EIKN[0,:,:], TargetTemp, Verbose = True)
  print "Running Python wham."
  FK3, LogwKN3 = PyWham(TempK, EIKN[0,:,:], NBin = 100, NIter = 50, Verbose=False)
  LogwKN1 -= LogwKN1.max()
  LogwKN2 -= LogwKN2.max()
  LogwKN3 -= LogwKN3.max()
  print "Free energies1: ", ", ".join(["%.4f" % x for x in FK1])
  print "Free energies2: ", ", ".join(["%.4f" % x for x in FK2])
  print "Free energies3: ", ", ".join(["%.4f" % x for x in FK3])
  print "Average differences of configuration weights: ", average(abs(LogwKN1 - LogwKN2)), \
        average(abs(LogwKN1 - LogwKN3))
  #SECOND TEST; WITH 1 RESTRAINT
  print "==ONE RESTRAINT=="
  print "Running wrapped whamlib wham."
  FK1, LogwKN1 = LibWham(TempK, EIKN[:2], TargetTemp, WeightIK = WeightIK[:2],
                         TargetWeightI = TargetWeightI[:2], Verbose = True)
  print "Running two-component whamlib wham."
  FK2 = whamlib.free_energy_2(EIKN[0,:,:], EIKN[1,:,:], WeightIK[0,:], WeightIK[1,:],
                              BetaK, 100, 1000, 50)
  LogwKN2 = whamlib.log_weight_2(EIKN[0,:,:], EIKN[1,:,:], WeightIK[0,:], WeightIK[1,:],
                                 BetaK, TargetWeightI[0], TargetWeightI[1], TargetBeta, FK2)
  LogwKN1 -= LogwKN1.max()
  LogwKN2 -= LogwKN2.max()
  print "Free energies1: ", ", ".join(["%.4f" % x for x in FK1])
  print "Free energies2: ", ", ".join(["%.4f" % x for x in FK2])
  print "Average differences of configuration weights: ", average(abs(LogwKN1 - LogwKN2))
  #THIRD TEST; WITH 2 RESTRAINTS
  print "==TWO RESTRAINTS=="
  print "Running wrapped whamlib wham."
  FK1, LogwKN1 = LibWham(TempK, EIKN, TargetTemp, WeightIK = WeightIK,
                         TargetWeightI = TargetWeightI, Verbose = True)
  LogwKN1 -= LogwKN1.max()
  print "Free energies1: ", ", ".join(["%.4f" % x for x in FK1])


  

