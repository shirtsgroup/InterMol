#!/usr/bin/env python

#last modified: 03-08-07

import os, math, cPickle, time, sys, glob, copy, shutil
import pdbtools, sequence, rexsock, makeconf, protein, assemble
import zamdata, zamrex

#GLOBALS
MaxInitConf = 50        #maximum number of initial conformations
MaxInitAssem = 10       #maximum number of pairwise assemblies to attempt
RexTimeBase = 5000.     #baseline time in ps to REX
RexTimePerAdd = 0.      #addl time in ps to REX per added residue
RexTimePerDOF = 0.      #addl time in ps to REX per degree of freedom
RexTimeMin = 4000.      #minimum total time to REX for any fragment
RexTimeSwapRest = 2000. #time for swapping restraints
AnalTime = 1000.        #time in ps to use for analysis

#timeouts
InitConfTimeout = None #timeout in seconds for the initial configuration making, or None


def PrepConfs(Frag):
  "Returns all files to be used for the initial configuration generation."
  #make the folder
  PdbPath = os.path.join(Frag.BasePath, "conf/")
  if not os.path.isdir(PdbPath):
    os.mkdir(PdbPath)
  #add any user-specified initial configurations to the top of the list
  PdbFiles = glob.glob(os.path.join(PdbPath, "seed*.pdb"))
  PdbFiles.sort()
  n = len(PdbFiles)
  Weights = [2. - i/float(n) for i in range(n)]
  #add parent fragments
  for (PFrag, Ind) in Frag.Parents:
    PSnapPdbs = [os.path.join(PFrag.BasePath, "conf/", f)
                 for f in PFrag.SnapPdbs]
    #get the weights of the parent snapshots;
    #need to accomodate old data before SnapWeight was added
    if len(PFrag.SnapWeight) == 0:
      m = len(PSnapPdbs)
      x = 2. / (m + 1.)
      PWeights = [x - x*i/float(m) for i in range(m)]
    else:
      PWeights = PFrag.SnapWeight
    #add the relevant snapshots
    if Ind < 0:
      PdbFiles.extend(PSnapPdbs)
      Weights.extend(PWeights)
    else:
      PdbFiles.append(PSnapPdbs[Ind])
      Weights.append(1.)
  return PdbFiles, Weights
    

def MakeInitConfs(Frag, PdbFiles, Weights):
  "Makes initial configurations for replica exchange."
  #make alignments
  LogFile = os.path.join(Frag.BasePath, "conf/seedlog.txt")
  Prefix = os.path.join(Frag.BasePath, "conf/init")
  if len(PdbFiles) > 0:
    #make initial conformations from old 
    if Frag.Stage == "grow":
      makeconf.MakePdbAlignments(Frag.Seq, PdbFiles, Prefix,
        MaxAlign = MaxInitConf, PdbWeights = Weights,
        CapN = Frag.CapN, CapC = Frag.CapC, Dehydrogen = True,
        LogFile = LogFile)
    elif Frag.Stage == "assem":
      makeconf.MakePdbAlignments(Frag.Seq, PdbFiles, Prefix,
        MaxAlign = MaxInitAssem, MaxConf = MaxInitConf, PdbWeights = Weights,
        CapN = Frag.CapN, CapC = Frag.CapC, Dehydrogen = True,
        LogFile = LogFile, MC = True)
  elif len(Frag.SS) > 0:
    #make conformations from secondary structure
    assemble.OptimizeMCFromSS(Frag.Seq, Frag.SS, Prefix, MaxInitConf,
                              Verbose = False, Summary = True)
    

def RunInitConfTask(sr, Frag, Task):
  """Runs initial configuration making as a socketring task."""
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
    #first prep the configuration making
    PdbFiles, Weights = PrepConfs(Frag)
    #make the task
    s = "import zamrun, zamdata"
    s += "\nFragData, PdbFiles, Weights = SockData"
    s += "\nFrag = zamdata.FragDataClass()"
    s += "\nFrag.LoadsData(FragData)"
    s += "\nzamrun.MakeInitConfs(Frag, PdbFiles, Weights)"
    s += "\nSockResult = Frag.DumpsData()"
    SockData = (Frag.DumpsData(), PdbFiles, Weights)
    Task = sr.AddTask(s, Data = SockData, Timeout = InitConfTimeout)
    return False, Task


def SetupFrag(rx, zd, Frag):
  """Sets up the replica exchange simulation based on Frag.RexTime."""
  #check to see if we should calculate the rex time
  if Frag.CalcRexTime:
    #calculate the simulation time
    Frag.RexTime = max(RexTimeMin, Frag.ResDOF() * RexTimePerDOF +
      Frag.NResAdd() * RexTimePerAdd + RexTimeBase)
  #see if we should calculate the initial restraints
  if Frag.CalcSwapRest:
    #make a list of the parent contacts, with best contacts first
    RestList = []
    for (f, ind) in Frag.Parents:
      RestList.extend(f.PropRest)
    RestList.sort(reverse=True)
    RestList = [(a,b) for (Score,a,b) in RestList
                if Frag.RestInRange(a,b)]
    Frag.SwapRest = RestList
  #see if we should calculate the initial rest time
  if Frag.CalcSwapRestTime and len(Frag.SwapRest) > 0:
    Frag.RexTimeSwapRest = RexTimeSwapRest  
  #calculate the time elapsed so far
  ElapsedTime = rx.NCycle * Frag.CycleTime
  #set up the appropriate REX
  if ElapsedTime < Frag.RexTimeSwapRest:
    Frag.RexStage = 0
    zamrex.SetupRex(rx, zd, Frag, RestListSort = Frag.SwapRest)
  else:
    Frag.RexStage = 1
    zamrex.SetupRex(rx, zd, Frag)
  #set the analysis time
  if Frag.CalcAnalTime: Frag.AnalTime = AnalTime


def RunFrag(rx, zd, Frag):
  """Runs the REX protocol.  Returns True if the REX
  simulations have finished, False otherwise."""
  #see if we need to initialize a new replica exchange
  if not rx.RexRunning:
    SetupFrag(rx, zd, Frag)
  #run some
  if zamrex.RunRex(rx, zd, Frag):
    #ok we finished this round.
    #if we are the last stage then we are done; otherwise setup next
    if Frag.RexStage == 1:
      return True
    else:
      #setup the next stage
      SetupFrag(rx, zd, Frag)
  #we're not yet done
  return False

