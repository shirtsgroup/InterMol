#!/usr/bin/env python

#LAST MODIFIED: 05-06-07

#TODO:

from numpy import *
import sys, os, gzip, time, random, math, copy, cPickle, glob, tempfile
import mdsim
import socketring

#use temporary (scratch) directory for running sims?
UseTempDir = True
#name of scratch directory to look for
TempDirBase = "/scratch"

#by default, uses gzip to save big files
fileMthd = gzip.GzipFile
Gzip = True
fileExt = ".gz"

#default path names
WorkPathName = "work"
DataPathName = "data"

#replica state constants (these must be negative
#since positive values correspond to proc numbers)
StateDone = -2
StateWaiting = -1

#processor constants (positive values correspond
#to state numbers)
ProcReady = -1
ProcMain = -2

#data to save
RestartData = ["MoveAtt", "MoveAcc", "NCycle", "UserData",
               "TempOptProb", "TempOptCount", "TempOptFlag",
               "Minimized", "Replica", "State", "Swaps",
               "LastSwapProbs"]

#wait interval in sec
WaitInterval = 0.1

#time limit in sec for updating auxiliary files (current.pdb, prmtop.parm7)
FullUpInterval = 10.*60.

#turn on temperature optimization
TempOptDflt = False
TempOptMix = 0.5

#debug flag
DBG = False

#timeout for any one md cycle in sec (or None)
CycleTimeoutDflt = 3600.



#======== COMMON ACCEPTANCE CRITERIA ========

def PaccFuncBasic(a, b):
  """Returns acceptance probability based on total EPOT alone,
and doesn't account for possible differences in Hamiltonians."""
  dBeta = 1./(mdsim.kB * a["TEMPSET"]) - 1./(mdsim.kB * b["TEMPSET"])
  dE = a["EPOT2"] - b["EPOT2"]
  return math.exp(min(0., dE*dBeta))

def PaccFuncRest(a, b):
  """Returns acceptance probability, accounting for differences
in restraint energies but assuming rest of Hamiltonian is the same."""
  if mdsim.SameRest(a,b):
    dBeta = 1./(mdsim.kB * a["TEMPSET"]) - 1./(mdsim.kB * b["TEMPSET"])
    dE = a["EPOT2"] - b["EPOT2"]
    return math.exp(min(0., dE*dBeta))
  else:
    #get the positions
    Posa, Posb = a.Pos, b.Pos
    if Posa is None: Posa = a.GetPos()
    if Posb is None: Posb = b.GetPos()
    #calculate energies before and after swaps
    Eaina = mdsim.RestEnergy(Posa, a.RestList)
    Eainb = mdsim.RestEnergy(Posa, b.RestList)
    Ebina = mdsim.RestEnergy(Posb, a.RestList)
    Ebinb = mdsim.RestEnergy(Posb, b.RestList)
    Betaa, Betab = 1./(mdsim.kB * a["TEMPSET"]), 1./(mdsim.kB * b["TEMPSET"])
    Term1 = ((a["EPOT2"] - a["EREST2"]) - (b["EPOT2"] - b["EREST2"])) * (Betaa - Betab)
    Term2 = Betaa * Eaina + Betab * Ebinb - Betab * Eainb - Betaa * Ebina
    return math.exp(min(0., Term1 + Term2))    


#======== THE REPLICA EXCHANGE CLASS ========

class RexClass:

  def __init__(self, SockRing = None, BasePath = None, N = None,
    DataList={"TEMPAVG":"tempavg.txt"}, Swaps = None,
    Restart = True):
    """Initializes a new replica exchange class and selects a main
server.  Optionally if BasePath and N are specified, will make
directories and initialize all variables so that a subsequent
call to Initialize() is not required.
* SockRing: a RingClass instance from socketring; automatically
  created if None is specified
* BasePath: string containing the directory to store class files
* N: integer number of states (e.g., temperatures)
* DataList: a dictionary whose keys are names of variables in
  the individual simulation instances (e.g., "PRES" as in
  Sim["PRES"]), and the values are the names of files to save
  the variables at every exchange point (e.g., "pressure.txt").
  Default is average temperature saved to "tempavg.txt".
* Swaps: a list of (a,b) pairs of allowed swaps between states.
  Default is neighboring state pairs.
* Restart: Boolean value.  True will look for a restart file
  and load old data if present.  If restart data is not there
  or if False the simulation starts anew.  Default is True."""
   
    #DETERMINE SERVER/CLIENT ROLES
    if SockRing is None:
      #default path is current
      self.sr = socketring.RingClass()
    else:
      self.sr = SockRing
    if self.sr.Finished():
      print """It appears the REX simulation has finished.
To add additional simulation cycles, first delete the file server.txt"""
      return

    if self.sr.IsServer:
      #WE ARE THE MAIN SERVER
      print "Starting REX simulation as main server"    

      #set the user dictionary
      self.UserData = {}

      #initialize the rest if desired
      if BasePath is not None and N is not None:
        self.InitServer(BasePath, N, DataList, Swaps, Restart)

    else:
      #we are not the main server; setup as client
      print "Server is " + self.sr.Server

    self.Initialized = False

    #flush std out
    sys.stdout.flush()


  #BASE SERVER ROUTINES

  def InitServer(self, BasePath, N, DataList={"TEMPAVG":"tempavg.txt"},
    Swaps = None, Restart = True):
    """Initializes variables and makes directories.  Run on the
main server only.
* BasePath: string containing the directory to store class files
* N: integer number of states (e.g., temperatures)
* DataList: a dictionary whose keys are names of variables in
  the individual simulation instances (e.g., "PRES" as in
  Sim["PRES"]), and the values are the names of files to save
  the variables at every exchange point (e.g., "pressure.txt").
  Default is average temperature saved to "tempavg.txt".
* Swaps: a list of (a,b) pairs of allowed swaps
  between states.  Default is neighboring state pairs.
* Restart: Boolean value.  True will look for a restart file
  and load old data if present.  If restart data is not there
  or if False the simulation starts anew.  Default is True."""
    
    if not self.sr.IsServer: return
    
    #SET PATHS
    #BasePath is the base path for REX simulation
    #must use absolute path here or temp dir won't work on clients
    self.BasePath = os.path.abspath(BasePath)
    if not os.path.isdir(self.BasePath): os.mkdir(self.BasePath)
    #DataPath is where all datafiles will go
    self.DataPath = os.path.join(self.BasePath, DataPathName)
    if not os.path.isdir(self.DataPath): os.mkdir(self.DataPath)
    #WorkPath is where the replica simulations will be run
    self.WorkPath = os.path.join(self.BasePath, WorkPathName)
    if not os.path.isdir(self.WorkPath): os.mkdir(self.WorkPath)

    #INITIALIZE REPLICAS
    self.NReplica = N
    #State[] gives the state number as a function of replica
    self.State = [i for i in range(0,self.NReplica)]
    #Replica[] gives the replica number as a fn of state number (e.g., temp)
    self.Replica = [i for i in range(0,self.NReplica)]
    #Tasks[] gives the task as a function of state number
    self.Tasks = [None]*self.NReplica
    #Swaps[] is an array of tuples (i,j) giving allowed swaps between states i & j
    if Swaps is None:
      self.Swaps = [(i,i+1) for i in range(0,self.NReplica-1)]
    else:
      self.Swaps = [(min([a,b]),max([a,b])) for (a,b) in Swaps]
    #False to allow replicas to be swapped only once each
    self.MultiSwap = True
    #Mode for swaps to occur.  0 means each swap is performed on disk before
    #proceeding to the next swap.  1 means all positions are loaded in memory
    #first, swapped in memory, and then saved to disk afterwards.  0 is more
    #disk intensive, 1 is more memory intensive.
    self.SwapMode = 1
    
    #initialize the system classes for each replica
    #Sim[] is indexed by state number, not replica
    self.Sim = [mdsim.SimClass(os.path.join(self.WorkPath, str(i)))
           for i in range(0, self.NReplica)]
    #go ahead and make directories for simulations
    for i in range(0,self.NReplica):
      if not os.path.isdir(self.Sim[i].RunPath): os.mkdir(self.Sim[i].RunPath)
    #add a user variable with the state number and a user variable for task time
    for (i,s) in enumerate(self.Sim):
      s["REXID"] = i
      s["TASKTIME"] = 0.
      
    #INITIALIZE COUNTS
    #number of exchange cycles
    self.NCycle = 0
    #number of cycles at the start of this run; may be
    #greater than zero if a restart was performed
    self.NCycleStart = 0
    #counts for the exchange move attempts and accepts
    self.MoveAtt = [0.]*len(self.Swaps)
    self.MoveAcc = [0.]*len(self.Swaps)
    #last swap probabilities
    self.LastSwapProbs = [0.]*len(self.Swaps)
    
    #time the class was initialized
    self.TimeStart = time.time()
    #time of penultimate update
    self.TimeLastUp = time.time()
    #time of most recent update
    self.TimeThisUp = time.time()
    #time of the last full update (current.prb, prmtop.parm7)
    self.TimeFullUp = 0.
    #do a full update this round?
    self.FullUp = True
    #number of cycles at last full update
    self.NCycleFullUp = 0
    #whether or not initial minimization has been performed
    self.Minimized = False

    #storage sizes
    self.DataStorageSize = 0
    self.WorkStorageSize = 0

    #TEMP OPTIMIZATION
    #use temperature optimization?
    self.TempOpt = TempOptDflt
    #cycles at which to perform update
    self.TempOptCycles = []
    self.TempOptReset()
    
    #load the restart file if it's there
    if Restart and os.path.isfile(os.path.join(self.WorkPath,"restart.dat")):
      self.RestartRead()
      self.Restart = True
      print "Restarting REX from last simulation data"
    else:
      self.Restart = False

    #this is the alternator variable for alternating
    #swap cascades up and down
    self.LastCascade = 1
    
    #SETUP DATA FILES
    #initialize the list of variable monitors
    self.SaveDataList = DataList
    #reset the data path
    if not self.Restart: self.ResetData()

    #we're initialized now
    self.Initialized = True

    #REX run variables   
    #is a REX in progress?
    self.RexRunning = False
    #reference to the acceptance function
    self.PaccFunc = None
    #total num of cycles requested
    self.NCycleTot = 0
    #progress file name
    self.ProgFile = None
    #perform swaps during REX?
    self.PerformSwaps = True
    #temperature (state) labels
    self.Labels = []
    #swap restraints as well as config?
    self.SwapRest = False

    #TIME USAGE 
    #update time
    self.SecUpdate = 0.
    #swapping time
    self.SecSwap = 0.
    #restart write time
    self.SecWrite = 0.
    #estimate of seconds remaining
    self.SecRemainEst = 0.
    #client times
    self.SecClient = [0. for i in range(self.NReplica)]

    #TIMEOUT
    self.CycleTimeout = CycleTimeoutDflt


  def __contains__(self, key):
    """Allows testing for a user variable."""
    return key in self.UserData

  def __getitem__(self, key):
    """Allows user to get variables using dictionary notation."""
    if self.UserData.has_key(key):
      return self.UserData[key]
    else:
      raise KeyError, "Key not found."

  def __setitem__(self, key, val):
    """Allows user to set variables using dictionary notation.
If user key does not exist, it will be added to the
list of user variables."""
    self.UserData[key] = val

  def __delitem__(self, key):
    """Allows user to use del to remove a user variable."""
    if self.UserData.has_key(key):
      del self.UserData[key]
    else:
      print "Key not found."    


  #RESTART FUNCTIONS
    
  def RestartWrite(self, SaveSimData = True):
    "Writes a restart file."
    self.SecWrite = -time.time()
    OldRstFile = os.path.join(self.WorkPath, "restart.dat")
    NewRstFile = os.path.join(self.WorkPath, "restartnew.dat")
    f = open(NewRstFile, "w")
    #dump a pickled string of restart data
    cPickle.dump([(itm,self.__dict__[itm]) for itm in RestartData], f,
                 cPickle.HIGHEST_PROTOCOL)
    f.close()
    #write over the current file
    if os.path.isfile(OldRstFile): os.remove(OldRstFile)
    os.rename(NewRstFile, OldRstFile)
    self.SecWrite += time.time()
    if SaveSimData:
      for s in self.Sim: s.SaveData()

  def RestartRead(self):
    "Reads all restart data."
    fn = os.path.join(self.WorkPath, "restart.dat")
    if os.path.isfile(fn):
      f = open(fn, "r")
      #load a pickled string of restart data
      alldata = cPickle.load(f)
      #parse items in the string
      for (itm, data) in alldata:
          self.__dict__[itm] = data       
      f.close()
    #update the NCycleStart variable with the restart data
    self.NCycleStart = self.NCycle
    #legacy in case didn't save LastSwapProbs:
    if not len(self.LastSwapProbs) == len(self.Swaps):
      self.LastSwapProbs = [0.] * len(self.Swaps)
    #load old simulation data for each state
    for s in self.Sim:
      s.LoadData()

  def RestraintsWrite(self):
    "Writes restraint data to a file."
    fn = os.path.join(self.DataPath, "restraints.dat.gz")
    cPickle.dump([s.RestList for s in self.Sim], gzip.GzipFile(fn, "w"))

  def RestraintsRead(self):
    "Reads restraint data from a file, for analysis purposes."
    fn = os.path.join(self.DataPath, "restraints.dat.gz")
    return cPickle.load(gzip.GzipFile(fn, "r"))    


  #ITERATORS OVER STATES
    
  def SimIter(self, func):
    """Iterates a function over all states.  Returns an
array of the results of the function applied to each state.
* func: a function to be iterated over the SimClass objects
  as func(s) where s is a SimClass instance"""
    return [func(self.Sim[i]) for i in range(0,self.NReplica)]

  def SimIterWithVar(self, Var, func):
    """Iterates a function over all states with a variable array.
Returns an array of the function results applied to each state.
* Var: an array of length NReplica containing variable values
* func: a function to be iterated over the SimClass objects,
  as func(s,Var[i)] where s is a SimClass instance and i is the
  respective index of s in the REX simulation."""
    return [func(self.Sim[i], Var[i]) for i in range(0,self.NReplica)]

  def SimSetVar(self, VarName, VarVal):
    """Sets dictionary variables in all states using a particular
variable name and an array of values.
* VarName: string name of the variable for assignment,
  which will be applied as s[VarName] = VarVal[i] where
  s is a SimClass object and i is its respective index
* VarVal: array of length NReplica containing the variable
  values for each state"""
    for i in range(0,self.NReplica):
      self.Sim[i][VarName] = VarVal[i]

  def UnifySettings(self):
    "Copies settings from state 0 to all other replicas."
    for i in range(1,self.NReplica):
      self.Sim[i].CopyData(self.Sim[0])


  #TASK MANAGEMENT
      
  def TasksSend(self, ExecStr, FullTempCopy = True):
    """Sends out a task to be run on all states.
* ExecStr: string representation of a function to be called
  by each state.  In this string, the token SIM will be
  replaced by a SimClass object.
* FullTempCopy: True will copy all files produced in a temporary
  scratch dir, if that option is enabled.  False will only copy
  current configuration (used during REX cycles)."""
    #create and add the tasks
    for i in range(0, self.NReplica):
      SimPath = os.path.join(self.BasePath, "work/%d" % i)
      s = "import mdsim, time\n"
      s += "SIM = mdsim.SimClass('%s')\n" % SimPath
      s += "SIM.LoadsData(SockData)\n"
      if UseTempDir:
        #copy data to scratch directory
        s += "if not 'RexSockTempDir' in globals():\n"
        s += "  import tempfile, os\n"
        s += "  if os.path.isdir('%s'):\n" % TempDirBase
        s += "    globals()['RexSockTempDir'] = tempfile.mkdtemp(dir = '%s')\n" % TempDirBase
        s += "  else:\n"
        s += "    globals()['RexSockTempDir'] = tempfile.mkdtemp()\n"
        s += "  SockTempFiles.append(globals()['RexSockTempDir'])\n"
        s += "  print 'Temporary REX folder is:', globals()['RexSockTempDir']\n"
        s += "SIM.MovePath(globals()['RexSockTempDir'], PrepOnly = True)\n"
      s += "SIM['TASKTIME'] = -time.time()\n"
      s += ExecStr + "\n"
      s += "SIM['TASKTIME'] = SIM['TASKTIME'] + time.time()\n"
      if UseTempDir:
        #copy data back to master work directory
        if FullTempCopy:
          s += "SIM.MovePath('%s')\n" % SimPath
        else:
          s += "SIM.MovePath('%s', CurrentOnly = True)\n" % SimPath
      s += "SIM.SaveData()\n"
      s += "SockResult = SIM.DumpsData()\n"
      Data = self.Sim[i].DumpsData()
      self.Tasks[i] = self.sr.AddTask(s, Data, Timeout = self.CycleTimeout)

  def TasksCheck(self, Timeout = 0.):
    """Checks to see if all recent tasks have finished.
If so, returns True.  If not, returns False.
* Timeout: float of number of seconds to wait during check."""
    StartTime = time.time()
    #check to see if no tasks have been assigned
    while time.time() - StartTime < Timeout or Timeout == 0.:
      #see if all the tasks are done
      TaskList = [t for t in self.Tasks if not t is None]
      if [t.Done for t in TaskList].count(False) == 0:
        #all of the tasks are done so clean them out
        for Task in TaskList:
          #update the simulation data
          Ind = self.Tasks.index(Task)
          if len(Task.Result) == 0:
            print "Error: did not receive simulation result."
          self.Sim[Ind].LoadsData(Task.Result)
        #get the elapsed times
        for i in range(0, self.NReplica):
          t = self.Tasks[i]
          if t is None:
            self.SecClient[i] = None
          else:
            self.SecClient[i] = t.ElapsedTime()
        return True
      else:
        time.sleep(WaitInterval)
      if Timeout == 0: break
    return False

  def TasksPurge(self):
    """Removes all done tasks from the task list."""
    for i in range(0, self.NReplica):
      if not self.Tasks[i] is None:
        if self.Tasks[i].Done:
          self.sr.RemoveTask(self.Tasks[i])
          self.Tasks[i] = None
          
  def TasksWait(self):
    """Waits for all tasks to complete."""
    while not self.TasksCheck():
      #run tasks
      self.sr.RunTasks(1.0)

  def TasksExist(self):
    """Returns True if any tasks exist, running or done."""
    return self.Tasks.count(None) < self.NReplica

  def SimFarm(self, ExecStr, FullTempCopy = True):
    """Farms out a task to all states and waits for completion.
* ExecStr: string representation of a function to be called
  by each client.  In this string, the token SIM will be
  replaced by a SimClass object.
* FullTempCopy: True will copy all files produced in a temporary
  scratch dir, if that option is enabled.  False will only copy
  current configuration (used during REX cycles)."""
    #make the tasks
    self.TasksSend(ExecStr, FullTempCopy)
    #wait for the tasks
    self.TasksWait()
    #purge the tasks
    self.TasksPurge()


  #DATA COLLECTION FUNCTIONS

  def UpdateStorageSize(self):
    "Updates the storage space in bytes for the data and work files."
    if os.path.exists(self.DataPath):
      self.DataStorageSize = sum([os.path.getsize(os.path.join(self.DataPath,f))
                  for f in os.listdir(self.DataPath)])
    else:
      self.DataStorageSize = 0
    self.WorkStorageSize = sum([s.StorageSize() for s in self.Sim])


  def Update(self):
    "Updates all REX data after one cycle."

    #increment the cycle count and update the times
    self.NCycle += 1
    self.TimeLastUp = self.TimeThisUp
    self.TimeThisUp = time.time()
    
    #update variable monitors
    cwd = os.getcwd()
    os.chdir(self.DataPath)
    line = " ".join([str(s)[0:3] for s in self.State]) + "\n"
    fileMthd("statenum.txt" + fileExt, "a").write(line)
    line = " ".join([str(s)[0:3] for s in self.Replica]) + "\n"
    fileMthd("replicanum.txt" + fileExt, "a").write(line)

    if self.TempOpt:
      #update the histograms
      self.TempOptCycle()
      #update the history file
      #see if we need to update the temperature
      if self.NCycle in self.TempOptCycles:
        Stat = self.TempOptUpdate()
        #write a new temperature file and update the history file
        Temps = ["%.3f" % s["TEMPSET"] for s in self.Sim]
        file("temps.txt", "w").write("\n".join(Temps) + "\n")
        s = "NCycle = %d\n%s\n" % (self.NCycle, Stat)
        if min(self.TempOptCycles) == self.NCycle:
          file("temphist.txt", "w").write(s)
        else:
          file("temphist.txt", "a").write(s)

    #first time only
    if self.NCycle == 1:
      #write a sequence file
      f = open("seq.txt", "w")
      f.write("\n".join(self.Sim[0].Seq) + "\n")
      f.close()
      #write a temperature file
      Temps = ["%.3f" % s["TEMPSET"] for s in self.Sim]
      file("temps.txt", "w").write("\n".join(Temps) + "\n")
    
    #write the rexstats file
    s = ""
    s += "Number of replicas  : %d\n" % self.NReplica
    s += "Number of cycles    : %d\n" % self.NCycle
    s += "Elapsed time in ps  : %.1f\n" % (self.NCycle * self.Sim[0]["STEPSMD"] * self.Sim[0]["STEPSIZE"])
    s += "Time in ps per cycle: %.1f\n" % (self.Sim[0]["STEPSMD"] * self.Sim[0]["STEPSIZE"])
    s += "Time in ps per frame: %.1f\n" % (self.Sim[0]["STEPSSAVE"] * self.Sim[0]["STEPSIZE"])
    s += "Frames per cycle    : %d\n" % (self.Sim[0]["STEPSMD"] / self.Sim[0]["STEPSSAVE"])
    s += "\nSequence:\n%s\n" % "\n".join(self.Sim[0].Seq)
    Temps = ["%.3f" % Sim["TEMPSET"] for Sim in self.Sim]
    s += "\nTemperatures:\n%s\n" % ("\n".join(Temps),)
    file("rexstats.txt", "w").write(s)
    
    #loop over the user-requested save data
    for k, v in self.SaveDataList.iteritems():
      line = " ".join(["%.4f" % s[k] for s in self.Sim]) + "\n"
      fileMthd(v + fileExt, "a").write(line)
    os.chdir(cwd)

    #update the storage space
    if self.FullUp: self.UpdateStorageSize()

    #update times
    self.SecUpdate = time.time() - self.TimeThisUp
    #estimated time remaining
    self.SecRemainEst = float(self.NCycleTot - self.NCycle) * (self.TimeThisUp - self.TimeLastUp)

    #set the update parameter for the next round
    if time.time() - self.TimeFullUp > FullUpInterval or self.NCycle == self.NCycleTot - 1:
      self.FullUp = True
      self.TimeFullUp = time.time()
      self.NCycleFullUp = self.NCycle + 1
    else:
      self.FullUp = False

    

  def ResetCount(self):
    "Resets move and cycle counts."
    self.MoveAtt = [0.]*len(self.Swaps)
    self.MoveAcc = [0.]*len(self.Swaps)
    self.NCycle = 0
    self.NCycleStart = 0
    self.TimeStart = time.time()
    self.TimeLastUp = time.time()
    self.TimeThisUp = time.time()

  def ResetData(self):
    "Resets the data directory to an initial state."
    #make the data path if it's not there
    if not os.path.isdir(self.DataPath): os.mkdir(self.DataPath)
    cwd = os.getcwd()
    os.chdir(self.DataPath)
    #initialize the statenum.txt and replicanum.txt files
    fileMthd("statenum.txt" + fileExt, "w")
    fileMthd("replicanum.txt" + fileExt, "w")
    for k, v in self.SaveDataList.iteritems():
      fileMthd(v + fileExt, "w")
    #remove any old files
    for spec in ["*mdtrj.crd", "*mdtrj.crd.gz", "*mdene.txt", "*mdene.txt.gz"]:
      for f in glob.glob(spec): os.remove(f)
    os.chdir(cwd)
    

  #SWAP FUNCTIONS

  def SetSwaps(self, Swaps):
    """Sets allowable swaps.
* Swaps is a list of (a,b) tuples where a and b are state indices."""
    #filter the list first
    Swaps = [(min(a,b), max(a,b)) for (a,b) in Swaps]
    #get mapping to old swaps
    l = []
    for (i, s) in enumerate(Swaps):
      if s in self.Swaps:
        l.append((i, self.Swaps.index(s)))
    #update swap variables
    N = len(Swaps)
    MoveAtt, MoveAcc, LastSwapProbs = [0.]*N, [0.]*N, [0.]*N
    for (NewInd, OldInd) in l:
      MoveAtt[NewInd] = self.MoveAtt[OldInd]
      MoveAcc[NewInd] = self.MoveAcc[OldInd]
      LastSwapProbs[NewInd] = self.LastSwapProbs[OldInd]
    #update the arrays
    self.MoveAtt = MoveAtt
    self.MoveAcc = MoveAcc
    self.LastSwapProbs = LastSwapProbs
    self.Swaps = Swaps

  def ReplicaSwap(self, a, b, Pacc = 1., Mode = 0):
    """Swaps configurations for two specified states.  Returns
a Boolean value indicating whether or not swap was accepted.
* a, b: integer indices of the states to swap
* Pacc: float acceptance probability for swap (default is 1.)
* Mode: 0 will perform the swap on disk; 1 will perform the
        swap in memory using Sim.Pos such that calls to
        Sim.SaveData() and Sim.SavePos() are needed later."""
    #get the replica numbers
    ra = self.Replica[a]
    rb = self.Replica[b]
    #form a tuple of the swap attempt
    p = (min([a,b]), max([a,b]))
    #update MoveAtt to reflect the swap attempt
    if p in self.Swaps: self.MoveAtt[self.Swaps.index(p)] += 1.
    #make the move swap according to Pacc
    if Pacc > random.random():
      #update acceptance count
      if p in self.Swaps: self.MoveAcc[self.Swaps.index(p)] += 1.
      #swap replica and state numbers
      self.Replica[a], self.Replica[b] = self.Replica[b], self.Replica[a]
      self.State[ra], self.State[rb] = self.State[rb], self.State[ra]
      #swap the configurations and related data
      mdsim.SwapConfig(self.Sim[a], self.Sim[b], Mode = Mode)
      #swap the restraints if desired
      if self.SwapRest: mdsim.SwapRest(self.Sim[a], self.Sim[b])
      #save the states
      if Mode == 0:
        self.Sim[a].SaveData()
        self.Sim[b].SaveData()
      return True
    else:
      return False
    

  def ReplicaSwapCascadeList(self, PaccFunc, SwapList):
    """From a list of swaps, attempts to swap according to the
list order.  Returns the swap probabilities.
* PaccFunc: function operating on two SimClass objects which gives
  an acceptance probability for a swap move
* SwapList: array of integer tuples (a,b) giving the indices of
  pairs of states which can be swapped
"""
    self.SecSwap = -time.time()
    #first load configurations for mode 1
    if self.SwapMode == 1:
      for s in self.Sim:
        s.LoadPos()
    #copy the swap list
    l = copy.deepcopy(SwapList)
    self.LastSwapProbs = [0.]*len(SwapList)
    while len(l) > 0:
      #remove the first swap
      (a,b) = l.pop(0)
      #get the acceptance probability
      Pacc = PaccFunc(self.Sim[a], self.Sim[b])
      #update the probabilitiy variable
      if (a,b) in self.Swaps:
        ind = self.Swaps.index((a,b))
        self.LastSwapProbs[ind] = Pacc
      #attempt a swap
      Swapped = self.ReplicaSwap(a, b, Pacc, self.SwapMode)
      if Swapped and not self.MultiSwap:
        #if accepted remove all other possible swaps which
        #use either of the just-swapped states
        l = [itm for itm in l if not (a in itm or b in itm)]
    #now save the positions
    if self.SwapMode == 1:
      for s in self.Sim:
        #save positions and class data
        s.SavePos()
        s.SaveData()
        #free up memory
        s.ClearPos()
    self.SecSwap += time.time()

  def ReplicaSwapCascadeUp(self, PaccFunc):
    """Iteratively swaps all replicas starting at first state
and moving up.
* PaccFunc: function operating on two SimClass objects which
  gives an acceptance probability for a swap move"""
    SwapList = copy.deepcopy(self.Swaps)
    self.ReplicaSwapCascadeList(PaccFunc, SwapList) 

  def ReplicaSwapCascadeDn(self, PaccFunc):
    """Iteratively swaps all replicas starting at last state
and moving down.
* PaccFunc: function operating on two SimClass objects which
  gives an acceptance probability for a swap move"""
    SwapList = copy.deepcopy(self.Swaps)
    SwapList.reverse()
    self.ReplicaSwapCascadeList(PaccFunc, SwapList) 

  def ReplicaSwapCascadeAlt(self, PaccFunc):
    """Atlernates swapping all replicas moving up and down.
* PaccFunc: function operating on two SimClass objects which
  gives an acceptance probability for a swap move"""
    if self.LastCascade == 1:
      self.ReplicaSwapCascadeUp(PaccFunc)
    else:
      self.ReplicaSwapCascadeDn(PaccFunc)
    self.LastCascade *= -1

  def ReplicaSwapCascadeRand(self, PaccFunc):
    """Swaps all replicas from a randomly ordered list of all
possible swaps.
* PaccFunc: function operating on two SimClass objects which
  gives an acceptance probability for a swap move"""
    SwapList = copy.deepcopy(self.Swaps)
    random.shuffle(SwapList)
    self.ReplicaSwapCascadeList(PaccFunc, SwapList)
    

  #STATS FUNCTIONS

  def NProcTot(self):
    return len(self.sr.Clients) + 1

  def NProcRun(self):
    return len(self.sr.Clients)

  def StatString(self, Labels = None):
    """Returns a summary string for the replica exchanges.
* Labels: a list of state names
"""
    def Avg(l):
      return sum(l) / float(len(l))
    #if no labels provided, use default
    if Labels is None: Labels = self.Labels
    #if default is none, use temp
    if Labels is None:
      Labels = ["%.fK" % (s["TEMPSET"],) for s in self.Sim]
    #overall acceptance counts
    MoveAttTot = sum(self.MoveAtt)
    MoveAccTot = sum(self.MoveAcc)
    ExFrac = 0.
    if MoveAttTot > 0: ExFrac = MoveAccTot / MoveAttTot
    #calculate elapsed time so far
    ElapTime = (time.time() - self.TimeStart) / 3600.
    #average time per cycle, in min
    AvgCycleTime = (time.time() - self.TimeStart) / (float(max([1, self.NCycle - self.NCycleStart])) * 60.)
    #time last cycle
    LastCycleTime = self.TimeThisUp - self.TimeLastUp
    #fraction of total simulation complete
    CycleFrac = float(self.NCycle) / float(self.NCycleTot)
    #estimated completion time, based on most recent update time
    EstFinishTime = time.localtime((self.TimeThisUp - self.TimeLastUp) \
                                   * float(self.NCycleTot - self.NCycle) \
                                   + time.time())
    #estimated time remaining
    EstTimeRemain = self.SecRemainEst / 3600.
    #average md time
    MDTimes = [s.ElapsedTime() for s in self.Sim]
    AvgMDTime = Avg(MDTimes)
    MaxMDTime = max(MDTimes)
    MinMDTime = min(MDTimes)
    #average write time
    AvgWriteTime = self.SecWrite + Avg([s['CONCATTIME'] for s in self.Sim])
    #client times
    AvgClientTime = Avg([t for t in self.SecClient if not t is None])
    #calculate work and data storage size (estimate based on most recent full update)
    WorkStorage = float(self.WorkStorageSize) / 1048576.
    DataStorage = (float(self.DataStorageSize) / 1048576.) \
                  * float(self.NCycle) / float(max(self.NCycleFullUp, 1))                  
    #write out string
    s = "This summary produced " + time.strftime("%c")
    s += "\nJob started " + time.strftime("%c",time.localtime(self.TimeStart))
    s += "\nBase path: " + os.path.abspath(self.BasePath)
    if self.Restart:
      s += "\nJob restarted from prior run."
    else:
      s += "\nJob started anew."
    s += "\n"
    s += "\nSimulation cycles completed : %.1f %% (%d/%d)" \
         % (CycleFrac*100., self.NCycle, self.NCycleTot)
    s += "\nElapsed/total simtime       : %.1f / %.1f ps" % \
         (self.NCycle * self.Sim[0]["STEPSMD"] * self.Sim[0]["STEPSIZE"],
          self.NCycleTot * self.Sim[0]["STEPSMD"] * self.Sim[0]["STEPSIZE"])
    s += "\nElapsed/remain/tot realtime : %.2f / %.2f / %.2f hrs" \
         % (ElapTime, EstTimeRemain, ElapTime + EstTimeRemain)
    s += "\nEstimated completion date   : " + time.strftime("%c", EstFinishTime)
    s += "\nTime for last cycle         : %.2f sec" % LastCycleTime
    s += "\nAvg times task/up/swap/write: %.2f / %.2f / %.2f / %.2f sec" \
         % (AvgClientTime, self.SecUpdate, self.SecSwap, AvgWriteTime)
    s += "\nMD times min/avg/max        : %.2f / %.2f / %.2f sec" \
         % (MinMDTime, AvgMDTime, MaxMDTime)
    s += "\nClient efficiency           : %.2f%%" % (100. * self.sr.Efficiency())
    s += "\nNum of available clients    : %d" % len(self.sr.Clients)
    s += "\nNumber of replicas          : %d" % self.NReplica
    s += "\n"
    s += "\nCurrent workspace storage requirement : %.3f MB" % WorkStorage
    s += "\nCurrent data storage requirement      : %.3f MB" % DataStorage
    if self.NCycle > 0:
      EstStorage = WorkStorage + DataStorage \
                   * float(self.NCycleTot) / float(max(self.NCycleFullUp, 1))
      s += "\nEstimated total storage requirement   : %.3f MB" % EstStorage
    s += "\n"    
    s += "\nOverall exchanges: %.3f (%d/%d)\n" \
         % (ExFrac, int(MoveAccTot), int(MoveAttTot))
    s += "\n".join(["  %s <> %s :  %.3f (%d/%d)" \
                    % (Labels[self.Swaps[i][0]], Labels[self.Swaps[i][1]],
                       min([self.MoveAcc[i] / max([self.MoveAtt[i], 0.1]), 1.]),
                       int(self.MoveAcc[i]), int(self.MoveAtt[i]))
                    for i in range(0,len(self.Swaps))])
    return s


  #REPLICA EXCHANGE RUNNING

  def GetNextRunStr(self):
    "Generates a string for running on clients."
    #concatenate data from most recent cycle to master files afterwards
    First = (self.NCycle == 1)
    s  = "SIM.RunMD()\n"
    s += "SIM['CONCATTIME'] = -time.time()\n"
    s += "SIM.ConcatData(str(SIM['REXID']), '%s', Gzip = %s, Current = %s, Params = %s)\n" \
         % (self.DataPath, str(Gzip), str(self.FullUp), str(First))
    s += "SIM['CONCATTIME'] = SIM['CONCATTIME'] + time.time()\n"
    return s

  def StartMin(self, ProgFile = None):
    """Starts a minimization run.
* ProgFile: name of the progress file to save (default none)
"""
    self.RexRunning = True
    self.ProgFile = ProgFile
    if not self.ProgFile is None:
      file(ProgFile, "w").write("Running minimization...")    

  def StartRex(self, PaccFunc, NCycle, Labels = None, ProgFile = None,
               PerformSwaps = True, dTemp = None, MultiSwap = None):  
    """Starts a replica exchange run.
* PaccFunc: function operating on two SimClass objects which
  gives an acceptance probability for a swap move
* NCycle: total number of cycles to run
* Labels: text labels for each state (default is current temp)
* ProgFile: name of the progress file to save (default none)
* PerformSwaps: True to perform swaps at updates (default is True)
* dTemp: amount to change each state's temperature at each cycle,
  an array of floats (default is none)
* MultiSwap: False to allow replicas to be swapped only once each
"""
    self.PaccFunc = PaccFunc
    self.NCycleTot = NCycle
    self.ProgFile = ProgFile
    self.PerformSwaps = PerformSwaps
    self.Labels = Labels
    self.RexRunning = True
    if dTemp is None:
      self.SimSetVar("DTEMPSET", [0.]*self.NReplica)
    else:
      self.SimSetVar("DTEMPSET", dTemp)
    if not MultiSwap is None: self.MultiSwap = MultiSwap
    #write out the restraints
    self.RestraintsWrite()

  def RunMin(self, Blocking = False, IntervalInSec = 5.):
    """Runs until minimization finishes.
Returns True if the simulation is done; False otherwise.
* Blocking: True will complete the minimization before
  returning; False will simply run the next set of tasks."""
    if Blocking:
      self.Minimized = False
      if not self.ProgFile is None:
        file(self.ProgFile, "w").write("Running minimization...") 
      self.SimFarm("SIM.RunMin()")
      if not self.ProgFile is None:
        file(self.ProgFile, "w").write("Finished minimization.")
      self.Minimized = True
    else:
      if self.RexRunning and self.TasksExist():
        if self.TasksCheck():
          #remove the done tasks
          self.TasksPurge()
          #produce a progress string and save to file
          if not self.ProgFile is None:
            file(self.ProgFile, "w").write("Finished minimization.")
          self.RexRunning = False
          self.Minimized = True
          return True
        else:
          return False
      else:
        self.RexRunning = True
        self.Minimized = False
        #send out the minimization tasks
        if not self.ProgFile is None:
          file(self.ProgFile, "w").write("Running minimization...")
        self.TasksSend("SIM.RunMin()")
        return False

  def RunRex(self, Blocking = False, IntervalInSec = 0.1):
    """Runs whatever current steps exist in the REX protocol.
Returns True if the REX simulation is done; False otherwise.
* Blocking: True will complete the REX simulation before
  returning; False will simply run the next set of tasks."""
    Done = False
    while not Done:
      #see if there are any tasks in the queue
      if self.TasksExist():
        if self.TasksCheck():
          #remove the done tasks
          self.TasksPurge()
          #all the tasks have finished; time to exchange
          self.Update()
          #perform the swap
          if self.PerformSwaps:
            self.ReplicaSwapCascadeRand(self.PaccFunc)
          #update the temperatures
          for s in self.Sim:
            s["TEMPSET"] = s["TEMPSET"] + s["DTEMPSET"]
          #write the restart data; skip saving sim data since we already did that
          self.RestartWrite(SaveSimData = False)
          #produce a statistics string and save to file
          if not self.ProgFile is None:
            s = self.StatString() + "\n"
            file(self.ProgFile,"w").write(s)
      #see if we are done
      if self.NCycle >= self.NCycleTot:
        self.RestraintsWrite()
        self.RexRunning = False
        if Blocking:
          Done = True
        else:
          return True
      #see if we need to farm more tasks
      if not self.TasksExist():
        RunStr = self.GetNextRunStr()
        self.TasksSend(RunStr, FullTempCopy = False)
      if Blocking:
        #let the socket ring run tasks for at least some time
        self.sr.RunTasks(IntervalInSec)
      else:
        #return a value if we're not blocking
        return False


  #TEMPERATURE OPTIMIZATION

  def TempOptReset(self):
    "Resets temperature optimization variables."
    ######## ADD ALL THIS TO RESTART DATA
    #make a variable containing the average swap prob
    N = len(self.Swaps)
    self.TempOptProb = zeros((2, N), float) + 0.00001
    self.TempOptCount = zeros((2, N), float) + 0.001
    self.TempOptFlag = ones(self.NReplica, int) * -1
    self.TempOptFlag[self.Replica[0]] = 1
    self.TempOptFlag[self.Replica[-1]] = 0
    
  def TempOptCycle(self):
    "Updates the temperature optimization data for a cycle."
    #update the swap probabilities
    for i in range(0, len(self.Swaps)):
      (a,b) = self.Swaps[i]
      Repa, Repb = self.Replica[a], self.Replica[b]
      for Flag in [self.TempOptFlag[Repa], self.TempOptFlag[Repb]]:
        if Flag >= 0:
          self.TempOptProb[Flag, i] += self.LastSwapProbs[i]
          self.TempOptCount[Flag, i] += 1
    #update the extrema flags
    self.TempOptFlag[self.Replica[0]] = 1
    self.TempOptFlag[self.Replica[-1]] = 0
    #now write a file summary if a full update
    if self.FullUp:
      Labels = ["%.1f " % (s["TEMPSET"],) for s in self.Sim]
      frac = self.TempOptFrac()
      s = "Average swap probabilities \n"
      s += "\n".join(["%s <-> %s :  %-8.3f" \
                      % (Labels[self.Swaps[i][0]],
                         Labels[self.Swaps[i][1]],
                         mean(self.TempOptProb[:,i]/self.TempOptCount[:,i]))
                      for i in range(0, len(self.Swaps))])
      s += "\n\nAverage low-temperature fraction\n"
      s += "\n".join(["%s %-8.3f" % (Labels[i], frac[i])
                      for i in range(self.NReplica)])      
      s += "\n"
      #write file; assumes we are in correct folder
      f = file("tempopt.txt","w")
      f.write(s)
      f.close()

  def TempOptFrac(self):
    "Calculates the fractional visit from the temperature opt probabilities."
    f = zeros(self.NReplica, float)
    #make the transition probability matrix
    A = zeros((self.NReplica, self.NReplica), float)
    for i in range(0, len(self.Swaps)):
      (a,b) = self.Swaps[i]
      val = mean(self.TempOptProb[:,i] / self.TempOptCount[:,i])
      A[a,b] = val
      A[b,a] = val
    #fill in the diagonal
    for i in range(0, self.NReplica):
      A[i,i] = 0.
      #make sure probs sum at most to 1
      A[i,:] = A[i,:] / min(1., sum(A[i,:]))
      A[i,i] = 1. - sum(A[i,:])
    #fix the first and last rows
    A[0,0] = 1.
    A[0,1:] = 0.
    A[-1,:] = 0.
    #now find the limiting distribution using the largest eigenvalue;
    #this is just a cheap way based on repeated multiplication
    for i in range(0, 25):
      A = dot(A, A)
    #f is now the first column
    f = A[:,0]
    return f

  def TempOptUpdate(self):
    "Updates the temperatures for optimal exchange and returns a data string."
    #define a function to do linear interpolation
    def Interp(xi, yi, x):
      #find the nearest two x indices and interpolate
      n = len(xi)
      i1, i2 = 0, 1
      #assume monotonic
      for i in range(0,n-1):
        if xi[i] < x and xi[i+1] >= x or xi[i+1] < x and xi[i] >= x:
          i1, i2 = i, i+1
          break
      return yi[i1] + (x - xi[i1]) * (yi[i2] - yi[i1]) / (xi[i2] - xi[i1])
    #calculate the fraction at each temperature
    f = self.TempOptFrac()
    #calculate the ideal delta fraction between temps
    df = 1./float(self.NReplica - 1)
    #calculate new temps at intermediate replicas
    OldTemps = [s["TEMPSET"] for s in self.Sim]
    NewTemps = [Interp(f, OldTemps, df*i) for i in range(self.NReplica-2,0,-1)]
    NewTemps = OldTemps[:1] + NewTemps + OldTemps[-1:]
    #mix with old temps
    NewTemps = [OldTemps[i]*(1.-TempOptMix) + NewTemps[i]*TempOptMix
                for i in range(0,self.NReplica)]
    #put the temperatures in the replicas
    self.SimSetVar("TEMPSET", NewTemps)
    #reset the optimization variables
    self.TempOptReset()
    #produce a data string
    s = "%-8s %-8s %-8s\n" % ("OldTemp", "NewTemp", "f")
    s += "\n".join(["%-8.3f %-8.3f %-8.3f" % (OldTemps[i], NewTemps[i], f[i]) \
                    for i in range(0, self.NReplica)])
    s += "\n"
    return s

    
    
