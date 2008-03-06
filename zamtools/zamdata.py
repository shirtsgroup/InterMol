#!/usr/bin/env python

#LAST MODIFIED: 05-02-07

from numpy import *
import os, sys, copy, shutil, glob, cPickle, gzip
import sequence, protein

#globals
PairMinCO = 3     #Minimum CO for pair analyses
PairMaxCO = 1000  #Maximum CO for pair analyses
PairMinECO = 3    #Minimum ECO for considering a pair for analysis
AllPairsDflt = True  #analyze all (true) or just hydrophobic (false) pairs
MinScore = -1.e200


#VALUES FOR Frag.Stage:
# "grow" =  frag is new or has been grown from a previous frag
#           and is flagged for simulation
# "assem" = frag has been assembled from previous frags and is
#           flagged for simulation
# "done" =  frag has finished being simulated; may be used for
#           future growth/assembly

#CONVENTIONS:
#-note that all FragDataClass member functions return information
# based on fragment prior restraints only (not global restraints),
# whereas corresponding functions in zamDataClass return information
# based on purely global restraints
#-internally, all atom numbers and residue numbers start at 0,
# but these are changed for output to base 1


def NumToAlpha(Num):
  "Converts a number to an alpha string."
  if Num <= 0: return ""
  Alpha = "abcdefghijklmnopqrstuvwxyz"
  n = Num
  r = ""
  while True:
    n, dig = divmod(n-1, len(Alpha))
    r = Alpha[dig] + r
    if n == 0:
      return r

def CalcECO(a, b, RestList, MaxECO = None, MaxDepth = 3,
            CurDepth = 0, Choices = None):
  "Returns the shortest route from residue nums a to b, incl restraints."
  a, b = min(a,b), max(a,b)
  if b - a <= 1: return b - a
  if not MaxECO is None: MaxDepth = min(MaxECO, MaxDepth)
  if CurDepth >= MaxDepth: return b - a
  if Choices is None:
    #start off with unique restraints
    Choices = []
    for (c,d) in RestList:
      c, d = min(c,d), max(c,d)
      if not (c,d) in Choices: Choices.append((c,d))
  #check for no choices
  if len(Choices) == 0: return b - a
  #recurse through the choices
  MinDist = b - a
  if not MaxECO is None: MinDist = min(MinDist, MaxECO)
  for Pair in Choices:
    (c, d) = Pair
    if abs(d - a) < abs(c - a): (c,d) = (d,c)
    Dist = abs(c - a) + 1
    if not MaxECO is None and Dist > MaxECO: continue
    DoneRange = range(min(c,a) + 1, max(c,a))
    NewChoices = [(x,y) for (x,y) in Choices if not (x,y) == Pair
                  and not x in DoneRange and not y in DoneRange]
    Dist += CalcECO(d, b, RestList, MaxECO, MaxDepth, CurDepth + 1, NewChoices)
    MinDist = min(MinDist, Dist)
  return MinDist

def CalcResWithinECO(a, b, RestList, MaxECO, MaxDepth = 1):
  "Returns true if a residue is within MaxECO."
  #shortcut for MaxDepth == 1
  if MaxDepth == 1:
    (a,b) = min(a,b), max(a,b)
    for (c,d) in RestList:
      (c,d) = min(c,d), max(c,d)
      if abs(a-c) + abs(b-d) + 1 <= MaxECO: return True
    return False
  else:
    return CalcECO(a, b, RestList, MaxECO + 1, MaxDepth) <= MaxECO
  


class SnapClass:
  
  def __init__(self, Pdb = None, Parent = (None,0),
               StartRes = 0, StopRes = 0):
    self.Pdb = Pdb            #Pdb file (relative to master base path) [string]
    self.Parent = Parent      #source fragment and structure index [(FragDataClass, Index)]
    self.StartRes = StartRes  #starting residue number [int]
    self.StopRes = StopRes    #stopping residue number [int]
    self.Contacts = []        #residue-residue contacts [list of (a,b) tuples]
    self.AvgCO = 0.
    if not Pdb is None: self.GetContacts()
    self.Score = MinScore*len(self.Contacts)     #score, based on contacts
    if not self.Parent == (None,0):
      Frag, ind = self.Parent
      self.Name = "%s:%d" % (Frag.Name(), ind)
    else:
      self.Name = ""

  def GetContacts(self):
    """Gets the contacts from the Pdb file."""
    p = protein.ProteinClass(Pdb = self.Pdb).Decap()
    self.Contacts = [(a + self.StartRes, b + self.StartRes)
                     for (a,b) in p.ResContactList()]
    m = len(self.Contacts)
    if m > 0:
      self.AvgCO = float(sum([abs(a-b) for (a,b) in self.Contacts])) / m
                     
  def UpdateScore(self, ContactScores):
    """Updates the structure score based on the contact scores."""
    self.Score = sum([ContactScores[a,b] for (a,b) in self.Contacts])

  def ContainsSnapRes(self, Snap):
    """Returns True if Snap is inside the residue range of self."""
    return self.StartRes <= Snap.StartRes and self.StopRes >= Snap.StartRes

  def ContainsSnap(self, Snap):
    """Returns True if all residues and contacts of Snap are contained within self."""
    if not self.ContainsSnapRes(Snap): return False
    for (a,b) in Snap.Contacts:
      if not (a,b) in self.Contacts: return False
    return True

  def SameContacts(self, Snap):
    """Returns True if Snap has the same contacts as Snap."""
    return self.ContainsSnap(Snap) and len(Snap.Contacts) == len(self.Contacts)

  def SubScore(self, ResNumList, ContactScores):
    """Returns the score only from contacts among residues in ResNumList."""
    return sum([ContactScores[a,b] for (a,b) in self.Contacts
                if a in ResNumList and b in ResNumList])


def SnapPairScore(Snap1, Snap2, ContactScores, COWeight = 1.):
  """Returns the sum scores of Snap1 and Snap2, where contacts in regions
  of overlap are averaged."""
  Snap1Res = range(Snap1.StartRes, Snap1.StopRes + 1)
  Snap2Res = range(Snap2.StartRes, Snap2.StopRes + 1)
  LoopRes = [x for x in Snap1Res if x in Snap2Res]
  Score = 0.
  for (a,b) in Snap1.Contacts + Snap2.Contacts:
    if a in LoopRes or b in LoopRes:
      Score += 0.5 * ContactScores[a,b]
    else:
      Score += ContactScores[a,b]
  Scale = (Snap1.StopRes - Snap1.StartRes + Snap2.StopRes - Snap2.StartRes + 2)
  Scale = (float(Scale) - len(LoopRes)) / float(Scale)
  Score += COWeight * Scale * (Snap1.AvgCO + Snap2.AvgCO)
  return Score


class SnapDatabaseClass:
  
  def __init__(self, zd, ContactScores = None, RemoveRedundant = False):
    "Searches through all fragments and gets snaps."
    self.Snaps = []
    for Frag in zd.Frags: self.AddSnaps(Frag)
    if not ContactScores is None: self.UpdateScores(ContactScores)
    if RemoveRedundant: self.RemoveRedundant()

  def __len__(self):
    return len(self.Snaps)

  def AddSnaps(self, Frag):
    "Adds a fragment's snaps to the database."
    SnapPdbs = [os.path.join(Frag.BasePath, "conf/", Pdb)
                for Pdb in Frag.SnapPdbs]
    for i in range(0, len(SnapPdbs)):
      Pdb = SnapPdbs[i]
      StartRes, StopRes = Frag.SnapResRange[i]
      Snap = SnapClass(Pdb = Pdb, Parent = (Frag,i),
                       StartRes = StartRes, StopRes = StopRes)
      self.Snaps.append(Snap)  

  def GetSnaps(self, ParentFrag):
    "Returns snaps with a given parent fragment."
    Snaps = [s for s in self.Snaps if s.Parent[0] == ParentFrag]
    Snaps.sort(cmp = lambda x, y: cmp(x.Parent[1], y.Parent[1]))
    return Snaps

  def UpdateScores(self, ContactScores):
    "Updates the snapshot scores based on contact scores."
    for Snap in self.Snaps:
      Snap.UpdateScore(ContactScores)

  def RemoveRedundant(self):
    """Removes redundant snapshots."""
    NewSnaps = []
    for Snap in self.Snaps:
      Same = False
      for s in NewSnaps:
        if s.SameContacts(Snap):
          Same = True
          break
      if not Same:
        NewSnaps.append(Snap)
    print "Removed %d redundant snapshots." % (len(self.Snaps) - len(NewSnaps))
    self.Snaps = NewSnaps


class FragDataClass:

  #Don't change
  AnalysisVars = ["MesoEntropy", "MesoGroundPop", "PunchStabTot", "PunchCoopTot",
                  "PairList", "CutProb", "Pmf", "PmfDist", "PmfScore", "BestnCont",
                  "PunchStab", "PunchCoop", "ClustPop"]
  PickleExcludeVars = ["rx", "Task"]
  DataExcludeVars = ["Parents", "Children"]  

  def __init__(self, StartRes = 0, StopRes = 0, Seq = None, Parents = []):
    #variables managed by this module:
    StartRes, StopRes = sorted([StartRes,StopRes])
    self.StartRes = StartRes    #starting residue number [int]
    self.StopRes = StopRes      #stopping residue number [int]
    self.Ver = 1                #version number, in case other fragments exist [int]
    self.Seq = Seq              #list of 3-letter sequence codes [list of strings]
    self.PriorRest = []         #list of restraints present upon fragment creation
                                #[list of (int,int) tuples where int is residue num]
    self.BasePath = None        #base path for this fragment [string]
    self.AllPairs = AllPairsDflt#consider all (True) or just hydrophobic (False) residue pairs [bool]
    self.RestPairs = True       #run analysis on restrained pairs too? [bool]
    self.Note = ""              #user note for reference [string]
    self.CapN = True            #cap the n-terminal of this fragment? [bool]
    self.CapC = True            #cap the c-terminal of this fragment? [bool]
    
    #variables managed by zamchoice module:
    self.Stage = ""             #stage of fragment [string]
    self.Parents = Parents      #list of parent fragments [list of (FragDataClass, Index)
                                #where the indices are the snapshot numbers or -1 for all
    self.Children = []          #list of child fragments [list of (FragDataClass, Index)
                                #where the indices are the snapshot numbers or -1 for all

    #variable managed by the zamcontact module    
    self.FragScore = 0.         #score based on fragment data [float]
    
    #variables managed by zamrex module:
    self.NRep = 0               #number of replicas to use during simulation [int]
    self.NFrame = 0             #total number of frames collected during last run [int]
    self.FrameTime = 0.         #time corresponding to one frame in ps [float]
    self.NCycle = 0             #total number of REX cycles [int]
    self.CycleTime = 0.         #time corresponding to one REX cycle in ps [float]
    self.RepelSalt = True       #True will add a repulsion between potential salt bridges [bool]
    self.ResRigid = []          #residues which are to be held rigid by phi-psi rest [list of int]
    self.ResAnchorBB = []       #list of residues whose backbone to anchor to space [list of int]
    self.ResAnchorAA = []       #list of residues whose all atoms to anchor to space [list of int]
    self.SS = ""                #secondary structure [string]
    self.ScaleRest = False      #scale restraint strength from 0 to 1 over replicas
    self.UmbrRest = []          #residues to umbrella sample [list of (a,b) tuples]
    self.UmbrDist = []          #distances in umbrella sampling [list of float]

    #variables managed by zamrun module:
    self.RexTime = 0.           #time in ps to run REX [float]
    self.AnalTime = 0.          #time in ps to use for analysis [float]
    self.SkipTime = -1.         #time in ps to skip before analysis; -1 puts at end [float]
    self.CalcRexTime = True     #whether or not to calculate the rex time [bool]
    self.CalcAnalTime = True    #whether or not to calculate the analysis time [bool]
    self.RexStage = 0           #used by grow and assem to remember where they are [int]
    self.MakeInitConf = True    #whether or not to make the initial conformations [bool]
    self.CalcSwapRest = False   #whether or not to automatically gen swapped rest [bool]
    self.CalcSwapRestTime = False#whether or not to automatically calc swap rest time [bool]
    self.RexTimeSwapRest = 0.   #time in ps (out of RexTime) to apply and swap different
                                #restraints in each replica [float]
    self.SwapRest = []          #list of swapped restraints, one per replica (in addition to
                                #Frag.PriorRest restraints) [list of (int,int) tuples]

    #variables managed by zamanal module:
    self.SnapPdbs = []          #snapshot pdbs [list of strings]
    self.SnapResRange = []      #residue range of the snapshot pdbs [list of (StartRes, StopRes)]
    self.SnapTrimEnds = False   #trim the snaps to remove loose ends [bool]
    self.SnapTrimSS = False     #remove snapshots with redundant secondary structure [bool]
    self.SnapWeight = []        #frac weight of each snapshot (cluster percentage) [list of float]
    self.MesoEntropy = 0.       #the mesostring entropy for this fragment [float]
    self.MesoGroundPop = 0.     #the fraction of the time spent in the ground mesostate,
                                #i.e., the most populous mesostate [float]
    self.PunchStabTot = 0.      #total stability assessed by punch analysis [float]
    self.PunchCoopTot = 0.      #total cooperativity assessed by punch analysis [float]
    self.RunPunch = True        #run the punch analysis? [bool]
    #the following are cleared after restraints are decided (to save space):
    self.PairList = []          #list of residue pairs to monitor for contacts
                                #[list of (int,int) tuples where int is residue num]
    self.CutProb = []           #cutoff probabilities for contacts [list of floats]
    self.Pmf = []               #list of pmf arrays [list of arrays of floats]
    self.PmfDist = []           #pmf distances [list of floats]
    self.PmfScore = []          #scores for pmf/contact prob [list of floats]
    self.PunchStab = []         #scores for punch stability [list of floats]
    self.PunchCoop = []         #scores for punch cooperativity [list of floats]
    self.BestnCont = []         #list of scores and best contacts for restraints [list of (float,int,int) tuples]
    self.ClustPop = []          #list of fractional cluster population [list of floats]

    #variables managed by zamrest module
    self.PropRest = []          #list of restraints proposed after fragment simulation;
                                #[list of (float,int,int) tuples where float is a
                                #confidence score and int is residue num]

    #variables managed by main zam module
    self.RestartInd = 0         #restart index for the current run [int]
    self.rx = None              #replica exchange class
    self.Task = None            #socketring task or list of tasks
    self.ProgPrefix = ""        #prefix for progress file updates [string]

    #update the parents
    self.UpdateParents()

  def __len__(self):
    return Frag.StopRes - Frag.StartRes + 1

  def __repr__(self):
    if len(self.Seq) > 0:
      return "(Frag %s : %s%d-%s%d)" % (self.Name(), self.Seq[0],
             self.StartRes+1, self.Seq[-1], self.StopRes+1)
    else:
      return "(Frag %s)" % self.Name()

  def __setattr__(self, key, val):
    """Sets a class attribute."""
    #set the value
    self.__dict__[key] = val
    #check the caps
    if key == "CapN" or key == "CapC" or key == "Seq":
      if self.__dict__.get("Seq", 0) > 0:
        a = self.Seq[0]
        if sequence.Terminal(a) or sequence.Cap(a):
          self.__dict__["CapN"] = False
        a = self.Seq[-1]
        if sequence.Terminal(a) or sequence.Cap(a):
          self.__dict__["CapC"] = False

  def UpdateVer(self):
    "Updates old versions to the latest FragDataClass."
    #first check snapshot
    if not "SnapPdbs" in self.__dict__:
      self.SnapPdbs = glob.glob(os.path.join(self.BasePath, "conf/snap*.pdb"))
      self.SnapResRange = [(self.StartRes, self.StopRes)]*len(self.SnapPdbs)
    #now add any other missing data elements
    Dummy = FragDataClass()
    for k, v in Dummy.__dict__.iteritems():
      if not k in self.__dict__:
        self.__dict__[k] = v
    #check cap
    if "Cap" in self.__dict__:
      self.CapN, self.CapC = self.Cap, self.Cap
    #update Parents and Children to tuples
    for i in range(0, len(self.Parents)):
      x = self.Parents[i]
      if not type(x) is tuple:
        self.Parents[i] = (x, -1)
    for i in range(0, len(self.Children)):
      x = self.Children[i]
      if not type(x) is tuple:
        self.Children[i] = (x, -1)    

  def UpdateParents(self):
    "Updates the parent fragments to reflect the current child."
    for (Frag, Ind) in self.Parents:
      if not (self, Ind) in Frag.Children:
        Frag.Children.append((self, Ind))

  def Unlink(self):
    "Unlinks from parent and child fragments."
    for (Frag, ind) in self.Parents:
      Frag.Children = [(f,i) for (f,i) in Frag.Children if not f is self]
    for (Frag, ind) in self.Children:
      Frag.Parents = [(f,i) for (f,i) in Frag.Parents if not f is self]

  def GetAllChildren(self):
    "Returns a list of all child fragments."
    l = copy.copy(self.Children)
    for (f,i) in l:
      l.extend(f.GetAllChildren())
    return l

  def GetAllParents(self):
    "Returns a list of all parent fragmentss."
    l = copy.copy(self.Parents)
    for (f,i) in l:
      l.extend(f.GetAllParents())
    return l

  def Name(self):
    """Returns a string name of the fragment.
    Version 0 has no suffix, all other have alpha suffices."""
    s = "%d-%d" % (self.StartRes+1, self.StopRes+1)
    if self.Ver > 0:
      s += NumToAlpha(self.Ver)
    return s

  def MakePath(self, MainPath = ".", AllPaths = False):
    "Makes the directory for this fragment."
    self.BasePath = "r%s/" % self.Name()
    self.BasePath = os.path.join(MainPath, self.BasePath)
    if not os.path.isdir(self.BasePath):
      os.mkdir(self.BasePath)
      if AllPaths:
        for SubPath in ["anal","conf","data"]:
          p = os.path.join(self.BasePath, SubPath)
          if not os.path.isdir(p): os.mkdir(p)

  def DelPath(self):
    "Removes the path and all data for this fragment."
    if os.path.isdir(self.BasePath):
      shutil.rmtree(self.BasePath)
    self.SnapPdbs, self.SnapResRange = [], []

  def DelAnalysis(self):
    "Removes all the analysis files."
    AnalPath = os.path.join(self.BasePath, "anal")
    if os.path.isdir(AnalPath): shutil.rmtree(AnalPath)
    Snaps = glob.glob(os.path.join(self.BasePath, "conf/snap*.pdb"))
    for fn in Snaps: os.remove(fn)
    self.SnapPdbs, self.SnapResRange = [], []

  def DelData(self, Work = True, Data = True, Conf = True):
    "Removes all the data, work, and initial configuration files."
    if Data:
      DataPath = os.path.join(self.BasePath, "data")
      if os.path.isdir(DataPath): shutil.rmtree(DataPath)
    if Work:
      WorkPath = os.path.join(self.BasePath, "work")
      if os.path.isdir(WorkPath): shutil.rmtree(WorkPath)
    if Conf:
      InitPdbs = glob.glob(os.path.join(self.BasePath, "conf/init*.pdb"))
      InitPdbs += glob.glob(os.path.join(self.BasePath, "conf/refinit*.pdb"))
      for fn in InitPdbs: os.remove(fn)      

  def GetOverlapRes(self, ResList):
    "Returns residues in ResList that are in the current fragment."
    return [a for a in ResList if a >= self.StartRes and a <= self.StopRes]
  
  def AddSS(self, ResList, SS):
    "Gives residues in ResList secondary structure restraints for SS."
    ResList = self.GetOverlapRes(ResList)
    if len(ResList) == 0: return
    n = self.StopRes - self.StartRes + 1
    if len(self.SS) == 0: self.SS = "-"*n
    if SS in " X": SS = "-"
    for i in ResList:
      j = i - self.StartRes
      self.SS = self.SS[:j] + SS[0:1] + self.SS[j+1:]
    if self.SS.count("-") == n: self.SS = ""

  def RestInRange(self, a, b):
    "Indicates whether or not the restraint (a,b) is within the fragment."
    return a >= self.StartRes and a <= self.StopRes and b >= self.StartRes and b <= self.StopRes
  
  def AddRest(self, RestList):
    "Adds relevant restraints to a fragment."
    for (a,b) in RestList:
      if self.RestInRange(a,b):
        a,b = min(a,b),max(a,b)
        if not (a,b) in self.PriorRest:
          self.PriorRest.append((a,b))

  def InPropRest(self, a, b):
    "Tells whether or not a restraint is among those proposed."
    return (a,b) in [(c,d) for (Score,c,d) in self.PropRest if Score > 0.]

  def ResConnected(self, a, b, OtherRest = []):
    "Tells if residues a and b are connected on backbone or by restraint."
    return abs(a-b) == 1 or (min(a,b), max(a,b)) in (self.PriorRest + OtherRest)

  def ResCO(self, a, b):
    "Returns contact order from residue nums a to b."
    return abs(a-b)

  def ResECO(self, a, b, OtherRest = [], MaxECO = None, MaxDepth = 10):
    "Returns the shortest route from residue nums a to b, incl restraints."
    return CalcECO(a, b, self.PriorRest + OtherRest, MaxECO, MaxDepth)
  
  def ResWithinECO(self, a, b, OtherRest = [], MaxECO = PairMinECO - 1, MaxDepth = 1):
    "Tells if a pair of residues is retrained within an ECO."
    return CalcResWithinECO(a, b, self.PriorRest + OtherRest, MaxECO, MaxDepth)

  def ResDOF(self):
    "Returns approx number of backbone degrees of freedom, incl restraints."
    n = 2*abs(self.StopRes - self.StartRes)
    #subtract rigid segments
    l = self.ResRigid + self.ResAnchorAA + self.ResAnchorBB
    l = [x for (i,x) in enumerate(l) if not x in l[i+1:]]
    n -= 2*len(l)
    #subtract prior restraints
    n -= len(self.PriorRest)
    return n

  def ContainsSeg(self, a, b):
    "Tells if a residue segment is contained within the current fragment."
    ResRange = range(self.StartRes, self.StopRes+1)
    if a in ResRange and b in ResRange:
      return True
    else:
      return False

  def ContainsFrag(self, Frag):
    "Tells if another fragment is contained within the current one."
    return self.ContainsSeg(Frag.StartRes, Frag.StopRes)

  def OverlapsSeg(self, a, b):
    "Tells if a residue segment overlaps with the current fragment."
    if (a < self.StartRes and b < self.StartRes) or (a > self.StopRes and b > self.StopRes):
      return False
    else:
      return True

  def OverlapsFrag(self, Frag):
    "Tells if another fragment overlaps with the current fragment."
    return self.OverlapsSeg(Frag.StartRes, Frag.StopRes)

  def SegDist(self, a, b):
    "Returns the number of residues between self and a segment."
    Dist1 = self.StartRes - b
    Dist2 = a - self.StopRes
    return max(Dist1, Dist2)
  
  def FragDist(self, Frag):
    "Returns the number of residues between self and Frag."
    return self.SegDist(Frag.StartRes, Frag.StopRes)
    
  def GetPairList(self, MinCO = PairMinCO, MaxCO = PairMaxCO, MinECO = PairMinECO):
    "Returns a list of pairs for monitoring contacts."
    PairList = []
    for gap in range(MinCO, MaxCO + 1):
      for i in range(self.StartRes, self.StopRes-gap+1):
        #skip pair if already restrained to be close
        if not self.RestPairs:
          if self.ResWithinECO(i, i+gap, MaxECO = MinECO): continue
        #check if pair is attractive
        if self.AllPairs or sequence.HydrophobicPair(self.Seq[i-self.StartRes], self.Seq[i+gap-self.StartRes]):
          PairList.append((i,i+gap))
    return PairList
  
  def NResOverlap(self, a, b):
    "Tells the number of overlapping residues with the segment a,b."
    return len([i for i in range(self.StartRes, self.StopRes+1) if i in range(a,b+1)])

  def NResAdd(self):
    "Tells the number of residues added relative to the parents."
    InParents = [False]*(self.StopRes - self.StartRes + 1)
    for (f,i) in self.Parents:
      for i in range(max(f.StartRes, self.StartRes), min(f.StopRes, self.StopRes)):
        InParents[i - self.StartRes] = True
    return InParents.count(False)

  def ParentsIn(self, Parents):
    """Returns True if Parents contains all the parents of self.
    False is returned if self has no parents."""
    if len(self.Parents) == 0:
      return False
    else:
      return [p in Parents for p in self.Parents].count(False) == 0

  #DATA PACKING

  def __getstate__(self):
    """Returns a dictionary for pickling."""
    NewDict = self.__dict__.copy()
    for v in FragDataClass.PickleExcludeVars:
      if v in NewDict: del NewDict[v]
    return NewDict

  def DumpsData(self, Exclude = []):
    "Puts all data except linked fragments into a pickle string."
    NewDict = self.__dict__.copy()
    #remove excluded variables
    Exclude = Exclude + FragDataClass.DataExcludeVars + FragDataClass.PickleExcludeVars
    for v in Exclude:
      if v in NewDict: del NewDict[v]
    return cPickle.dumps(NewDict, cPickle.HIGHEST_PROTOCOL)

  def LoadsData(self, s):
    "Updates data with data from a pickled string."
    NewDict = cPickle.loads(s)
    self.__dict__.update(NewDict)

  def SaveAnalData(self):
    "Saves all analysis data into a pickled file for easy retrieval."
    AnalPath = os.path.join(self.BasePath, "anal/")
    if not os.path.isdir(AnalPath): os.mkdir(AnalPath)
    dat = dict([(k,v) for (k,v) in self.__dict__.iteritems()
                if k in FragDataClass.AnalysisVars])
    cPickle.dump(dat, gzip.GzipFile(os.path.join(AnalPath, "fraganal.dat.gz"), "w"))

  def LoadAnalData(self):
    "Loads all analysis data from a pickled file."
    f = os.path.join(self.BasePath, "anal/fraganal.dat.gz")
    if os.path.isfile(f):
      dat = cPickle.load(gzip.GzipFile(f,"r"))
      self.__dict__.update(dat)
      return
    f = os.path.join(self.BasePath, "anal/fraganal.dat")
    if os.path.isfile(f):
      dat = cPickle.load(file(f,"r"))
      self.__dict__.update(dat)
      return

  def ClearAnalData(self):
    "Clears all analysis data."
    for k in FragDataClass.AnalysisVars:
      if type(self.__dict__[k]) is list:
        self.__dict__[k] = []
      elif type(self.__dict__[k]) is float:
        self.__dict__[k] = 0.
      elif type(self.__dict__[k]) is int:
        self.__dict__[k] = 0       
      else:
        self.__dict__[k] = []


class ZamDataClass:

  def __init__(self, Seq = []):  
    self.Seq = sequence.SeqToList(Seq)         #list of 3-letter sequence codes [list of strings]
    self.Frags = []                            #list of fragments [list of FragDataClass instances]
    self.Rest = []                             #master list of restraints [list of (int,int) tuples
                                               # where int is residue num]
    self.PreBuildString = ""                   #string run in tleap before system is built [string]
    self.PostBuildString = ""                  #string run in tleap after system is built [string]
    self.UpdateFrags = False                   #indicates whether or not to automatically create new fragments [boolean]
    self.UpdateRest = False                    #whether or not to update master restraints automatically [boolean]
    self.Stage = ""                            #the stage of the zam simulation [string]
    self.Verbose = False                       #verbose mode for decision making [bool]
    self.SS = ""                               #secondary structure [string]
    self.AssembleByFrag = True                 #False to assemble by individual snapshots instead of whole fragments at a time [bool]

  def UpdateVer(self):
    "Updates old versions to the latest ZamDataClass."
    #add any other missing data elements
    Dummy = ZamDataClass()
    for k, v in Dummy.__dict__.iteritems():
      if not k in self.__dict__:
        self.__dict__[k] = v
    for Frag in self.Frags: Frag.UpdateVer()

  def SetSeq(self, Seq):
    "Sets the sequence."
    self.Seq = sequence.SeqToList(Seq)

  def AddFrag(self, StartRes, StopRes, Parents = [], Stage = "grow", Versioned = False):
    "Adds a fragment to the zam simulation and returns it."
    Seq = self.Seq[StartRes:StopRes+1]
    NewFrag = FragDataClass(StartRes, StopRes, Seq, Parents)
    #first add the master restraints and secondary structure
    NewFrag.AddRest(self.Rest)
    if len(self.SS) > 0:
      NewFrag.SS = self.SS[StartRes:StopRes+1]
      if NewFrag.SS.replace("-","") == "": NewFrag.SS = ""
    NewFrag.Stage = Stage
    #check to see if this frag already exists
    SameFrags = self.GetFrags(StartRes = StartRes, StopRes = StopRes)
    if len(SameFrags) > 0:
      #increment the version number
      NewFrag.Ver = max([f.Ver for f in SameFrags]) + 1
      if Versioned:
        NewFrag.Ver = max(NewFrag.Ver, 1)
    self.Frags.append(NewFrag)
    return NewFrag

  def DelFrag(self, FragNum):
    "Deletes a fragment from the zam data structure."
    Frag = self.Frags[FragNum]
    Frag.Unlink()
    del self.Frags[FragNum]

  def AddSS(self, ResList, SS):
    "Gives residues in ResList secondary structure restraints for SS."
    if len(ResList) == 0: return
    n = len(self.Seq)
    if len(self.SS) == 0: self.SS = "-"*n
    if SS in " X": SS = "-"
    for i in ResList:
      self.SS = self.SS[:i] + SS[0:1] + self.SS[i+1:]
    if self.SS.count("-") == n: self.SS = ""

  def AddRest(self, a, b, MinECO = PairMinECO, Force = False):
    "Adds a restraint to the zam simulation."
    (a,b) = (min(a,b), max(a,b))
    if (a,b) in self.Rest: return
    if Force or not self.ResWithinECO(a, b, MaxECO = MinECO):
      self.Rest.append((a,b))
    
  def ResConnected(self, a, b, OtherRest = []):
    "Tells if residues a and b are connected on backbone or by restraint."
    return abs(a-b) == 1 or (min(a,b), max(a,b)) in (self.Rest + OtherRest)

  def ResCO(self, a, b):
    "Returns contact order from residue nums a to b."
    return abs(a-b)

  def ResECO(self, a, b, OtherRest = [], MaxECO = None, MaxDepth = 10):
    "Returns the shortest route from residue nums a to b, incl restraints."
    return CalcECO(a, b, self.Rest + OtherRest, MaxECO, MaxDepth)

  def ResWithinECO(self, a, b, OtherRest = [], MaxECO = PairMinECO - 1, MaxDepth = 1):
    "Tells if a pair of residues is retrained within a minimum ECO."
    return CalcResWithinECO(a, b, self.Rest + OtherRest, MaxECO, MaxDepth)

  def ResDOF(self, a, b):
    "Returns approx number of backbone degrees of freedom, incl restraints."
    n = 2*abs(b-a)-2
    for (c,d) in self.Rest:
      if c >= a and c <= b and d >= a  and d <= b:
        n -= 1
    return n
 
  def GetFrags(self, Stage = [], StartRes = None, StopRes = None,
               Length = None):
    """Returns a list of fragments, optionally with values in Stage
    or with given starting and/or stopping residues."""
    if len(Stage) == 0:
      FragList = copy.copy(self.Frags)
    else:
      FragList = [f for f in self.Frags if f.Stage in Stage]
    if not StartRes is None:
      FragList = [f for f in FragList if f.StartRes == StartRes]
    if not StopRes is None:
      FragList = [f for f in FragList if f.StopRes == StopRes]
    if not Length is None:
      FragList = [f for f in FragList if f.StopRes - f.StartRes + 1 == Length]
    return FragList

  def GetOrderedFrags(self, Stage = []):
    "Returns a list of fragments in order of starting residue."
    FragList = [(f.StartRes, f) for f in self.GetFrags(Stage)]
    FragList.sort()
    return [f for (StartRes, f) in FragList]

  def GetRunFrags(self):
    "Returns a list of fragments to be simulated (with stage grow or assem)."
    return self.GetFrags(["grow", "assem"])

  def GetDoneFrags(self):
    "Returns a list of fragments that have been simulated (with stage done)."
    return self.GetFrags(["done"])

  def AnyFrags(self, Stage = []):
    "Returns True if there are any fragments with specified stages."
    return len(self.GetFrags(Stage)) > 0

  def AnyRunFrags(self):
    "Returns True if there are any fragments to be simulated (with stage grow or assem)."
    return len(self.GetRunFrags()) > 0

  def SegInFrags(self, a, b, Stage = []):
    "Returns True if a segment is an existing fragment."
    return len(self.GetFrags(Stage, StartRes=a, StopRes=b)) > 0

  def SegWithinFrags(self, a, b, Stage = []):
    "Returns True if a segment is within an existing fragment."
    return [f.ContainsSeg(a,b) for f in self.GetFrags(Stage)].count(True) > 0  

  def SegOverlapsFrags(self, a, b, Stage = []):
    "Returns True if a segment overlaps with an existing fragment."
    return [f.OverlapsSeg(a,b) for f in self.GetFrags(Stage)].count(True) > 0

  def GetFragsWithParents(self, Parents):
    "Returns fragments whose parents are contained in Parents."
    return [f[0] for f in Parents[0][0].Children if f[0].ParentsIn(Parents)]

  def AnyFragsWithParents(self, Parents):
    "Returns True if fragments exist whose parents are contained in Parents."
    return len(self.GetFragsWithParents(Parents)) > 0
  