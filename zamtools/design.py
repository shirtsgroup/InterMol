#!/usr/bin/env python

#LAST UPDATED: 02-06-07

from numpy import *
import sys, os, shutil, random, tempfile, copy, time, math, cPickle, glob
import mdsim, protein, sequence, pdbtools, proteinalg
import socketring


TempDirBase = ""
if len(sys.argv) > 1: TempDirBase = os.path.abspath(sys.argv[1])


#GLOBALS
TargetPdb = "target.pdb"
CEFile = "contactenergies.dat"
MutTemp = 0.005
NMut = 25000
NBest = 10
RerunFreq = 0.10
HybridBestFreq = 0.10
SeqTol = 0.4
RandomStartSeq = True
SkipResStart = ['PHE', 'TYR', 'TRP', 'HIS', 'PRO', 'CYS']
SkipResMut = ['PRO', 'CYS', 'HIS']
AllResMut = ['ALA', 'CYS', 'ASP', 'GLU', 'GLY', 'PHE', 'HIS', 'ILE', 'LYS', 'LEU',
             'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']


#MD PROTOCOL
#these are list of (NStepsMin1, NStepsMin2, NStepsMD, Stepsize, Temp1, Temp2,
#                   FConst, NBScale1, NBScale2, RadScale1, RadScale2, VelLimit)
MDStages = [(400, 0, 1000, 0.001, 270, 270, 1.0, 1.0, 1.0, 1.0, 1.0, 40),
            (  0, 0, 1000, 0.001, 270, 270, 1.0, 1.0, 1.0, 1.0, 1.0, 40)]


#database type
class BestClass:
  def __init__(self, Seq = [], Score = 1.e200, Prot = None):
    self.Score = Score
    self.Seq = copy.deepcopy(Seq)
    if Prot is None:
      self.Prot = None
    else:
      self.Prot = Prot.Copy()
      if not self.Prot.Seq == self.Seq:
        self.Prot.MutateSeq(self.Seq)


#RESTART DATA
RestartData = ["MutStart", "MutDone", "BestList", "ProtMut"]


#CLIENT INSTRUCTIONS
ClientRunStr = """import design
if not 'MDTempDir' in globals():
  import tempfile, os
  if not design.TempDirBase == "" and os.path.isdir(design.TempDirBase):
    globals()['MDTempDir'] = tempfile.mkdtemp(dir = design.TempDirBase)
  else:
    globals()['MDTempDir'] = tempfile.mkdtemp()
  SockTempFiles.append(globals()['MDTempDir'])
  print 'Temporary MD folder is:', globals()['MDTempDir']
Seq, ProtMut, ProtBB, ContactList, MDStages, First = SockData
Score, Pdb = design.RunMD(Seq, ProtMut, ProtBB, MDTempDir, ContactList, MDStages, First)
SockResult = (Seq, Score, Pdb)
"""




#FUNCTION DEFINITIONS

def GetMutationProbs(NRes):
  #look for a profile file
  if os.path.isfile("resprofile.txt"):
    f = file("resprofile.txt", "r")
    AAs = f.readline().split()
    AAs = [sequence.AAInst(x).AA3 for x in AAs]
    Ind = [i for i in range(len(AAs)) if not AAs[i] in SkipResMut]
    AAs = [AAs[i] for i in Ind]
    MutProbs = []
    while True:
      s = f.readline()
      if s == '': break
      #remove skipped residues and normalize
      Probs = array([float(x) for x in s.split()], float)
      Probs = Probs.take(Ind)
      Probs = Probs / sum(Probs)
      Probs = zip(AAs, Probs)
      MutProbs.append(Probs)
    f.close()
    return MutProbs
  else:  
    MutRes = [x for x in AllResMut if not x in SkipResMut]
    MutProb = 1. / len(MutRes)
    MutationList = [(x, MutProb) for x in MutRes]
    MutProbs = [MutationList] * NRes
  return MutProbs

def RandMutation(MutationList):
  "Returns a random mutation."
  Prob = random.random()
  CumProb = 0.
  for (Res, ResProb) in MutationList:
    CumProb += ResProb
    if Prob < CumProb:
      return Res

def RandResMutation(Seq):
  "Returns a sequence mutated by one residue."
  while True:
    ind = min(int(random.random() * len(Seq)), len(Seq) - 1)
    NewSeq = copy.deepcopy(Seq)
    NewSeq[ind] = RandMutation(MutProbs[ind])
    if not NewSeq[ind] == Seq[ind]:
      return NewSeq

def RandInitSeq(N):
  "Returns a random initial sequence."
  Seq = ['']*N
  for i in range(N):
    while Seq[i] == '' or Seq[i] in SkipResStart:
      Seq[i] = RandMutation(MutProbs[i])
  return Seq

def RandBest(Boltz = False):
  "Returns a best item, optionally with probability proportional to exp(-score/temp)."
  global BestList
  if Boltz:
    SeqProbs = array([b.Score for b in BestList], float) / MutTemp
    SeqProbs = exp(SeqProbs - SeqProbs.max())
    SeqProbs = SeqProbs / SeqProbs.sum()
    Prob = random.random()
    CumProb = 0.
    for i in range(NBest):
      CumProb += SeqProbs[i]
      if Prob < CumProb:
        return BestList[i]
  else:
    ind = int(random.random() * NBest)
    return BestList[ind]

def RandHybrid(Best1, Best2):
  "Returns a hybrid sequence and the closest BestClass."
  N1, N2 = 0, 0
  NewSeq = copy.deepcopy(Best1.Seq)
  for i in range(len(NewSeq)):
    if random.random() > 0.5:
      NewSeq[i] = Best1.Seq[i]
    else:
      NewSeq[i] = Best2.Seq[i]
    if NewSeq[i] == Best1.Seq[i]: N1 += 1
    if NewSeq[i] == Best2.Seq[i]: N2 += 1
  if N1 > N2:
    return NewSeq, Best1
  else:
    return NewSeq, Best2

def GetMutations(Seq, ProtMut):
  "Returns which residues were mutated."
  return [i for i in range(len(Seq)) if not Seq[i] == ProtMut.Seq[i]]

def GetRearrangeRes(Seq, ProtMut, Sys, ContactList):
  "Returns a list of mutated residues and their neighbors and a list of others."
  MutRes = GetMutations(Seq, ProtMut)
  BBRes = copy.deepcopy(MutRes)
  #add neighboring residues
  for i in MutRes:
    for j in ContactList[i]:
      if not j in BBRes:
        BBRes.append(j)
  #sort and make all list
  BBRes.sort()
  AARes = [i for i in range(len(Seq)) if not i in BBRes]
  return AARes, BBRes

def CopyFiles(MDTempDir):
  "Copies files back to the work directory."
  WorkDir = MDTempDir.strip()
  if WorkDir[-1] in ["\\","/"]: WorkDir = WorkDir[:-1]
  WorkDir = os.path.join("work/", os.path.basename(WorkDir))
  if not os.path.isdir(WorkDir): os.mkdir(WorkDir)
  for fn in glob.glob(os.path.join(MDTempDir, "*.pdb")):
    shutil.copy(fn, os.path.join(WorkDir, os.path.basename(fn)))

def RunMD(Seq, ProtMut, ProtBB, MDTempDir, ContactList,
          MDStages, First):
  "Runs the MD protocol for a sequence."
  s = mdsim.SimClass(MDTempDir)
  #make the backbone template (dehydrogened) and save
  Pdb = ProtBB.MutateSeq(Seq).GetPdb()
  Pdb = pdbtools.Dehydrogen(Pdb)
  PdbFn = os.path.join(MDTempDir, "target.pdb")
  file(PdbFn, "w").write(Pdb)
  s.SysInitPdbUsingSeq(PdbFn, Seq)
  s.SysBuild()
  s.SaveConfig(os.path.join(MDTempDir, "target.crd"))
  #mutate the target
  t = ProtMut.Copy()
  Muts = GetMutations(Seq, ProtMut)
  for i in Muts:
    t = t.MutateRes(i, Seq[i])
  #get the atoms we want to hold
  if First:
    AARes, BBRes = [], range(len(Seq))
  else:
    AARes, BBRes = GetRearrangeRes(Seq, ProtMut, s, ContactList)
  #now optimize the side chains for non-overlap and use as an initial config
  t.OptimizeSC(ResInd = BBRes, Ti = 1.e5, Tf = 1.e-5)
  Pdb = t.GetPdb()
  Pdb = pdbtools.Dehydrogen(Pdb)
  PdbFn = os.path.join(MDTempDir, "init.pdb")
  file(PdbFn, "w").write(Pdb)
  s.SysInitPdbUsingSeq(PdbFn, Seq)
  s.SysBuild()
  #set the temperature frequency
  s["STEPSTEMP"] = 100
  #turn on berendsen to remove energy generated
  s["TEMPMODE"] = 1
  #turn off center of mass resetting
  s.RecenterMD = False
  #add salt repulsion
  s.RestAddIonRepulsion()
  #set the restraints
  s.PosRestRefFile(os.path.join(MDTempDir, "target.crd"))
  #run stages
  Ct = 0
  for (NStepsMin1, NStepsMin2, NStepsMD, Stepsize, Temp1, Temp2, FConst,
       NBScale1, NBScale2, RadScale1, RadScale2, VelLimit) in MDStages:
    s["STEPSIZE"] = Stepsize
    if Stepsize < 0.002:
      s["SHAKEMODE"] = 1
    else:
      s["SHAKEMODE"] = 2
    s["TEMPSET1"], s["TEMPSET2"] = Temp1, Temp2
    s["NONBONDSCALE1"], s["NONBONDSCALE2"] = NBScale1, NBScale2
    s["RADIUSSCALE1"], s["RADIUSSCALE2"] = RadScale1, RadScale2
    s["VELOCITYLIMIT"] = VelLimit
    s.PosRestSetRes(AAResList = [], BBResList = "*", FConst = FConst)
    if NStepsMin1 > 0 or NStepsMin2 > 0:
      try:
        s.RunMin(NStepsMin1, NStepsMin2)
      except mdsim.SimError:
        CopyFiles(MDTempDir)
        return 1.e201, ""
      s.WritePdb(os.path.join(s.RunPath, "stage%dmin.pdb" % Ct))
    if NStepsMD > 0:
      try:
        s.RunMD(NSteps = NStepsMD)
      except mdsim.SimError:
        CopyFiles(MDTempDir)
        return 1.e201, ""
      s.WritePdb(os.path.join(s.RunPath, "stage%dmd.pdb" % Ct))
    Ct += 1
    EAvg = s["EPOTAVG"] 
    LastFConst = FConst
  #minimize and get the current pdb
  s.RunMin()
  Pdb = os.path.join(MDTempDir, "current.pdb")
  Pdb = file(Pdb, "r").read()
  CopyFiles(MDTempDir)
  return EAvg, Pdb
  
def AccMut(NewE, OldE):
  global MutTemp
  return exp(min(0, (OldE - NewE) / MutTemp)) > random.random()

def WriteBestSummary():
  global BestList
  s = ""
  for i in range(NBest):
    b = BestList[i]
    s += "%-3d %-10.2f %s\n" % (i+1, b.Score, sequence.SeqToAA1(b.Seq))
  file("bestseqs.txt", "w").write(s)

def NDiff(Seq1, Seq2):
  return [Seq1[i] == Seq2[i] for i in range(len(Seq1))].count(False)

def UpdateSeq(Seq, Score, Pdb, Index):
  "Updates the best list and output files."
  global BestList, SeqTol
  #check for an error
  if len(Pdb) == 0:
    s = "%-10d e %-10.2f %s\n" (Index, Score, sequence.SeqToAA1(Seq))
    file("progress.txt", "a").write(s)
  else:
    #first see if there is a low-homology sequence
    SortList = [(NDiff(Seq, BestList[i].Seq), i) for i in range(len(BestList))]
    SortList.sort()
    n, ind = SortList[0]
    if n > int(SeqTol * len(Seq)): ind = argmax([b.Score for b in BestList])
    #now check if it has a better score
    Acc = AccMut(Score, BestList[ind].Score) or n == 0
    s = "%-10d" % Index
    if Acc:
      BestList[ind].Score = Score
      BestList[ind].Seq = Seq
      BestList[ind].Prot = protein.ProteinClass(Pdb = Pdb)
      fn = "best%0*d.pdb" % (int(math.log(NBest, 10) + 1), ind + 1)
      file(fn, "w").write(Pdb)
      WriteBestSummary()
      s += " a"
    else:
      s += " r"
    #write to progress file
    s += " %-10.2f %s\n" % (Score, sequence.SeqToAA1(Seq)) 
    file("progress.txt", "a").write(s)


#ZSCORE STUFF

def AvgPCont(CO):
  "Returns the average contact probability as a function of contact order."
  if CO < 3:
    return 1.
  else:
    return exp(-3.52248 + 584.59511 * (CO - 3.78480) / (CO**4.10343 - 249.61419))
  
def AvgNNatCont(NRes):
  "Returns the average number of CO >= 3 native contacts for given length."
  return NRes * (4.2607 * NRes - 2.6029) / (NRes + 41.33)

def PrepZScore(NRes):
  "Preps the z-score"
  global ContEne, ContFreq
  ContEne = cPickle.load(file(CEFile, "r"))
  ContFreq = zeros((NRes, NRes), float)
  n = 0.
  for i in range(NRes):
    for j in range(i+1,NRes):
      ContFreq[i,j] = 1. #AvgPCont(j - i)
      ContFreq[j,i] = ContFreq[i,j]
      if j > i+2: n += ContFreq[i,j]
  n = AvgNNatCont(NRes) / n
  for i in range(NRes):
    for j in range(i+3,NRes):
      ContFreq[i,j] = ContFreq[i,j] * n
      ContFreq[j,i] = ContFreq[j,i] * n

def GetEne(a, b):
  if (a,b) in ContEne:
    return ContEne[(a,b)]
  elif (b,a) in ContEne:
    return ContEne[(b,a)]
  else:
    print sorted(ContEne.keys())
    raise IndexError, "Can't find " + repr((a,b))  
      
def ZScore(E, Seq, Pdb):
  "Returns the mean-field zscore for a sequence and energy."
  global ContEne, ContFreq
  EAvg, EVar = 0., 0.
  #self energies
  for i in range(NRes):
    U1, U2 = GetEne(Seq[i], Seq[i])
    EAvg += U2 / 2.
  #contact energies
  if False:
    for i in range(NRes):
      for j in range(i+1,NRes):
        U1, U2 = GetEne(Seq[i], Seq[j])
        U = U1 - U2
        f = ContFreq[i,j]
        EAvg += U * f
        EVar += U*U * f*(1. - f)
    return (E - EAvg) / sqrt(EVar)
  elif False:
    p = protein.ProteinClass(Pdb=Pdb)
    for (i,j) in p.ResContactList(MinCO=0):
      U1, U2 = GetEne(Seq[i], Seq[j])
      EAvg += U1 - U2
    return (E - EAvg) / len(p.Atoms)
  else:
    p = protein.ProteinClass(Pdb=Pdb)
    NPhobic = [sequence.Hydrophobic(x) for x in Seq].count(True)
    NCharge = [sequence.Charged(x) for x in Seq].count(True)
    return E * float(len(Seq)) / float(len(p.Atoms)) + 0.25*NPhobic + 5.0*NCharge   



#RESTART STUFF

def WriteRestart():
  dat = dict([(k,v) for (k,v) in globals().iteritems() if k in RestartData])
  cPickle.dump(dat, file("restart.dat","w"))

def LoadRestart():
  if os.path.isfile("restart.dat"):
    print "Loading restart data."
    dat = cPickle.load(file("restart.dat","r"))
    globals().update(dat)


if __name__ == '__main__':

  #MAIN LOOP

  #initialize the socketring class, which will decide who is the main server
  sr = socketring.RingClass()

  if sr.IsServer:

    MutStart = 0
    MutDone = 0

    #make the work directory
    if not os.path.isdir("work/"): os.mkdir("work")

    #Get the target structure
    Prot = protein.ProteinClass(Pdb = TargetPdb)
    Prot.ResAtom = 'CB'
    NRes = len(Prot)
    #make a backbone-only version
    ProtBB = Prot.Copy()

    #initialize the z-score
    PrepZScore(NRes)
        
    #make a contact list
    ContactMap = Prot.ResContactMap(MinCO = 0)
    ContactList = []
    for i in range(0, NRes):
      ContactList.append(flatnonzero(ContactMap[i]))

    #make the mutation probabilities
    MutProbs = GetMutationProbs(NRes)

    #best sequences and scores
    #initialize with random sequences if desired
    if RandomStartSeq:
      BestList = [BestClass(Seq = RandInitSeq(NRes), Prot = Prot) for i in range(NBest)]
    else:
      BestList = [BestClass(Seq = Prot.Seq, Prot = Prot) for i in range(NBest)]

    file("progress.txt", "w").write("Started simulation on " + time.strftime("%c") + "\n")

    #load any restart data
    LoadRestart()

    while MutDone < NMut:
      
      #first see if anything is done
      DoneTasks = [t for t in sr.Tasks if t.Done]
      for Task in DoneTasks:
        NewSeq, EAvg, NewPdb = Task.Result
        NewScore = ZScore(EAvg, NewSeq, NewPdb)
        sr.RemoveTask(Task)
        MutDone += 1
        UpdateSeq(NewSeq, NewScore, NewPdb, MutDone)
      if len(DoneTasks) > 0: WriteRestart()
        
      #now send out new mutations
      NRemain = NMut - MutDone
      NSend = min(len(sr.FreeClients()), NRemain)
      for i in range(NSend):
        #make a mutation; if this is the first one do the orig seq
        r = random.random()
        First = False
        if MutStart < NBest:
          First = True
          BestMut = BestList[MutStart]
          NewSeq = copy.deepcopy(BestMut.Seq)
        elif r < RerunFreq:
          #rerun an old mutation
          BestMut = RandBest(Boltz = False)
          NewSeq = copy.deepcopy(BestMut.Seq)
        elif r < HybridBestFreq + RerunFreq:
          #hybridize two best sequences
          Best1 = RandBest(Boltz = False)
          Best2 = RandBest(Boltz = False)
          while Best1 == Best2:
            Best2 = RandBest(Boltz = False)
          NewSeq, BestMut = RandHybrid(Best1, Best2)
        else:
          #make a new mutation from one of the best sequences
          BestMut = RandBest(Boltz = False)
          NewSeq = RandResMutation(BestMut.Seq)
        #send the task
        Data = (NewSeq, BestMut.Prot, ProtBB, ContactList, MDStages, First)
        sr.AddTask(ClientRunStr, Data)
        MutStart += 1

      #run some tasks
      sr.RunTasks(0.1)      

    sr.Finish()

  else:

    #run client and stop if something unexpected happened
    if not sr.RunClient():
      print "There were errors connecting with the socketring server."
      print "Terminating job."
      sys.exit()

       