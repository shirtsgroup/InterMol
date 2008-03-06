#!/usr/bin/env python

#LAST MODIFIED: 05-17-07

from numpy import *
import os, math, cPickle, time, sys, glob
import rexsock, pdbtools, sequence, makeconf, protein, mdsim, geometry
import zamdata

#global replica exchange run variables
NRepIncr = 5       #will only have integer multiples of this many replicas;
                   #a negative value uses the current number of processors
                   #with a minimum of abs(NRepIncr)
MinTemp = 270.0    #minimum replica temperature 
MaxTemp = 600.0    #maximum replica temperature
SwapRest = True    #swap restraints when they are sorted?

#cycle duration; first number is fragment size, second is cycle time in ps
#for all fragments of that size or less. fragments must be in order
#smallest to largest.  
CycleTimes = [(8, 20.), (12, 10.), (16, 5.), (24, 2.), (30, 1.)]

#rigidifying restraint params
RestRigidFConst = 0.1

#secondary structure restraint params
#these are (Phi, Psi, PhiTol, PsiTol)
RestSSFConst = 0.5
RestSSParams = {"H":(-58., -47., 10., 10.),
                "E":(-130., 125., 40., 40.)}

#contact restraint params; the distances here are absolute
RestContactParams = {"DIST1":0., "DIST2":0., "DIST3":7.5, "DIST4":8.0,
    "FCONST2":0.0, "FCONST3":0.5}
RestContAtom = "*"  #use residue centroid

#positional restraints for anchoring things to space
AnchorRestFConst = 1.0

#ANCHORING RESTRAINTS FOR MULTIPLE CHAINS:
#this is the distance added to Rsurf1+Rsurf2 at which a restraint starts to
#hold chains near each other so they don't fly off into space
ChainHoldDist = 10.
ChainHoldFConst = 1.0

#UMBRELLA SAMPLING RESTRAINTS
#umbrella sampling force constant
UmbrFConst = 0.5
#width of umbrella bins
UmbrDelta = 5.0
#umbrella atom
UmbrAtom = "*"  #centroid

#temperature optimization
TempOptTimes = []    #times in ps at which temperatures are optimized

#minimum fractional restraint strength when scaling
ScaleRestMin = 0.


def CalcNTemp(Frag):
  "Calculates the number of temperatures for a fragment."
  #this is an empirical relation based on past experience
  NRes = len(Frag.Seq)
  NTemp = log(MaxTemp / MinTemp) * 8.14018 * NRes**0.4
  #round up, and make sure greater than 0
  NTemp = max(int(NTemp + 0.5), 1)
  return NTemp

def CalcNRep(Frag):
  "Calculates the number of replicas to use for a fragment."
  #First calculate the ideal num of processors
  N = CalcNTemp(Frag)
  #now find the nearest multiple of NRepIncr
  N = NRepIncr * int(N / float(NRepIncr) + 0.5)
  N = max(N, abs(NRepIncr))
  if len(Frag.UmbrRest) > 0:
    N = 2*N + len(Frag.UmbrDist) - 2
  return N

def CalcTempScale(Frag):
  "Returns an array of temperatures for each replica."
  if Frag.NRep == 1:
    Temp = [MinTemp]
  else:
    if len(Frag.UmbrRest) > 0:
      NDist = len(Frag.UmbrDist)
      NTemp = (Frag.NRep - (NDist - 2)) / 2
    else:
      NTemp = Frag.NRep
    dlnT = (math.log(MaxTemp) - math.log(MinTemp)) / float(NTemp-1)
    Temp = [MinTemp*math.exp(float(i)*dlnT) for i in range(NTemp)]    
  return Temp

def CalcStates(Frag):
  "Calculates temperature and chain states."
  if len(Frag.UmbrRest) > 0:
    NDist = len(Frag.UmbrDist)
    NTemp = (Frag.NRep - (NDist - 2)) / 2
    States = [(i,0) for i in range(NTemp)] + \
             [(i,NDist-1) for i in range(NTemp)] + \
             [(0,i) for i in range(1,NDist-1)]
    States.sort()
  else:
    States = [(i,0) for i in range(Frag.NRep)]
  return States

def CalcSwaps(States):
  "Calculates allowed swaps."
  Swaps = []
  #find neighboring swaps
  for (i, (TempInd1, ChainInd1)) in enumerate(States):
    for (j, (TempInd2, ChainInd2)) in enumerate(States[i+1:]):
      if abs(TempInd1 - TempInd2) <= 1 and abs(ChainInd1 - ChainInd2) <= 1:
        Swaps.append([i, j+i+1])
  return Swaps

def CalcCycleTime(rx, Frag):
  "Calculates the number of ps between exchanges."
  n = len(Frag.Seq)
  for (FragLen, CycleTime) in CycleTimes:
    if n <= FragLen: return CycleTime
  #use the last one
  (FragLen, CycleTime) = CycleTimes[-1]
  return CycleTime

def FormatBuildString(s, Frag):
  """Replaces residue numbers with numbers for this fragment.
  The build string can have tokens like RES001, where numbers are
  relative to the full chain (starting at residue 1)."""
  t = s
  if Frag.CapN:
    Shift = 1 - Frag.StartRes
  else:
    Shift = -Frag.StartRes
  #fill in fragment-relative residue numbers
  while "RES" in t:
    n = t.index("RES") + 3
    ResTag, ResNum = "RES", ""
    for tn in t[n:]:
      if not tn in "0123456789": break
      ResNum += tn
      ResTag += tn
    ResNum = int(ResNum)
    if ResNum < Frag.StartRes + 1 or ResNum > Frag.StopRes + 1:
      t = t.replace(ResTag, "DELETELINE", 1)
    else:
      t = t.replace(ResTag, str(ResNum + Shift), 1)
  #remove non-pertinent lines
  t = t.split("\n")
  u = []
  for l in t:
    if not "DELETELINE" in l:
      u.append(l)
  u = "\n".join(u)
  return u


def PrepRefPdbs(Frag):
  "Sets reference file for rigidifying restraints in fragment."
  #check to see if rigid segments have already been proposed
  if len(Frag.ResRigid + Frag.ResAnchorBB + Frag.ResAnchorAA) > 0:
    #first see if any references have been specified
    RefPdbs = glob.glob(os.path.join(Frag.BasePath, "conf/ref*.pdb"))
    if len(RefPdbs) == 0:
      #check the seeds
      RefPdbs = glob.glob(os.path.join(Frag.BasePath, "conf/seed*.pdb"))
      if len(RefPdbs) > 0:
        #just use the first seed
        RefPdbs = RefPdbs[:1]
      else:
        #raise an error... no reference found
        raise IOError, "No reference pdb files for fragment %s found." % Frag.Name()
    #now make the reference files pretty for amber
    for (i, Pdb) in enumerate(RefPdbs):
      NewPdb = os.path.join(Frag.BasePath, "conf/refinit%03d.pdb" % (i+1))
      mdsim.PrepPdb(Pdb, OutPdb = NewPdb, Seq = Frag.Seq,
                    CapN = Frag.CapN, CapC = Frag.CapC)
      RefPdbs[i] = NewPdb
    return RefPdbs
  else:
    return None


def PrepChains(rx, zd, Frag, RefPdbs):
  """Locates CA atom indices closest to the centroid of each chain,
  along with the Rg of each chain."""
  #make a protein class
  if RefPdbs is None:
    #just take the first configuration
    p = protein.ProteinClass(Pdb = rx.Sim[0].GetPdb())
  else:
    #take the first reference
    p = protein.ProteinClass(Pdb = RefPdbs[0])
  #find the atom and residue ranges for each chain
  NChain = len(p.Chains)
  AtomRange = [p.ChainAtomNums[i:i+2] for i in range(NChain)]
  ResRange = [p.ChainResNums[i:i+2] for i in range(NChain)]
  #get the centroid and rg for each chain
  Centroid = [p.Centroid(ChainNum = i) for i in range(NChain)]
  Rg = [p.RadiusOfGyration(ResInd = range(a,b)) for (a,b) in ResRange]
  Rsurf = [sqrt(5./3.)*x for x in Rg]
  #get the CA index nearest to the centroid of each chain
  CAInd = [p.AtomInd(ResNum = range(a,b), AtomName = "CA")
           for (a,b) in ResRange]
  ChainAtomInd = [Ind[argmin(sum((p.Pos[Ind] - Centroid[i])**2, axis=1))]
               for (i,Ind) in enumerate(CAInd)]
  #write a summary file
  s = "CHAIN NUMBER, RESIDUE RANGE, ATOM RANGE, CENTROID ATOM, RG, Rsurf\n"
  Res, Atoms = p.Seq, p.AtomNames()
  for i in range(NChain):
    r1, r2 = ResRange[i]
    a1, a2 = AtomRange[i]
    s += "%d %s%d-%s%d %s(%d)-%s(%d) %s(%d) %.2f %.2f\n" % \
         (i, Res[r1], r1+1, Res[r2-1], r2,
          Atoms[a1], a1+1, Atoms[a2-1], a2,
          Atoms[ChainAtomInd[i]], ChainAtomInd[i]+1, Rg[i], Rsurf[i])
  file(os.path.join(Frag.BasePath, "data/chains.txt"), "w").write(s)
  #set the variables
  rx["ChainAtomInd"], rx["ChainRg"], rx["ChainRsurf"], rx["ChainResNums"] = \
                      ChainAtomInd, Rg, Rsurf, p.ChainResNums


def PrepInitPdbs(rx, zd, Frag, RefPdbs):
  """Prepares initial pdbs for running."""
  InitPdbs = glob.glob(os.path.join(Frag.BasePath, "conf/init*.pdb"))
  #get any anchored res
  ResInd = Frag.ResAnchorAA + Frag.ResAnchorBB
  #remove duplicates and shift res
  a = -Frag.StartRes
  if Frag.CapN: a += 1
  ResInd = [x+a for (i,x) in enumerate(ResInd) if not x in ResInd[i+1:]]
  if not RefPdbs is None and len(ResInd) > 0:
    for (i, InitPdb) in enumerate(InitPdbs):
      #find the corresponding reference pdb
      RefPdb = RefPdbs[i % len(RefPdbs)]
      #get positions of reference
      pref = protein.ProteinClass(Pdb = RefPdb)
      Ind1 = pref.AtomInd(AtomName = ["N", "CA", "C"], ResNum = ResInd)
      Pos1 = pref.Pos.take(Ind1, axis=0)
      #now align initial conf to reference
      p = protein.ProteinClass(Pdb = InitPdb)
      Ind2 = p.AtomInd(AtomName = ["N", "CA", "C"], ResNum = ResInd)
      Pos2 = p.Pos.take(Ind2, axis=0)
      Vec1, Vec2, RotMat, Resid = geometry.AlignmentRMSD(Pos1, Pos2, Center = True)
      p.Pos = dot(p.Pos + Vec2, RotMat) - Vec1
      p.WritePdb(InitPdb)
  #prep them for running using mdsim (dehydrogenating, etc for amber)
  for fn in InitPdbs:
    mdsim.PrepPdb(fn, CapN = Frag.CapN, CapC = Frag.CapC)
  #return the file list
  return InitPdbs

def PrepStates(rx, Frag):
  """Prepares temperature and umbrella sampling states."""
  States = CalcStates(Frag)
  Swaps = CalcSwaps(States)
  TempScale = CalcTempScale(Frag)
  rx.SetSwaps(Swaps)
  #set the states
  for (i, (TempInd, UmbrInd)) in enumerate(States):
    rx.Sim[i]["TempInd"] = TempInd
    rx.Sim[i]["UmbrInd"] = UmbrInd
  #set temperatures
  for s in rx.Sim:
    s["TEMPSET"] = TempScale[s["TempInd"]]  
  #write out the state indices
  fn = os.path.join(rx.DataPath, "states.txt")
  s = "REPLICA STATES\nreplica index, temperature index, umbrella index\n"
  s += "\n".join(["%d %d %d" % (i, t, c) for (i, (t,c)) in enumerate(States)])
  file(fn, "w").write(s)
  
    
def PrepFrag(rx, zd, Frag):
  "Prepares the REX simulation class for a particular fragment."
  #initialize the rx object; timings and restraints can change
  #for each run, so reset them each time this is called, even
  #if from restart file
  if Frag.NRep == 0: Frag.NRep = CalcNRep(Frag)
  rx.InitServer(BasePath = Frag.BasePath, N = Frag.NRep, Restart = True)

  #assign temperature and umbrella sampling states to fragment
  PrepStates(rx, Frag)

  #check to set default values  
  if not "ChainAtomInd" in rx: rx["ChainAtomInd"] = []
  if not "ChainResNums" in rx: rx["ChainResNums"] = []
  if not "ChainRg" in rx: rx["ChainRg"] = []
  if not "ChainRsurf" in rx: rx["ChainRsurf"] = []

  #setup
  rx["STAGE"] = "START"
  if not rx.Restart:
    #update the progress file
    ProgFile = os.path.join(Frag.BasePath, "progress.txt")
    file(ProgFile,"w").write("starting up...\n")

    #load in system setup variables
    PreBuildString = FormatBuildString(zd.PreBuildString, Frag)
    PostBuildString = FormatBuildString(zd.PostBuildString, Frag)
    rx.SimIter(lambda s: s.SetPreBuildString(PreBuildString))
    rx.SimIter(lambda s: s.SetPostBuildString(PostBuildString))

    #check for rigid segment reference file
    RefPdbs = PrepRefPdbs(Frag)
    #make the reference point
    if not RefPdbs is None:
      for (i, s) in enumerate(rx.Sim):
        #get this reference pdb file
        Pdb = RefPdbs[i % len(RefPdbs)]
        s.SysInitPdbUsingSeq(Pdb, Frag.Seq, CapN = Frag.CapN, CapC = Frag.CapC)
        s.SysBuild()
        s.PosRestRefCurrent()
        #set the reference pdb
        s["RefPdb"] = Pdb
        
    #check for initial configs
    InitPdbs = PrepInitPdbs(rx, zd, Frag, RefPdbs)
    #if there are initial configs, spread over all replicas
    if len(InitPdbs) > 0:
      for (i, s) in enumerate(rx.Sim):
        Pdb = InitPdbs[i % len(InitPdbs)]
        s.SysInitPdbUsingSeq(Pdb, Frag.Seq, CapN = Frag.CapN, CapC = Frag.CapC)
      #build the systems
      rx.SimIter(lambda s: s.SysBuild())
    #otherwise start with extended, provided that there are no reference files
    #(if there are reference files, this will just used the systems already built)
    elif RefPdbs is None:
      for s in rx.Sim:
        s.SysInitSeq(Frag.Seq, CapN = Frag.CapN, CapC = Frag.CapC)
      #build the systems
      rx.SimIter(lambda s: s.SysBuild())
    
    #get anchoring indices
    PrepChains(rx, zd, Frag, RefPdbs)
    
    #add restraints
    for s in rx.Sim:
      AddRest(s, Frag, rx["ChainAtomInd"], rx["ChainResNums"], rx["ChainRsurf"])

    #save the rex and sim data
    rx.RestartWrite()

  #set the exchange timing
  if Frag.CycleTime == 0:
    Frag.CycleTime = CalcCycleTime(rx, Frag)
  #round to nearest multiple
  NStep = int(Frag.CycleTime / rx.Sim[0]["STEPSIZE"] + 0.5)
  Frag.CycleTime = float(NStep) * rx.Sim[0]["STEPSIZE"]
  #set the steps variable
  rx.SimSetVar("STEPSMD", [NStep]*rx.NReplica)

  #turn off recentering if there are cartesian restraints
  AutoRecenter = True
  if len(Frag.ResAnchorAA + Frag.ResAnchorBB) > 0: AutoRecenter = False
  for s in rx.Sim:
    s.AutoRecenter = AutoRecenter  


def AddPhiPsiRest(Sim, ResNum, IndN, IndCA, IndC, Phi, Psi, PhiTol, PsiTol, FConst):
  """Adds a phi-psi restraint to a simulation class."""
  Res = Sim.Seq[ResNum]
  if not ResNum == 0 and not Res in protein.FixedPhi:
    a1, a2, a3, a4 = IndC[ResNum-1], IndN[ResNum], IndCA[ResNum], IndC[ResNum]
    if a1 > 0 and a2 > 0 and a3 > 0 and a4 > 0:
      Sim.RestSetAtoms(a1, a2, a3, a4,
                       FConst2 = FConst, FConst3 = FConst,
                       Dist1 = Phi - 179., Dist2 = Phi - PhiTol,
                       Dist3 = Phi + PhiTol, Dist4 = Phi + 179.)
  if not ResNum == len(Sim.Seq) - 1 and not Res in protein.FixedPsi:
    a1, a2, a3, a4 = IndN[ResNum], IndCA[ResNum], IndC[ResNum], IndN[ResNum+1]
    if a1 > 0 and a2 > 0 and a3 > 0 and a4 > 0:
      Sim.RestSetAtoms(a1, a2, a3, a4,
                       FConst2 = FConst, FConst3 = FConst,
                       Dist1 = Psi - 179., Dist2 = Psi - PsiTol,
                       Dist3 = Psi + PsiTol, Dist4 = Psi + 179.)


def AddRest(Sim, Frag, ChainAtomInd, ChainResNums, Radius, RestList = None):
  "Adds both contact and rigidifying restraints to a Sim."
  Sim.RestClear()
  if Frag.CapN:
    ResShift = 1
  else:
    ResShift = 0
  #calculate offset
  Offset = ResShift - Frag.StartRes

  #GET SCALING
  RestScale = 1.0
  if "RestScale" in Sim: RestScale = Sim["RestScale"]
  
  #FIRST GET ATOM INDICES
  IndN = [-1 for i in Sim.Seq]
  IndCA = [-1 for i in Sim.Seq]
  IndC = [-1 for i in Sim.Seq]
  for (i,a) in enumerate(Sim.Atoms):
    a, rn = a.strip(), Sim.AtomRes[i]
    if a == "N":
      IndN[rn] = i
    elif a == "CA":
      IndCA[rn] = i
    elif a == "C":
      IndC[rn] = i

  #KEEP TRACK OF WHICH CHAINS HAVE BEEN CONNECTED
  def GetChainNum(ResNum):
    for (i, a) in enumerate(ChainResNums[1:]):
      if ResNum < a: return i
    return 0
  #first make a maxrix of chains
  NChain = len(ChainAtomInd)
  ChainsLinked = zeros((NChain,NChain), bool)

  #ADD RIGIDIFYING RESTRAINTS
  if "RefPdb" in Sim:
    RefPdb = Sim["RefPdb"]
    if os.path.isfile(RefPdb):
      #read pdb and phi psi values
      p = protein.ProteinClass(RefPdb)
      PhiPsi = [p.PhiPsi(i) for i in range(len(p))]
      #add the restraints
      for ResNum in Frag.ResRigid:
        rn = ResNum + Offset
        Phi, Psi = PhiPsi[rn]
        AddPhiPsiRest(Sim, rn, IndN, IndCA, IndC, Phi, Psi, 0., 0.,
                      RestScale * RestRigidFConst)

  #ADD SECONDARY STRUCTURE RESTRAINTS
  for (i, s) in enumerate(Frag.SS):
    if s in RestSSParams:
      Phi, Psi, PhiTol, PsiTol = RestSSParams[s]
      rn = i + ResShift
      AddPhiPsiRest(Sim, rn, IndN, IndCA, IndC, Phi, Psi, PhiTol, PsiTol,
                    RestScale * RestSSFConst)
          
  #ADD SALT BRIDGE REPULSIONS
  if Frag.RepelSalt: Sim.RestAddIonRepulsion()
  
  #ADD CONTACT RESTRAINTS
  #change parameters
  for (k,v) in RestContactParams.iteritems():
    if k.startswith("FCONST"):
      Sim[k] = v * RestScale
    else:
      Sim[k] = v
  #check to see if specific restraints were sent
  if RestList is None: RestList = Frag.PriorRest
  for (a,b) in RestList:
    if Frag.ContainsSeg(a,b):
      Sim.RestSetRes(a + Offset, b + Offset, RestContAtom, RestContAtom,
                     Strength = RestScale)
      
  #ADD POSITIONAL / ANCHORING RESTRAINTS
  if len(Frag.ResAnchorBB + Frag.ResAnchorAA) > 0:
    AAResList = [a + Offset for a in Frag.ResAnchorBB]
    BBResList = [a + Offset for a in Frag.ResAnchorAA]
    Sim.PosRestSetRes(AAResList = AAResList, BBResList = BBResList,
                      FConst = AnchorRestFConst, Strength = RestScale)

  #ADD UMBRELLA SAMPLING RESTRAINTS
  UmbrInd = Sim["UmbrInd"]
  for (i, (a,b)) in enumerate(Frag.UmbrRest):
    r1 = max(0., Frag.UmbrDist[UmbrInd] - UmbrDelta * 0.5)
    if UmbrInd == 0: r1 = 0.
    r2 = Frag.UmbrDist[UmbrInd] + UmbrDelta * 0.5
    Sim.RestSetRes(a + Offset, b + Offset, UmbrAtom, UmbrAtom,
                   FConst2 = UmbrFConst, FConst3 = UmbrFConst,
                   Dist1 = 0., Dist2 = r1, Dist3 = r2, Dist4 = 10. * r2)
    #update which chains we added it to
    aChain, bChain = GetChainNum(a), GetChainNum(b)
    ChainsLinked[aChain,bChain] = True
    ChainsLinked[bChain,aChain] = True
    
  #ADD CHAIN RESTRAINTS
  #add restraints between anchoring atoms so that chains don't fly off,
  #only to chains that don't have umbrella sampling
  for (i, a1) in enumerate(ChainAtomInd):
    for (j, a2) in enumerate(ChainAtomInd[i+1:]):
      if ChainsLinked[i,j+i+1]: continue
      r = Radius[i] + Radius[j+i+1] + ChainHoldDist
      Sim.RestSetAtoms(a1, a2, FConst2 = 0., FConst3 = ChainHoldFConst,
                       Dist1 = 0., Dist2 = 0., Dist3 = r, Dist4 = 10. * r)


def SetupRex(rx, zd, Frag, PerformSwaps = True,
  RestListAll = None, RestListSort = None):
  "Sets up a REX simulation for speficied number of picoseconds."

  ProgFile = os.path.join(Frag.BasePath, "progress.txt")

  #update fragment variables
  (StepSize, NStep, StepsSave) = (rx.Sim[0]["STEPSIZE"],
    rx.Sim[0]["STEPSMD"], rx.Sim[0]["STEPSSAVE"])
  Frag.FrameTime = StepSize * StepsSave
  Frag.NCycle = int(Frag.RexTime / (StepSize * NStep) + 0.5)
  Frag.NFrame = Frag.NCycle * NStep / StepsSave

  #setup restraint strengths
  if Frag.ScaleRest:
    n = float(len(rx.Sim) - 1)
    for (i,s) in enumerate(rx.Sim):
      s["RestScale"] = ScaleRestMin + (1. - ScaleRestMin) * float(i) / n
  else:
    for s in rx.Sim:
      s["RestScale"] = 1.

  #update restraints; some restraints will be applied to all replicas;
  #some will be sorted across the replicas
  if RestListAll is None:
    RestListAll = Frag.PriorRest
  if RestListSort is None or RestListSort == []:
    for s in rx.Sim:
      AddRest(s, Frag, rx["ChainAtomInd"], rx["ChainResNums"], 
              rx["ChainRsurf"], RestListAll)
  else:
    for (i,s) in enumerate(rx.Sim):
      #assign restraints based on replica number (so we can recover
      #the same swapped restraints in each replica if we restart)
      RestNum = rx.Replica[i] % len(RestListSort)
      AddRest(s, Frag, rx["ChainAtomInd"], rx["ChainResNums"],
              rx["ChainRsurf"], RestListAll + [RestListSort[RestNum]])

  #cycle loop
  rx["STAGE"] = "CYCLES"
  if len(Frag.UmbrRest) > 0:
    Labels = ["%.fK,U%d" % (s["TEMPSET"], s["UmbrInd"]) for s in rx.Sim]
  else:
    Labels = ["%.fK" % s["TEMPSET"] for s in rx.Sim]

  #update the temperature optimization cycles
  TempOptCycles = [int(float(x) / Frag.CycleTime) for x in TempOptTimes]
  rx.TempOptCycles = TempOptCycles

  #initialize the Rex
  if RestListSort is None or RestListSort == []:
    rx.StartRex(rexsock.PaccFuncRest, Frag.NCycle, Labels = Labels,
      ProgFile = ProgFile, PerformSwaps = PerformSwaps)
    rx.SwapRest = False
  else:
    rx.StartRex(rexsock.PaccFuncBasic, Frag.NCycle, Labels = Labels,
      ProgFile = ProgFile, PerformSwaps = PerformSwaps)
    rx.SwapRest = SwapRest


def RunRex(rx, zd, Frag):
  "Runs the next elements in the REX simulation."
  if rx.Minimized:
    rx["STAGE"] = "CYCLES"
    if rx.RunRex():
      rx["STAGE"] = "DONE"
      return True
    else:
      return False
  else:
    rx["STAGE"] = "MIN"
    rx.RunMin()


def RestartRex(rx, zd, Frag):
  "Erases counts and accumulated trajectory data but keeps latest configs."
  rx.ResetData()
  rx.ResetCount()

def RestartRexData(rx, zd, Frag):
  "Erases accumulated trajectory data but keeps latest configs."
  rx.ResetData() 

def CleanUp(rx, zd, Frag, ConserveDisk):
  "Removes workspace and data files from a REX simulation."
  if ConserveDisk:
    for DelPath in ["work/", "data/"]:
      Path = os.path.join(Frag.BasePath, DelPath)
      for f in glob.glob(Path):
        try:
          os.remove(f)
        except OSError:
          print "Could not remove file %s" % f
      try:
        os.rmdir(Path)
      except OSError:
        print "Could not remove path %s" % Path

    
