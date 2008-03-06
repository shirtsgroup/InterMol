#!/usr/bin/env python


#LAST MODIFIED: 09-01-06


Usage = """Assembles two pdb structures.

Usage     : assemble.py PDB1 PDB2 LOOPSEQ PREFIX [OPTIONS]

PDB*      : input pdb files
LOOPSEQ   : one-letter sequence codes
PREFIX    : prefix of output pdb files; will append 001, 002, etc
OPTIONS   : "--maxconf=N" to set maximum number of configs (dflt 50)
            "--debug" to keep temporary files and produce debug output
            "--overlap=R" to set overlap distance to R (dflt 1.0)
            "--mc" to use the Monte Carlo method
"""

#check for instructions
import sys
if __name__ == "__main__" and len(sys.argv) < 3:
  print Usage
  sys.exit()


from numpy import *
import os, sys, math, shutil, glob, rmsd
import sequence, pdbtools, protein, geometry
import scripttools

#GLOBALS:
OverlapDist = 1.0
DEBUG = False
SplatFiles = ["amino_acids.mol2", 'bbdep02.May.sortlib',
              'ramachandran.lib', "splat", "assemble"]
AlignFiles = ["align"]

#SPLAT GLOBALS
MaxSample = 4000
MaxPerturb = 15
NRotSample = 20

#MSS GLOBALS
StericWeight = 1.
RgAllWeight = 0.025
RgPhobicWeight = 0.1
LoopWeight = 1.
AnchorWeight = 1.
HBondWeight = 10.
MaxdAng1Dflt, MaxdAng2Dflt, MaxdRDflt = 20., 20., 1.0
MaxdAng1, MaxdAng2, MaxdR = MaxdAng1Dflt, MaxdAng2Dflt, MaxdRDflt
MaxdDihDflt = 10.
MaxdDih = MaxdDihDflt
AnnealTempInit, AnnealTempFinal = 10., 0.1
NStepAnnealDflt = 1000
AnnealSched1 = [(100., 100., 100)]
AnnealSched2 = [(100., 10., 400), (10., 1., 400), (1., 0.5, 100)]
NRot = 10
AlignBuffer = 8.
AnchorBuffer = 9.
DistPerRes = 3.3


#-------- CHECKING FOR FILES --------

Path = os.environ["PATH"]
if os.name == 'nt':
  Path = Path.split(';')
else:
  Path = Path.split(':')

#CHECK FOR SPLAT PATH
SplatPath = ""
if "SPLATPATH" in os.environ:
  SplatPath = os.environ["SPLATPATH"]
  x = [os.path.isfile(os.path.join(SplatPath, f)) for f in SplatFiles]
  if x.count(False) > 0:
    print "Error: could not find all splat files"
    print "Required files are " + " ".join(SplatFiles)
    #sys.exit(1)
else:
  #look in paths
  Found = False
  for p in ["./"] + Path:
    x = [os.path.isfile(os.path.join(p, f)) for f in SplatFiles]
    if x.count(False) == 0:
      Found = True
      SplatPath = p
      break
  if Found == False:
    print "Error: could not find all splat files"
    print "Required files are " + " ".join(SplatFiles)
    #sys.exit(1)
AssembleExec = os.path.join(SplatPath, "assemble")
SplatExec = os.path.join(SplatPath, "splat")
    
#CHECK FOR ALIGN PATH
AlignPath = ""
if "ALIGNPATH" in os.environ:
  AlignPath = os.environ["ALIGNPATH"]
  x = [os.path.isfile(os.path.join(AlignPath, f)) for f in AlignFiles]
  if x.count(False) > 0:
    print "Error: could not find all align files"
    print "Required files are " + " ".join(AlignFiles)
    #sys.exit(1)
else:
  #look in paths
  Found = False
  for p in ["./"] + Path:
    x = [os.path.isfile(os.path.join(p, f)) for f in AlignFiles]
    if x.count(False) == 0:
      Found = True
      AlignPath = p
      break
  if Found == False:
    print "Error: could not find all align files"
    print "Required files are " + " ".join(AlignFiles)
    #sys.exit(1)
AlignExec = os.path.join(AlignPath, "align")


#-------- MSS ROUTINES --------

def PrepScore(p, AnchorResInd):
  "Prepares a ProteinClass for the scoring function."
  p.PhobicResInd = p.ResInd(Hydrophobic = True)
  p.AnchorAtom = p.Res[AnchorResInd].AtomNum("CB") \
                 + p.Res[AnchorResInd].StartAtom

def Score(p1, p2, LoopLen):
  #get residue positions
  ResPos1, ResPos2 = p1.ResPos(), p2.ResPos()
  #get Rg of everything; cube so it's extensive in NRes
  RgAll = protein.RadiusOfGyration(p1, p2,
          ResPos1 = ResPos1, ResPos2 = ResPos2)
  RgAll = RgAll * RgAll * RgAll
  #get Rg of hydrophobics; cube so it's extensive in NRes
  RgPhobic = protein.RadiusOfGyration(p1, p2,
             ResInd1 = p1.PhobicResInd, ResInd2 = p2.PhobicResInd,
             ResPos1 = ResPos1, ResPos2 = ResPos2)
  RgPhobic = RgPhobic * RgPhobic * RgPhobic
  #get the Sterics
  StericScore = protein.StericScore(p1, p2)
  #get the hbond energy
  HBondScore = protein.HBondScore(p1, p2)
  #get the loop energy
  l = geometry.Length(p1.Pos[-1] - p2.Pos[0])
  LoopE = max(l - DistPerRes * LoopLen + 2, 0)**12
  #get the anchor energy
  l = geometry.Length(p1.Pos[p1.AnchorAtom] - p2.Pos[p2.AnchorAtom])
  AnchorE =  max(l - AnchorBuffer + 1, 0)**12
  return AnchorWeight * AnchorE \
         + LoopWeight * LoopE \
         + RgAllWeight * RgAll \
         + RgPhobicWeight * RgPhobic \
         + StericWeight * StericScore \
         + HBondWeight * HBondScore


def RunMCSteps(p1, p2, LoopLen, NStep, Temp1, Temp2,
               StepsUpdateMax = None,
               StepsUpdateTemp = None,
               Verbose = False):
  global MaxdAng1, MaxdAng2, MaxdR
  #initialize variables
  OldE = Score(p1, p2, LoopLen)
  p1Centroid = p1.Centroid()
  dTemp = (Temp2 - Temp1) / (NStep - 1)
  Temp = Temp1
  
  #run the steps
  NAcc, NAtt = zeros(3,float), ones(3,float) * 0.1
  for i in range(0, NStep):
    #random move
    AccFact = 0.
    Move = int(random.random()*3)
    OldPos = p2.Pos.copy()
    if Move == 0:
      #random translation towards or away from p1
      dR = (2*random.random() - 1) * MaxdR
      Vec = p2.Centroid() - p1.Centroid()
      p2.Translate(dR * geometry.UnitVec(Vec))
      OldRSq = dot(Vec, Vec)
      NewRSq = (sqrt(OldRSq) + dR)**2
      AccFact = log(NewRSq/OldRSq)
    elif Move == 1:
      #random rotation of p2 around p1
      Vec = geometry.RandVec()
      p2.Rotate(Vec, (2*random.random() - 1) * MaxdAng1, p1Centroid)
    else:
      #random rotation of p2 around its centroid
      p2.Rotate(geometry.RandVec(), (2*random.random() - 1) * MaxdAng2)
    NAtt[Move] += 1

    #accept/reject    
    NewE = Score(p1, p2, LoopLen)
    Acc = exp(min(0., -(NewE-OldE) / Temp) + AccFact) > random.random()
    if Acc:
      OldE = NewE
      NAcc[Move] += 1
    else:
      p2.Pos = OldPos
    if i % 20 == 0 and Verbose:
      print "round %d  %.2f %.2f %.2f %.2f %.2f %.2f %.2f" % \
       (i, 100*NAcc[0]/NAtt[0], 100*NAcc[1]/NAtt[1], 100*NAcc[2]/NAtt[2],
        Temp, MaxdR, MaxdAng1, MaxdAng2)
    if i % 20 == 0:
      Pdb = "../mc%05d.pdb" % i
      #p1.Copy().Concat(p2).WritePdb(Pdb)
      
    #update the maximum moves
    if not StepsUpdateMax is None and (i+1) % StepsUpdateMax == 0:
      MaxdAng1 = min(MaxdAng1 * (NAcc[1]/NAtt[1] + 0.1)/0.6, 180)
      MaxdAng2 = min(MaxdAng2 * (NAcc[2]/NAtt[2] + 0.1)/0.6, 180)
      MaxdR = MaxdR * (NAcc[0]/NAtt[0] + 0.1)/0.6
      NAcc, NAtt = zeros(3,float), ones(3,float) * 0.1
      
    #change temp
    Temp += dTemp


def AlignRMSD(p1, p2):
  Pos1 = p1.Pos.take(p1.AtomInd(AtomName = ["N", "CA", "C"]), axis=0)
  Pos2 = p2.Pos.take(p2.AtomInd(AtomName = ["N", "CA", "C"]), axis=0)
  return sqrt(((Pos1 - Pos2)**2).sum() / len(Pos1))


def Align(p1, p2, LoopLen, ResInd1, ResInd2, MaxDist, NStepAnneal):

  #align the two residues
  Pos1 = p1.ResPos(ResAtom = "*")
  Pos2 = p2.ResPos(ResAtom = "*")
  Vec1 = geometry.UnitVec(Pos1[ResInd1] - p1.Centroid())
  Vec2 = geometry.UnitVec(Pos2[ResInd2] - p2.Centroid())
  Vec, Ang = geometry.GetVecMapping(Vec1, Vec2)
  p2.Rotate(-Vec, 180 - Ang)
  
  #place target residues on top of each other
  Pos1 = p1.ResPos([ResInd1], "*")[0]
  Pos2 = p2.ResPos([ResInd2], "*")[0]
  p2.Translate(Pos1 - Pos2)
  
  #move so just a small buffer between
  Ind1 = dot(p1.Pos, Vec1).argmax()
  Ind2 = dot(p2.Pos, -Vec1).argmax()
  Buffer = (AlignBuffer - dot(p2.Pos[Ind2] - p1.Pos[Ind1], Vec1)) * Vec1
  p2.Translate(Buffer)

  #set the anchoring atoms and hydrophobic residues
  PrepScore(p1, ResInd1)
  PrepScore(p2, ResInd2)
  
  #sort through rotations and pick the best
  dAng = 360. / NRot
  E = []
  for i in range(0, NRot):
    p2.Rotate(Vec1, dAng)
    E.append(Score(p1, p2, LoopLen))
  p2.Rotate(Vec1, (argmax(E) - NRot - 1) * dAng)
  
  #optimize with monte carlo
  RunMCSteps(p1, p2, LoopLen, NStepAnneal, AnnealTempInit, AnnealTempInit)
  RunMCSteps(p1, p2, LoopLen, NStepAnneal, AnnealTempInit, AnnealTempFinal)
  RunMCSteps(p1, p2, LoopLen, NStepAnneal, AnnealTempFinal, AnnealTempFinal)

    
def AssembleMC(PdbFile1, PdbFile2, LoopSeq, Prefix, MaxConf,
              MaxConfAlign = None, NStepAnneal = NStepAnnealDflt,
              MinRMSD = 2.0,
              SaveAlign = False, Verbose = True, CleanUp = True):

  if MaxConfAlign is None:
    MaxConfAlign = max(MaxConf, 100)
  
  #change to a temp directory
  n = 0
  while os.path.isdir("pdbassem%d" % n):
    n += 1
  if n > 49:
    raise IOError, "Error: too many pdbassem temporary folders."
  Dir = "pdbassem%d" % n
  os.mkdir(Dir)
  os.chdir(Dir)

  #prep the pdb files
  PrepPdb(os.path.join("../", PdbFile1), "pep1.pdb")
  PrepPdb(os.path.join("../", PdbFile2), "pep2.pdb")
  
  #make protein classes
  p1 = protein.ProteinClass(Pdb = "pep1.pdb")
  p2 = protein.ProteinClass(Pdb = "pep2.pdb")
  p1.ResAtom, p2.ResAtom = "CB", "CB"
  
  #find the maximum radii
  MaxR1 = sqrt(((p1.Pos - p1.Centroid())**2).sum(axis=1).max())
  MaxR2 = sqrt(((p2.Pos - p2.Centroid())**2).sum(axis=1).max())
  MaxDist = MaxR1 + MaxR2
  
  #get the loop
  LoopSeq = sequence.SeqToList(LoopSeq)
  LoopLen = len(LoopSeq)
   
  #make a data holder
  class TempClass:
    def __init__(self, Score, Pdb, Pair):
      self.Score = Score
      self.Pdb = Pdb
      self.Pair = Pair
      self.Rank = -1
  Generated = []
  
  #sort through hydrophobic pairs on the exterior and align
  #if there are no hydrophobics, just choose every third residue
  p1ResInd = p1.ResInd(Hydrophobic = True, MaxCoord = protein.IntMinCoordDflt - 1)
  if len(p1ResInd) == 0: p1ResInd = range(0, len(p1), 3)
  p2ResInd = p2.ResInd(Hydrophobic = True, MaxCoord = protein.IntMinCoordDflt - 1)
  if len(p2ResInd) == 0: p2ResInd = range(0, len(p2), 3)
  NAlign = len(p1ResInd) * len(p2ResInd)
  NDone = 0
  for i in p1ResInd:
    for j in p2ResInd:
      NDone += 1
      if Verbose: print "Performing alignment %d/%d..." % (NDone, NAlign)
      p2new = p2.Copy()
      Align(p1, p2new, LoopLen, i, j, MaxDist, NStepAnneal)
      #compare to the scores
      s = Score(p1, p2new, LoopLen)
      if Verbose: print "...final score is %.3f" % s
      if len(Generated) < MaxConfAlign:
        Pdb = "pep2_%03d_temp.pdb" % (len(Generated) + 1)
        Generated.append(TempClass(s, Pdb, (i,j)))
        p2new.WritePdb(Pdb)
      else:
        Scores = array([g.Score for g in Generated], float)
        ind = Scores.argmax()
        if s < Scores[ind]:
          Pdb = "pep2_%03d_temp.pdb" % (ind+1)
          Generated[ind] = TempClass(s, Pdb, (i,j))
          p2new.WritePdb(Pdb)
          
  #sort and reorder files
  Generated.sort(cmp=lambda x,y: cmp(x.Score, y.Score))
  #check rmsd values
  i = 1
  while i < len(Generated):
    g = Generated[i]
    p = protein.ProteinClass(Pdb = g.Pdb)
    #see if the rmsd is too close to another of better score
    RMSDFlag = [AlignRMSD(p, protein.ProteinClass(Pdb = g2.Pdb)) > MinRMSD
                for g2 in Generated[:i]]
    if False in RMSDFlag:
      if Verbose: print "Removing alignment %d; too similar to another alignment." % i
      os.remove(g.Pdb)
      del Generated[i]
    else:
      i += 1
  #rename the files as sorted
  i = 0
  for g in Generated:
    NewPdb = "pep2_%03d.pdb" % (i+1)
    shutil.move(g.Pdb, NewPdb)
    g.Pdb = NewPdb
    g.Rank = i
    i += 1

  #see if we need to save the alignments
  if SaveAlign:
    for i in range(0, len(Generated[:MaxConf])):
      g = Generated[i]
      p2new = protein.ProteinClass(Pdb = g.Pdb)
      NewPdb = "%s_align_%03d.pdb" % (Prefix, i+1)
      p1.Copy().Concat(p2new).WritePdb(os.path.join("../", NewPdb))
    
  #make loops
  if Verbose: print "Making loops"
  MakeLoop("pep1.pdb", "pep2.pdb", LoopSeq, "done",
           NPdb = len(Generated), DoPdbPrep = False)
  
  #find loop files and match up with generated
  Done = glob.glob("done*.pdb")
  for g in Generated: g.Pdb = ""
  for d in Done:
    i = int(d.replace("_001.pdb","").split("_")[-1]) - 1
    Generated[i].Pdb = d
  Generated = [g for g in Generated if not g.Pdb == ""]

  #trim to the right number
  Generated = Generated[:MaxConf]
  
  #move files and rename
  for i in range(0, len(Generated)):
    NewPdb = "%s_%03d.pdb" % (Prefix, i+1)
    shutil.move(Generated[i].Pdb, os.path.join("../", NewPdb))
    Generated[i].Pdb = NewPdb
    
  #change back
  os.chdir("../")
  #write a summary file
  if Verbose: print "Writing summary file"
  s = "ALIGNMENT RESULTS\n"
  s += "Number, Alignment rank, Pdb, Aligned pair, Score\n"
  for i in range(0, len(Generated)):
    g = Generated[i]
    (ResInd1, ResInd2) = g.Pair
    Res1 = p1.Res[ResInd1].Name
    Res2 = p2.Res[ResInd2].Name
    ResInd2 += len(p1) + LoopLen
    PairName = "%s%d-%s%d" % (Res1, ResInd1+1, Res2, ResInd2+1)
    s += "%-3d %-3d %-20s %-15s %.2f\n" % (i+1, g.Rank+1, g.Pdb, PairName, g.Score)
  file("%s-alignments.txt" % Prefix, "w").write(s)
  
  #remove temp dir
  if CleanUp: shutil.rmtree(Dir)


#-------- MSS ROUTINES 2: no loop building ---------

def PrepScore2(p, LoopRes):
  "Prepares a ProteinClass for the scoring function."
  p.PhobicResInd = [x for x in p.ResInd(Hydrophobic = True)
                    if not x in LoopRes]
  p.PhobicResInd = array(p.PhobicResInd, int)
  p.RgResInd = [x for x in range(0, len(p))
                if not x in LoopRes]
  p.RgResInd = array(p.RgResInd, int)
  
def Score2(p):
  #get residue positions
  ResPos = p.ResPos()
  #get Rg of everything; cube so it's extensive in NRes
  RgAll = p.RadiusOfGyration(ResInd = p.RgResInd, ResPos = ResPos)
  RgAll = RgAll * RgAll * RgAll
  #get Rg of hydrophobics; cube so it's extensive in NRes
  RgPhobic = p.RadiusOfGyration(ResInd = p.PhobicResInd, ResPos = ResPos)
  RgPhobic = RgPhobic * RgPhobic * RgPhobic
  #get the Sterics
  StericScore = p.StericScore()
  #get the hbond energy
  HBondScore = p.HBondDipoleScore()
  return   RgAllWeight * RgAll \
         + RgPhobicWeight * RgPhobic \
         + StericWeight * StericScore \
         + HBondWeight * HBondScore


def RunMCSteps2(p, LoopRes, NStep, Temp1, Temp2,
               StepsUpdateMax = None,
               StepsUpdateTemp = None,
               Verbose = False):
  global MaxdDih
  #initialize variables
  OldE = Score2(p)
  dTemp = (Temp2 - Temp1) / (NStep - 1)
  Temp = Temp1

  #find the phi-angle and psi-angle residues
  PhiRes = [i for i in range(len(p)) if i in LoopRes or i-1 in LoopRes]
  PsiRes = [i for i in range(len(p)) if i in LoopRes or i+1 in LoopRes]
  
  #run the steps
  NAcc, NAtt = zeros(1,float), ones(1,float) * 0.1
  for i in range(0, NStep):
    #random move
    AccFact = 0.
    Move = int(2. * random.random())
    Ang = (2*random.random() - 1) * MaxdDih
    OldPos = p.Pos.copy()
    if Move == 0:
      #modify phi
      ResNum = PhiRes[int(random.random() * len(PhiRes))]
      p.RotatePhi(ResNum, Ang) 
    else:
      #modify psi
      ResNum = PsiRes[int(random.random() * len(PsiRes))]
      p.RotatePsi(ResNum, Ang)
    NAtt[0] += 1

    #accept/reject    
    NewE = Score2(p)
    Acc = exp(min(0., -(NewE-OldE) / Temp) + AccFact) > random.random()
    if Acc:
      OldE = NewE
      NAcc[0] += 1
    else:
      p.Pos = OldPos
    if i % 20 == 0 and Verbose:
      print "round %d  %.2f %.2f %.2f  %.2f" % \
       (i, 100*NAcc[0]/NAtt[0], Temp, MaxdDih, OldE)
    if i % 20 == 0:
      Pdb = "../mc%05d.pdb" % i
      #p1.Copy().Concat(p2).WritePdb(Pdb)
      #file("../ene.txt","a").write("%d %.2f %.4f\n" % (i, Temp, OldE))
      
    #update the maximum moves
    if not StepsUpdateMax is None and (i+1) % StepsUpdateMax == 0:
      MaxdPhi = min(MaxdDih * (NAcc[0]/NAtt[0] + 0.1)/0.6, 180)
      NAcc, NAtt = zeros(3,float), ones(3,float) * 0.1
      
    #change temp
    Temp += dTemp


def AlignRMSD2(p1, p2, LoopRes):
  NonLoopRes = [x for x in range(len(p1)) if not x in LoopRes]
  Pos1 = p1.Pos.take(p1.AtomInd(AtomName = ["N", "CA", "C"],
                                ResNum = NonLoopRes), axis=0)
  Pos2 = p2.Pos.take(p2.AtomInd(AtomName = ["N", "CA", "C"],
                                ResNum = NonLoopRes), axis=0)
  return rmsd.RMSD(Pos1, Pos2)

    
def AssembleMC2(PdbFiles, LoopSeqs, Prefix, MaxConf,
              NGen = None, MinRMSD = 2.0, Verbose = True, CleanUp = True):
  if NGen is None: NGen = MaxConf

  #change to a temp directory
  n = 0
  while os.path.isdir("pdbassem%d" % n):
    n += 1
  if n > 49:
    raise IOError, "Error: too many pdbassem temporary folders."
  Dir = "pdbassem%d" % n
  os.mkdir(Dir)
  os.chdir(Dir)

  #make the filename templates
  Tag1 = "pep_%0" + "%d" % int(math.log(NGen, 10) + 1) + "d.pdb"
  Tag2 = "%s_%0" + "%d" % int(math.log(MaxConf, 10) + 1) + "d.pdb"

  #prep the pdb files and make the protein class
  p = protein.ProteinClass()
  LoopRes = []
  for i in range(len(PdbFiles)):
    Pdb = "pep%d.pdb" % (i+1)
    PrepPdb(os.path.join("../", PdbFiles[i]), Pdb)
    p = p - protein.ProteinClass(Pdb = Pdb)
    if i < len(LoopSeqs):
      Seq = LoopSeqs[i]
      LoopRes.extend(range(len(p), len(p) + len(Seq)))
      p = p - protein.ProteinClass(Seq = Seq)
  
  #set protein class options and prep for score
  p.ResAtom = "CB"
  PrepScore2(p, LoopRes)
  
  #make a data holder
  class TempClass:
    def __init__(self, Score, Pdb):
      self.Score = Score
      self.Pdb = Pdb
  Generated = []
  
  #run mc
  NDone = 0
  NLoop = len(LoopSeqs)
  NResLoop = len(LoopRes)
  Scale = max(NResLoop * NLoop**0.5, 10)
  for i in range(NGen):
    if Verbose: print "Performing alignment %d/%d..." % (i, NGen)
    #run a high temperature randomization, starting with last high T config
    for (Temp1, Temp2, NStep) in AnnealSched1:
      RunMCSteps2(p, LoopRes, NStep * Scale, Temp1, Temp2)
    #copy to a new class and anneal
    pnew = p.Copy()
    for (Temp1, Temp2, NStep) in AnnealSched2:
      RunMCSteps2(pnew, LoopRes, NStep * Scale, Temp1, Temp2)
    #compare to the scores
    s = Score2(pnew)
    if Verbose: print "...final score is %.3f" % s
    Pdb = Tag1 % (len(Generated) + 1)
    Generated.append(TempClass(s, Pdb))
    pnew.WritePdb(Pdb)
          
  #sort and reorder files
  Generated.sort(cmp=lambda x,y: cmp(x.Score, y.Score))
  #check rmsd values
  i = 0
  while i < min(len(Generated), MaxConf):
    g = Generated[i]
    p = protein.ProteinClass(Pdb = g.Pdb)
    #filter out other too close RMSDs
    Generated = Generated[:i+1] + [g2 for g2 in Generated[i+1:]
                if AlignRMSD2(p, protein.ProteinClass(Pdb = g2.Pdb), LoopRes) > MinRMSD]
    i += 1
      
  #rename and move the files as sorted
  Generated = Generated[:MaxConf]
  for i in range(0, len(Generated)):
    NewPdb = Tag2 % (Prefix, i+1)
    shutil.move(Generated[i].Pdb, os.path.join("../", NewPdb))
    Generated[i].Pdb = NewPdb

  #change back
  os.chdir("../")
  #write a summary file
  if Verbose: print "Writing summary file"
  s = "ALIGNMENT RESULTS\n"
  s += "Number, Pdb, Score\n"
  for i in range(0, len(Generated)):
    g = Generated[i]
    s += "%-3d %-20s %.2f\n" % (i+1, g.Pdb, g.Score)
  file("%s-alignments.txt" % Prefix, "w").write(s)
  
  #remove temp dir
  if CleanUp: shutil.rmtree(Dir)


  

#-------- ROUTINES --------

def PrepPdb(PdbFile, NewFile):
  """Prepares files for align and splat.
  Returns either the original file name or the temp file name."""
  Pdb = pdbtools.ReturnPdbData(PdbFile)
  Pdb = pdbtools.Standardize(Pdb)
  Pdb = pdbtools.Decap(Pdb)
  Pdb = pdbtools.Renumber(Pdb)
  file(NewFile, "w").write(Pdb)


def MakeLoop(PdbFile1, PdbFile2, LoopSeq, Prefix,
             NPdb = 1, MaxConf = 1, Verbose = True,
             CleanUp = True, DoPdbPrep = True):
  """Uses Splat to make a loop.  PdbFile1 and PdbFile2
  are the names of the pdb files.  LoopSeq is the sequence
  of amino acids (any format).  Prefix is the output pdb
  prefix.  NPdb is the number of Pdb2 files (with names
  like pdb2_001.pdb, etc).  MaxConf is the maximum number
  of loops to generate."""
    
  #prepare the configurations
  Pdb2Prefix = os.path.basename(PdbFile2).replace(".pdb", "")
  if DoPdbPrep:
    PrepPdb(PdbFile1, "temppep1.pdb")
    PrepPdb(PdbFile2, "temppep2.pdb")
    PdbFile1, PdbFile2 = "temppep1.pdb", "temppep2.pdb"
    if NPdb > 1:
      for i in range(0, NPdb):
        PdbFile = PrepPdb("%s_%03d.pdb" % (Pdb2Prefix, i+1),
                          "temppep2_%03d.pdb" % (i+1))

  #prepare loop vars
  LoopSeq = sequence.SeqToList(LoopSeq)
  LoopLen = len(LoopSeq)

  #setup for splat
  Pivots = [i+1 for i in range(LoopLen) if not LoopSeq[i] == "PRO"]
  NPivots = len(Pivots)
  if NPivots < 3:
    print "Error: found less than three non-proline pivots in loop."
    sys.exit(1)
    
  #par down the pivots to three and make the rest perturb
  Pivots = [Pivots[0], Pivots[(NPivots-1)/2], Pivots[-1]]
  Perturb = [i for i in range(1, LoopLen + 1) if not i in Pivots]
  #limit perturb residues
  if len(Perturb) > MaxPerturb:
    dx = float(len(Perturb)) / float(MaxPerturb)
    Ind = [int(dx*i) for i in range(0, MaxPerturb)]
    Ind = [i for i in Ind if i < len(Perturb)]
    Perturb = [Perturb[i] for i in Ind]    
  NPerturb = len(Perturb)
    
  #make the assembly script
  s = '\n'
  s += '%s\n' % SplatPath
  s += '%s\n' % PdbFile1
  s += '%s\n' % PdbFile2
  s += '%s\n' % Prefix
  s += '%d\n' % MaxConf  # closure solutions
  s += '%d, %d, %d\n' % tuple(Pivots)  #pivot res numbers
  s += '%d\n' % NPerturb  #number of perturb res
  s += ", ".join([str(x) for x in Perturb]) + "\n"   #perturb res numbers
  if DEBUG:
    s += '1\n'  # debug info
  else:
    s += '-1\n' # no output
  s += 'C, U\n'
  s += '%d, 0, %d, %d\n' % (MaxSample, NPdb, NRotSample) # sampling size, Perturb hinge (0,none;1,first;2,second;3,both)
  s += '0, 10, 10\n'  # seed ir ir
  s += '%.2f, .95\n' % OverlapDist  # steric cutoff: loop-other; loop-self
  s += "\n".join(LoopSeq) + "\n"  #sequence
  file("assemble.in", "w").write(s)

  #run the code
  if Verbose: print "Running assembly" 
  os.system(AssembleExec + ' assemble.in')
  if Verbose: print "Running SPLAT" 
  os.system(SplatExec)

  #remove excess files
  if CleanUp and not DEBUG:
    TempFiles = ["assemble.in", "assemble.tmp",
                "temppep1.pdb", "temppep2.pdb"]
    for f in TempFiles + glob.glob("temppep2*.pdb"):
      if os.path.isfile(f): os.remove(f)



def Assemble(PdbFile1, PdbFile2, LoopSeq, Prefix, MaxConf = None,
             MinConf = 10, SortByRgPhobic = True, Verbose = True):
  """Assembles PdbFile1 and PdbFile2 by wrapping Albert's and Vageli's code.
  PdbFile1 and PdbFile2 are the names of the pdb files.  LoopSeq is the
  sequence of amino acids (any format).  Prefix is the output pdb prefix.
  MaxConf is the maximum number of configurations to generate.
  Returns a list of the filenames generated."""

  #check prefix
  if Prefix.lower() == "pep":
    print "Error: cannot use 'pep' as prefix."
    sys.exit(1)
  
  #make a temporary path
  n = 0
  while os.path.isdir("pdbassem%d" % n):
    n += 1
  if n > 49:
    raise IOError, "Error: too many pdbassem temporary folders."
  Dir = "pdbassem%d" % n
  os.mkdir(Dir)

  #prepare and copy the configurations over
  PrepPdb(PdbFile1, Dir + "/pep1.pdb")
  PrepPdb(PdbFile2, Dir + "/pep2.pdb")

  #change dirs  
  cwd = os.getcwd()
  os.chdir(Dir)
  
  #make the loop file
  LoopSeq = sequence.SeqToList(LoopSeq)
  file("loop.txt", "w").write(sequence.SeqToAA3(LoopSeq))
  LoopLen = len(LoopSeq)

  #generate alignments
  if Verbose: print "Running alignment"
  os.system(AlignExec + " pep1.pdb pep2.pdb %d %d" % (LoopLen, MinConf))

  #sort through configurations
  PdbFiles = glob.glob("pep2_*.pdb")
  if SortByRgPhobic:
    if Verbose: print "Sorting slignments by hydrophobic radius of gyration"
    l = []
    p1 = protein.ProteinClass(Pdb = "pep1.pdb")
    for f in PdbFiles:
      p = protein.ProteinClass(Pdb = f).Concat(p1)
      RgPhobic = p.RadiusOfGyration(ResInd = p.ResInd(Hydrophobic=True))
      l.append((RgPhobic, f))
    l.sort()
    PdbFiles = [f for (RgPhobic, f) in l]
    #rename old files
    for f in PdbFiles:
      os.rename(f, f.replace("pep2_", "oldpep2_"))
    #trim the list
    if not MaxConf is None: PdbFiles = PdbFiles[:MaxConf]
    #make new files
    for i in range(0, len(PdbFiles)):
      f = PdbFiles[i]
      shutil.copy(f.replace("pep2_", "oldpep2_"), "pep2_%03d.pdb" % (i+1))
  else:    
    #just trim the list
    if not MaxConf is None: PdbFiles = PdbFiles[:MaxConf]
  
  #get number of pdbs
  NPdb = len(PdbFiles)

  #make the loops
  MakeLoop("pep1.pdb", "pep2.pdb", LoopSeq, Prefix,
            NPdb = NPdb, MaxConf = 1, Verbose = True,
            CleanUp = False, DoPdbPrep = False)

  #get the output
  PdbFiles = glob.glob(Prefix + "*.pdb")
  PdbFiles.sort()

  #move the files
  if Verbose: print "Moving final files"
  for i in range(0, len(PdbFiles)):
    shutil.move(PdbFiles[i], "../%s%03d.pdb" % (Prefix, i+1))    

  #cleanup
  os.chdir(cwd)
  if not DEBUG: shutil.rmtree(Dir)

  if Verbose: print "Produced %d assembled structures" % len(PdbFiles)

  return PdbFiles  



#command line running
if __name__ == "__main__":
  PdbFile1, PdbFile2 = sys.argv[1:3]
  LoopSeq = sequence.SeqToList(sys.argv[3])
  Prefix = sys.argv[4]
  Args = scripttools.ParseArgs(sys.argv,
         {"maxconf":50, "debug":False, "overlap":OverlapDist})
  DEBUG = Args["debug"]
  OverlapDist = Args["overlap"]
  MaxConf = Args["maxconf"]
  if "mc" in Args["FLAGS"]:
    AssembleMC(PdbFile1, PdbFile2, LoopSeq, Prefix, MaxConf)
  else:
    Assemble(PdbFile1, PdbFile2, LoopSeq, Prefix, MaxConf)

