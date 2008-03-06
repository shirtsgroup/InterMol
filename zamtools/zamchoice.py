#!/usr/bin/env python

#last modified: 04-29-2007

from numpy import *
import os, math, cPickle, time, sys
import zamdata, sequence

#GLOBALS
NResGap = 1             #stride in num of residues between initial 8mer frags
MaxFragsDflt = 10       #maximum number of new fragments to create in one round after 16mers
MaxChoicesDflt = 100    #maximum number of assembly choices to report
NPredictions = 10       #maximum number of full chain predictions to make
NResLoop = 4            #won't penalize loops of this length or less
AssemResPenalty = 5.0   #penalty per res per num of previous assemblies in score 
COWeight = 10.75        #weight multiplied by contact order and added to score

#Status strings:
# "stuck" - simulation cannot proceed any further
# "done" - simulation is done
# "run" - run new frags
# "stop" - user has indicated to stop running new fragments


def BestCmp(x,y):
  """Defines a comparison function for sorting assembly choices."""
  LoopLen1 = max(x[0] - NResLoop, 0)
  LoopLen2 = max(y[0] - NResLoop, 0)
  c1 = cmp(LoopLen1, LoopLen2)
  if c1 == 0:
    c2 = -cmp(x[1], y[1])
    return c2
  else:
    return c1

def BestSort(zd, Best):
  """Sorts assemblies and prunes them for redundant attempts covering the same residues."""
  Best.sort(cmp = BestCmp)
  NAssem = zeros(len(zd.Seq), int)
  NewBest = []
  while len(Best) > 0:
    #run through and find the minimum number of overlapping res
    (LoopLen, Score, x1, x2, a, b, Desc) = Best.pop(0)
    if average(NAssem[a:b+1]) < MaxAssemPerRes:
      NewBest.append((LoopLen, Score, x1, x2, a, b, Desc))
      NAssem[a:b+1] += 1
  return NewBest

#this is the old version, with maximum number of assemblies per res
MaxAssemPerRes = 4.0
def BestSortOld(zd, Best):
  """Sorts assemblies and prunes them for redundant attempts covering the same residues."""
  Best.sort(cmp = BestCmp)
  NewBest = []
  while len(Best) > 0:
    #add the top scoring to the new list
    NewBest.append(Best.pop(0))
    #update the scores for the other fragments
    a, b = NewBest[-1][4], NewBest[-1][5]
    ResPenalty = zeros(len(zd.Seq), int)
    ResPenalty[a:b+1] = AssemResPenalty
    for x in Best: x[1] += sum(ResPenalty[x[4]:x[5]+1])
    #resort
    Best.sort(cmp = BestCmp)
  return NewBest
 

def GetNextActions(zd, cd, sd, MaxFrags = None, MaxChoices = None,
                   AssembleByFrag = None):
  "Creates new fragments for growth and/or assembly."

  if MaxFrags is None: MaxFrags = MaxFragsDflt
  if MaxChoices is None: MaxChoices = MaxChoicesDflt
  if AssembleByFrag is None: AssembleByFrag = zd.AssembleByFrag

  #CHECK TO SEE IF FRAGS ARE ALREADY RUNNING
  if len(zd.GetRunFrags()) > 0: return "run"
 
  #CHECK TO SEE IF THE FULL PROTEIN HAS BEEN SIMULATED
  #********NEED TO FILL THIS IN ********
  #if SOMETHING: return "done"  

  #MAKE 8MERS
  if len(zd.GetFrags(Length = 8)) == 0:
    #make 8mers; ensure both end frags are there
    zd.Stage = "8mers"
    SegList = [(i, i+7) for i in range(0, len(zd.Seq) - 7, NResGap)]
    LastSeg = (len(zd.Seq) - 8, len(zd.Seq) - 1)
    if not LastSeg in SegList: SegList.append(LastSeg)
    NewFrags = [zd.AddFrag(a, b, [], "grow") for (a,b) in SegList
                if not zd.SegInFrags(a, b)]
    #return to run
    if len(NewFrags) > 0: return "run"

  #MAKE 12MERs
  if len(zd.GetFrags(Length = 12)) == 0:
    #make 12mers from top 2/3 of 8mers
    zd.Stage = "12mers"
    #get 8mer fragments
    Frags = zd.GetFrags(Length = 8)
    Frags.sort(reverse=True, cmp=lambda x,y: cmp(x.FragScore, y.FragScore))
    N = int(float(len(Frags)) * 2./3.)
    Frags = Frags[:N]
    NewFrags = []
    for Frag in Frags:
      a, b = Frag.StartRes - 2, Frag.StopRes + 2
      if a < 0:
        a, b = 0, b - a
      elif b >= len(zd.Seq):
        a, b = a + len(zd.Seq) - b - 1, len(zd.Seq) - 1
      a, b = max(a, 0), min(b, len(zd.Seq)-1)
      if not zd.SegInFrags(a, b):
        NewFrag = zd.AddFrag(a, b, [(Frag, -1)], "grow")
        NewFrags.append(NewFrag)
    #return to run
    if len(NewFrags) > 0: return "run"

  #MAKE 16MERS
  if len(zd.GetFrags(Length = 16)) == 0:
    #make 16mers from top 2/3 of 12mers
    zd.Stage = "16mers"
    #get 12mer fragments
    Frags = zd.GetFrags(Length = 12)
    Frags.sort(reverse=True, cmp=lambda x,y: cmp(x.FragScore, y.FragScore))
    N = int(float(len(Frags)) * 2./3.)
    Frags = Frags[:N]
    NewFrags = []
    for Frag in Frags:
      a, b = Frag.StartRes - 2, Frag.StopRes + 2
      if a < 0:
        a, b = 0, b - a
      elif b >= len(zd.Seq):
        a, b = a + len(zd.Seq) - b - 1, len(zd.Seq) - 1
      a, b = max(a, 0), min(b, len(zd.Seq)-1)
      #make the new fragments
      if not zd.SegInFrags(a, b):
        NewFrag = zd.AddFrag(a, b, [(Frag, -1)], "grow")
        #enable swapping restraints for 16mers to generate different topologies
        NewFrag.CalcSwapRest = True
        NewFrag.CalcSwapRestTime = True
        NewFrags.append(NewFrag)
    #return to run
    if len(NewFrags) > 0: return "run"      

  zd.Stage = "assembly"

  if AssembleByFrag:
    
    #ASSEMBLE TWO FRAGMENTS (ALL SNAPSHOTS OF EACH)
    Best = []
    NewFrags = []
    print "Analyzing pairs of %d fragments." % len(zd.Frags)
    #sort through fragment pairs

    for (m, Frag1) in enumerate(zd.Frags):
      Snaps1 = sd.GetSnaps(Frag1)
      
      for (n, Frag2) in enumerate(zd.Frags[m+1:]):
        #check to see if this is useful
        if Frag1.ContainsFrag(Frag2) or Frag2.ContainsFrag(Frag1): continue
        #check to see if we've tried this before
        if zd.AnyFragsWithParents([(Frag1, -1), (Frag2, -1)]): continue
        Snaps2 = sd.GetSnaps(Frag2)

        #sort through pairs of snapshots
        Score = 0.
        for (i, Snap1) in enumerate(Snaps1):
          for (j, Snap2) in enumerate(Snaps2):
            #get the score
            Weight = Frag1.SnapWeight[i] * Frag2.SnapWeight[j]
            Score += zamdata.SnapPairScore(Snap1, Snap2, cd.Scores, COWeight) * Weight

        #get the residue lengths
        LoopLen = max(Frag1.StartRes - Frag2.StopRes - 1, Frag2.StartRes - Frag1.StopRes - 1)
        OldFragLen = Frag1.StopRes - Frag1.StartRes + Frag2.StopRes - Frag2.StartRes
        NewFragLen = OldFragLen + LoopLen

        #get the text descriptor
        if Frag1.StartRes < Frag2.StartRes:
          f1, f2 = Frag1, Frag2
        else:
          f1, f2 = Frag2, Frag1
        Desc = "%-10s %-5d %-5d %-5d %-10s %-5d %-5d %-5d %-5d %.4e" \
          % (f1.Name(), f1.StartRes+1, f1.StopRes+1, 1 + f1.StopRes - f1.StartRes,
             f2.Name(), f2.StartRes+1, f2.StopRes+1, 1 + f2.StopRes - f2.StartRes,
             LoopLen, Score)
        #add to the best list:
        Best.append((LoopLen, Score, f1, f2, f1.StartRes, f2.StopRes, Desc))
        #sort and prune list
        Best = BestSort(zd, Best)
        Best = Best[:MaxChoices]
        
    #write the choices out
    s = "%-10s %-5s %-5s %-5s %-10s %-5s %-5s %-5s %-5s %s" % \
        ("Frag1", "Start", "Stop", "Len", "Frag2", "Start", "Stop", "Len",
         "NLoop", "Score")
    for (LoopLen, Score, Frag1, Frag2, StartRes, StopRes, Desc) in Best:
      s += "\n" + Desc
    s = s + "\n"
    file("fragchoices.txt", "w").write(s)
    #trim the list
    Best = Best[:MaxFrags]
    #make the fragments
    for (LoopLen, Score, Frag1, Frag2, StartRes, StopRes, s) in Best:
      NewFrag = zd.AddFrag(StartRes, StopRes, [(Frag1, -1), (Frag2, -1)], "assem")
      NewRest = [(a,b) for (s,a,b) in Frag1.PropRest + Frag2.PropRest]
      NewRest = [x for (i, x) in enumerate(NewRest) if not x in NewRest[i+1:]]
      NewFrag.PriorRest = sorted(NewRest)
      NewFrags.append(NewFrag)
      #check for full sequence
      if StartRes == 0 and StopRes == len(zd.Seq) - 1:
        NewFrag.CapN = False
        NewFrag.CapC = False
    #return to run
    if len(NewFrags) > 0: return "run"
        

  else:
    
    #ASSEMBLE TWO SPECIFIC SNAPSHOTS FROM THE STRUCTURE DATABASE
    Best = []
    NSnap = len(sd.Snaps)
    NewFrags = []
    print "Analyzing pairs of %d snapshots." % NSnap
    #sort through snapshot pairs
    t = "%-10s %-5s %-5s %-5s %s" % ("Snap", "Start", "Stop", "Len", "Score")
    for i in range(0, NSnap):
      Snap1 = sd.Snaps[i]
      t += "\n%-10s %-5s %-5s %-5s %s" % (Snap1.Name, Snap1.StartRes+1, Snap1.StopRes+1,
                                          Snap1.StopRes - Snap1.StartRes + 1, Snap1.Score)
      Frag1 = Snap1.Parent[0]
      for j in range(i+1, NSnap):
        Snap2 = sd.Snaps[j]
        Frag2 = Snap2.Parent[0]
        #skip if same parent
        if Frag1 == Frag2: continue
        #make sure we haven't tried this before
        if zd.AnyFragsWithParents([Snap1.Parent, Snap2.Parent]): continue
        #make sure one fragment isn't inside the other
        if Snap1.ContainsSnapRes(Snap2) or Snap2.ContainsSnapRes(Snap1): continue
        #get the residue lengths
        LoopLen = max(Snap1.StartRes - Snap2.StopRes - 1, Snap2.StartRes - Snap1.StopRes - 1)
        OldSnapLen = Snap1.StopRes - Snap1.StartRes + Snap2.StopRes - Snap2.StartRes
        NewSnapLen = OldSnapLen + LoopLen
        #get the score
        Score = zamdata.SnapPairScore(Snap1, Snap2, cd.Scores, COWeight)
        #get the text descriptor
        if Snap1.StartRes < Snap2.StartRes:
          s1, s2 = Snap1, Snap2
        else:
          s1, s2 = Snap2, Snap1
        Desc = "%-10s %-5d %-5d %-5d %-10s %-5d %-5d %-5d %-5d %.4e" \
          % (s1.Name, s1.StartRes+1, s1.StopRes+1, 1 + s1.StopRes - s1.StartRes,
             s2.Name, s2.StartRes+1, s2.StopRes+1, 1 + s2.StopRes - s2.StartRes,
             LoopLen, Score)
        #add to the best list:
        Best.append((LoopLen, Score, s1, s2, s1.StartRes, s2.StopRes, Desc))
        #sort and prune list
        Best = BestSort(zd, Best)
        Best = Best[:MaxChoices]
    #write the choices out
    s = "%-10s %-5s %-5s %-5s %-10s %-5s %-5s %-5s %-5s %s" % \
        ("Snap1", "Start", "Stop", "Len", "Snap2", "Start", "Stop", "Len",
         "NLoop", "Score")
    for (LoopLen, Score, Snap1, Snap2, StartRes, StopRes, Desc) in Best:
      s += "\n" + Desc
    s = s + "\n\n" + t
    file("fragchoices.txt", "w").write(s)
    #trim the list
    Best = Best[:MaxFrags]
    #make the fragments
    for (LoopLen, Score, Snap1, Snap2, StartRes, StopRes, s) in Best:
      NewFrag = zd.AddFrag(StartRes, StopRes, [Snap1.Parent, Snap2.Parent], "assem")
      NewRest = [(a,b) for (s,a,b) in Snap1.Parent.PropRest if (a,b) in Snap1.Contacts] \
              + [(a,b) for (s,a,b) in Snap2.Parent.PropRest if (a,b) in Snap2.Contacts]
      NewRest = [x for (i, x) in enumerate(NewRest) if not x in NewRest[i+1:]]
      NewFrag.PriorRest = sorted(NewRest)
      NewFrags.append(NewFrag)
      #check for full sequence
      if StartRes == 0 and StopRes == len(zd.Seq) - 1:
        NewFrag.CapN = False
        NewFrag.CapC = False
    #return to run
    if len(NewFrags) > 0: return "run"

  #WE'RE STUCK
  return "stuck"


def FinalizeActions(zd, cd, sd):
  "Finishes up fragments after they have been run."
  for Frag in zd.GetRunFrags():
    Frag.Stage = "done"
    
