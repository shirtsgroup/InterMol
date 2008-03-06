#!/usr/bin/env python

#LAST MODIFIED: 03-07-07

import os, math, cPickle, time, sys, copy
import rexsock, socketring, sequence, mdsim
import zamdata, zamchoice, zamanal, zamrun, zamcontact, zamrex
import protein, sequence

#GLOBALS
ConserveDisk = False   #True will delete large files to save space
ShowRest = False       #whether or not to show proposed restraints
ExcessClients = 1      #excess number of socket clients to keep

#global parameters (don't change)
zd = None
sd = None
cd = None
RestartData = ["zd"]
NewStage = False


#-------- PROGRESS FILE --------
def Prog(s):
  "Updates a progress file."
  file("progress.txt", "a").write(s + "\n")

def ProgInit():
  "Initializes the progress file."
  file("progress.txt", "w").write("Started zam simulation on " + time.strftime("%c") + "\n")
  

#------- RESTART FUNCTIONS --------

def RestartWrite():
  "Writes the restart data."
  global zd, RestartData
  #first write the new file
  Data = [(itm, globals()[itm]) for itm in RestartData]
  cPickle.dump(Data, file("restartnew.dat", "w"))
  #now copy to the old one
  if os.path.isfile("restart.dat"): os.remove("restart.dat")
  os.rename("restartnew.dat", "restart.dat")

def RestartUpdate(Frag, Ind = None):
  "Updates a fragment restart index and writes the restart file."
  global zd, RestartData, NewStage
  if Ind is None:
    Frag.RestartInd += 1
  else:
    Frag.RestartInd = Ind
  RestartWrite()
  NewStage = True
 
def RestartLoad():
  "Reads the restart data."
  global zd, RestartData
  if os.path.isfile("restart.dat"):
    Data = cPickle.load(file("restart.dat","r"))
    for itm, dat in Data:
      globals()[itm] = dat
    #fix for old data
    zd.UpdateVer()
    return True
  else:
    return False


#-------- MAX NUM CLIENTS NEEDED FOR A FRAG --------
def ClientsNeeded(Frag):
  "Returns the number of socket clients needed."
  if Frag.RestartInd > 3:
    return 0
  elif Frag.RestartInd > 1:
    return 1
  elif Frag.RestartInd == 1:
    return Frag.rx.NReplica
  else:
    return 100
  

#-------- MAIN PROGRAM LOOP --------

#initialize zam data variable
zd = zamdata.ZamDataClass()

#check for a restart file
RestartLoad()

#check for a setup file
if os.path.isfile("zamsetup.py"):
  execfile("zamsetup.py")

#check that a sequence has been specified
if len(zd.Seq) == 0:
  print "Sequence has not been specified."
  sys.exit()

#initialize the socketring class, which will decide who is the main server
sr = socketring.RingClass()

if sr.IsServer:

  ProgInit()
  Prog("Sequence is %d residues: %s" % (len(zd.Seq), sequence.SeqToAA1(zd.Seq)))

  #initializd the contact and snapshot databases, in that order
  Prog("Initializing contact database.")
  cd = zamcontact.ContactDatabaseClass(zd)
  Prog("Initializing snapshot database.")
  sd = zamdata.SnapDatabaseClass(zd, cd.Scores)

  #loop over fragments
  Cont = True
  while Cont:
  
    #figure out what to do next
    Stat = ""
    if not zd.AnyRunFrags():
      if zd.UpdateFrags:
        Prog("\nGetting next set of fragments")
        Stat = zamchoice.GetNextActions(zd, cd, sd)
        if Stat == "stuck":
          Cont = False
          Prog("Stuck in ZAM simulation")
        elif Stat == "done":
          Cont = False
          Prog("Done with ZAM simulation")
          #REMOVE ALL REST AND RELAX??    
        elif Stat == "stop":
          Cont = False
          Prog("ZAM simulation has been stopped")
      else:
        Cont = False
        Prog("No more fragments to run")
 
    #get the list of new fragments
    FragList = zd.GetRunFrags()   
    for Frag in FragList:
      #make a string prefix for progress file
      Frag.ProgPrefix = Frag.Name() + ": "
      #make the path
      Frag.MakePath()
      #update the progress file
      s = ", ".join([f.Name() for (f,ind) in Frag.Parents])
      if len(s) > 0: s = "; parents are " + s
      s = "INIT" + s
      Prog(Frag.ProgPrefix + s)


    #FRAGMENT RUN SCHEDULE
    Update = True
    TimeLastUpdate = 0.
    RunFragList = copy.copy(FragList)
    FirstRound = True
    while True:
      #update the progress file
      if Update or time.time() - TimeLastUpdate > 3600:
        l = []
        for Frag in RunFragList:
          if Frag.rx is None or Frag.rx.SecRemainEst == 0:
            l.append(Frag.Name())
          else:
            l.append("%s (%.1f hrs)" % (Frag.Name(), Frag.rx.SecRemainEst/3600))
        s = "Active fragments: " + ", ".join(l)
        Prog(s)
        Update = False
        TimeLastUpdate = time.time()
        
      #loop through fragments
      for Frag in RunFragList:
        sr.RunTasks(0.1)
        NewStage = False
        
        #INITIAL SETUP AND CONFIGURATION MAKING
        if Frag.RestartInd == 0:
          #make the initial config pdb files
          if NewStage or FirstRound: Prog(Frag.ProgPrefix + "Put initial configuration making in task queue")
          InitConfDone, Frag.Task = zamrun.RunInitConfTask(sr, Frag, Frag.Task)
          if InitConfDone: RestartUpdate(Frag)

        #REX
        if Frag.RestartInd == 1:
          if Frag.rx is None:
            #determine the num of proc and set the REX properties for the new fragment
            Prog(Frag.ProgPrefix + "Initializing replica exchange")
            Frag.rx = rexsock.RexClass(SockRing = sr)
            zamrex.PrepFrag(Frag.rx, zd, Frag)
          if NewStage or FirstRound: Prog(Frag.ProgPrefix + "Put replica exchange in task queue")
          RexDone = zamrun.RunFrag(Frag.rx, zd, Frag)
          if RexDone: RestartUpdate(Frag)

        #ANALYSIS
        if Frag.RestartInd == 2:
          if NewStage or FirstRound: Prog(Frag.ProgPrefix + "Put analysis in task queue")
          AnalDone, Frag.Task = zamanal.RunAnalTask(sr, Frag, Frag.Task)
          if AnalDone: RestartUpdate(Frag)

        #SNAPSHOTS
        if Frag.RestartInd == 3:
          if NewStage or FirstRound: Prog(Frag.ProgPrefix + "Put snapshot extraction in task queue")
          SnapDone, Frag.Task = zamanal.RunSnapTask(sr, Frag, Frag.Task)
          if SnapDone: RestartUpdate(Frag)

        #CLEANUP
        if Frag.RestartInd == 4:
          #erase old data
          Prog(Frag.ProgPrefix + "DONE; cleaning up")
          zamanal.CleanUp(Frag, ConserveDisk)
          zamrex.CleanUp(Frag.rx, zd, Frag, ConserveDisk)
          RestartUpdate(Frag)
          Update = True
          
      #now sort through and remove the done fragments
      RunFragList = [Frag for Frag in FragList if not Frag.RestartInd > 3]
      #if there are none left, exit the loop
      if len(RunFragList) == 0: break
      #run tasks if they're there
      sr.RunTasks(0.1)
      FirstRound = False     
      #now update the total number of clients needed
      sr.MaxNeeded = sum([ClientsNeeded(Frag) for Frag in FragList]) + ExcessClients
      

    #PROCESS RESULTS FROM ALL FRAGMENTS
    if Cont:
      Prog("Updating contact scores")
      for Frag in FragList: cd.UpdateFrag(zd, Frag)
      cd.UpdateScores()
      cd.WriteScores("logit.txt")
      cd.WriteHtml("logit.html")
      Prog("Looking for new snaps")
      for Frag in FragList: sd.AddSnaps(Frag)
      Prog("Updating snap scores")
      sd.UpdateScores(cd.Scores)
      Prog("Removing redundant snaps")
      sd.RemoveRedundant()
      Prog("Snapshot database has %d structures" % len(sd))
      #finalize everything
      zamchoice.FinalizeActions(zd, cd, sd)

  print "End of run."
  Prog("End of run.")
  RestartWrite()
  sr.Finish()


else:

  #run client and stop if something unexpected happened
  if not sr.RunClient():
    print "There were errors connecting with the REX server."
    print "Terminating job."
    sys.exit(0)

 