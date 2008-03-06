#!/usr/bin/env python

#last modified: 04-05-06

#SHORTCUT SETTINGS:
#Modify these variables without having to change code below:

Seq = "ACE GLY GLY GLY NHE"  #capped sequence
MinTemp = 270.0    #minimum replica temperature 
MaxTemp = 700.0    #maximum replica temperature
NCycle = 1000      #total number of md/exchange cycles
NSteps = 5000      #number of md steps per cycle
NRep = 18          #total number of replicas
CapSeq = False     #whether or not to cap the sequence with ACE and NME
Restart = True     #False means ignore any restart data if it is there;
                   #True means load if found; otherwise, start anew


#MAIN PROGRAM

import os, math
import socketring, sequence
import rexsock

#program globals (don't change)
kB= (1.381 * 6.02214) / 4184.0
Temp = []

#convert seq to a list
Seq = sequence.SeqToList(Seq)

#This defines the acceptance function for swaps, used by the rexsock module.
#a and b are the simulations of the individual "states" (e.g., temperatures).
#The example's variables TEMPSET and EPOT are the setpoint temperature and
#potential energy.  The suffixes "1" and "2" indicate those variables at the
#beginning and end of the trajectory segment.  Other variables available are:
#EKIN (kinetic energy), EEL (electrostatic energy), EREST (restraint energy),
#EVDW (van der Waals energy), ESURF (surface energy), EBOND (bond energy),
#EANG (angle energy), EDIH (dihedral energy), TIME (simulation time), TEMP
#(kinetic temperature), PRES (virial pressure), all of which have suffixes
#of either "1" or "2".
def PaccFunc(a, b):
  dBeta = 1./(kB * a["TEMPSET"]) - 1./(kB * b["TEMPSET"])
  dE = a["EPOT2"] - b["EPOT2"]
  return math.exp( min([0., dE*dBeta]) )


#This is the main loop of the program.
if __name__ == "__main__":

  #This creates a socket ring with a master server and clients; the
  #file containing ring descriptors is saved in path "."
  sr = socketring.RingClass(".")

  #See if we are the main server

  if sr.IsServer:  
  
    #This initializes the replica exchange program in the path "." and
    #designates the first computer to run the program as the main server
    #which will dole out job tasks.  If you want to monitor the values
    #of specific variables in each state (possibilities listed above), you
    #can add a third argument to the RexClass constructor which is a dict
    #of variable name (key) + file to save in (value) pairs, for example:
    #rx = rexsock.RexClass(NRep, ".", {"TEMP":"temp.txt", "PRES":"pres.txt"})

    rx = rexsock.RexClass(SockRing = sr, BasePath = ".", N = NRep, Restart = Restart)

    #tell the socket ring how many client processors we need
    sr.MaxNeeded = NRep

    #Setup the progress.txt file
    ProgFile = "progress.txt"    
    file(ProgFile,"w").write("starting up...\n")

    #Setup the temperatures.
    if rx.NReplica == 1:
      Temp = [MinTemp]
    else:
      dlnT = (math.log(MaxTemp)-math.log(MinTemp)) / float(rx.NReplica-1)
      Temp = [MinTemp*math.exp(float(i)*dlnT) for i in range(0,rx.NReplica)]
    
    #If a restart was not detected, build the systems up from scratch.
    #The function SimIter will iterate a function over each state
    #replica.  SimSetVar will set a simulation variable in each replica
    #using an array whose size is the number of replicas.  SimFarm will
    #iterate a function string over each replica by allocating the task
    #to the available processors.
    if not rx.Restart:
      rx.SimIter(lambda s: s.SysInitSeq(Seq, Cap = CapSeq))
      rx.SimSetVar("STEPSMD", [NSteps for i in range(0,rx.NReplica)])
      rx.SimSetVar("TEMPSET", Temp)
      rx.SimIter(lambda s: s.SysBuild())
      rx.SimFarm("SIM.RunMin()")
    
    #Run the cycles.  Labels is used to identify the states in
    #the progress file.  rx.StartRex sets up everything necessary
    #for the rex simulation.
    Labels = ["%.1f K" % (t,) for t in Temp]
    rx.StartRex(PaccFunc, NCycle, Labels, ProgFile)

    #Continue running until the REX simulation is finished
    while not rx.RunRex():
      #let the socket ring run tasks for 1 second
      sr.RunTasks(1.)
      
    #Close down the socket ring and tell clients to go home.
    sr.Finish()

  #If the current processor is not the main server, run as a client and
  #connect to the server.
  else:

    sr.RunClient()
    



