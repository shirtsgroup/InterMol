#!/usr/bin/env python

#LAST MODIFIED: 08-10-06

import os, sys, getpass, socket

Usage = """Submits or deletes groups of jobs at once.

Usage    : multijob.py sub N JOBFILE
      OR : multijob.py delnums STARTNUM STOPNUM
      OR : multijob.py deljob JOBNAME
      
JOBFILE  : name of job script to submit
N        : number of times to submit job
STARTNUM : starting number of jobs to delete
STOPNUM  : stopping number of jobs to delete
JOBNAME  : name of job to delete 
"""

#strings and token numbers for different computers
#numbers are (substr, killstr, statstr, jobnum, jobname, user)
HostData = {"baton2.compbio.ucsf.edu":("qsub ", "qdel ", "qstat", 0, 2, 3),
            "chef.compbio.ucsf.edu":("qsub ", "qdel ", "qstat", 0, 2, 3),
            "zwicky.compbio.ucsf.edu":("qsub ", "qdel ", "qstat", 0, 1, 2),
            "tuna":("bsub < ", "bkill ", "bjobs", 0, 6, 1),
            "tunb":("bsub < ", "bkill ", "bjobs", 0, 6, 1),
            "tunc":("bsub < ", "bkill ", "bjobs", 0, 6, 1),
            "tund":("bsub < ", "bkill ", "bjobs", 0, 6, 1),
            "tune":("bsub < ", "bkill ", "bjobs", 0, 6, 1) }

Host = socket.gethostname().strip()
if Host in HostData:
  SubCmd, KillCmd, StatCmd, tnu, tnm, tu = HostData[Host]
else:
  print "I do not know how to parse job information for host %s." % Host
  sys.exit()
  

if len(sys.argv) > 1:
  
  if sys.argv[1] == "sub":
    N = int(sys.argv[2])
    Job = sys.argv[3].strip()
    if os.path.isfile(Job):
      for i in range(0,N):
        os.system(SubCmd + Job)
    else:
      print "Job file not found."

  elif sys.argv[1] == "deljob":
    Job = sys.argv[2].strip()
    Usr = getpass.getuser()
    f = os.popen(StatCmd)
    s = f.readlines()
    f.close()
    for l in s:
      Tokens = l.split()
      try:
        nm = Tokens[tnm]
        nu = Tokens[tnu]
        u = Tokens[tu]
      except IndexError:
        continue
      if u==Usr and nm==Job:
        os.system(KillCmd + nu)

  elif sys.argv[1] == "delnums":
    for i in range(int(sys.argv[2]), int(sys.argv[3])+1):
      os.system(KillCmd + str(i))

  else:
    print "Command not recognized."
    
else:
  print Usage



