#!/usr/bin/env python

#LAST MODIFIED: 05-26-06

Usage = """Checks a rex simulation for bad files and tries to fix.

Usage     : checkrex.py MODE REXPATH

MODE      : either 'q' for quick or 'f' for full
REXPATH   : directory of the replica exchange path (use "all" for r*/)
"""

ExceptErrors = (IOError, ValueError, LookupError, EOFError)


import sys
#check for instructions
if len(sys.argv) == 1:
  print Usage
  sys.exit()

import rexsock, mdsim, coords
import os, sys, cPickle, shutil, gzip

if sys.argv[1].lower() == "f":
  Full = True
else:
  Full = False

if sys.argv[2].lower() == "all":
  Paths = [p for p in os.listdir(".") if os.path.isdir(p) and p.startswith("r")]
else:
  Paths = sys.argv[2:]

for p in Paths:
  print "Checking path %s" % p
  WorkPath = os.path.join(p,"work")
  DataPath = os.path.join(p,"data")
  #first get stats
  try:
    f  = file(os.path.join(DataPath, "rexstats.txt"))
    while True:
      s = f.readline()
      if len(s) == 0: break
      if "Number of replicas" in s:
        NRep = int(s.split()[-1])
      elif "Number of cycles" in s:
        NCycle = int(s.split()[-1])
      elif "Sequence" in s:
        Seq = []
        while True:
          s = f.readline().strip()
          if len(s) == 0: break
          Seq.append(s)
      elif "Temperatures" in s:
        Temps = []
        while True:
          s = f.readline().strip()
          if len(s) == 0: break
          Seq.append(float(s))
  except:
    print "ERROR: Cannot find rexstats.txt.  Please enter values."
    NRep = int(raw_input("NRep: "))
    NCycle = int(raw_input("NCycle: "))
    Seq = raw_input("Seq (3 letter codes): ").split()
    Temps = [float(x) for x in raw_input("Temps: ").split()]
  try:
    f.close()
  except:
    pass

  #first check worspace
  print "Checking workspace..."
  #check restart file
  RstFile = os.path.join(WorkPath, "restart.dat")
  try:
    dat = cPickle.load(file(RstFile,"r"))
  except:
    print "ERROR: Problem with restart.dat file."
    cPickle.dump([('NCycle', NCycle)], file(RstFile, "w"))
    print "FIXED: restart.dat"
  #check each work path
  for i in range(0,NRep):
    ThisPath = os.path.join(WorkPath, str(i))
    if not os.path.isdir(ThisPath):
      print "ERROR: work path %s does not exist" % ThisPath
      continue
    try:
      dat = cPickle.load(file(os.path.join(ThisPath, "sim.dat")))
    except:
      print "ERROR: problem with sim.dat in %s" % ThisPath
    fn = os.path.join(ThisPath, "current.crd")
    dat = None
    if os.path.isfile(fn):
      if os.path.getsize(fn) > 0:
        dat = coords.GetRstCoords(os.path.join(ThisPath, "current.crd"))
    if dat is None:
      print "ERROR: problem with current.crd in %s" % ThisPath
      #try to use another crd
      for OtherFile in ["mdout.crd","minout.crd"]:
        f = os.path.join(ThisPath, OtherFile)
        if os.path.isfile(f):
          if os.path.getsize(f) > 0:
            shutil.copy(f, os.path.join(ThisPath,"current.crd"))
            print "FIXED: current.crd"
            break

  #check data path
  for i in range(0,NRep):
    sys.stdout.flush()
    print "Checking data for state %d..." % i
    fn = os.path.join(DataPath, "%d.prmtop.parm7" % i)
    if not os.path.isfile(fn):
      print "ERROR: could not find %s" % fn
      fn2 = os.path.join(WorkPath, str(i), "prmtop.parm7")
      if os.path.exists(fn2):
        shutil.copy(fn2, fn)
        print "FIXED: %s" % fn
    for suff in [".mdtrj.crd.gz",".mdene.txt.gz"]:
      fn = os.path.join(DataPath, str(i)+suff)
      if os.path.isfile(fn):
        try:
          f = gzip.GzipFile(fn)
          if Full:
            while len(f.readline()) > 0:
              pass
        except:
          print "ERROR: problem with %s" % fn
        try:
          f.close()
        except:
          pass
      else:
        print "ERROR: could not find %s" % fn

  print "Done checking.\n"
    
      
    
  

