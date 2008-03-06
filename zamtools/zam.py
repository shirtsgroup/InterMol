#!/usr/bin/env python

#last modified: 05-11-07

Usage = """Displays or modifies data from a zam simulation.
Capitalized tokens are parameters.  Brackets indicate optional.

Usage     : zam.py help
   OR     : zam.py initializerun SEQ
   OR     : zam.py restartrun
   OR     : zam.py show [FRAGNUM]
   OR     : zam.py showglobal
   OR     : zam.py showall [FRAGNUM]
   OR     : zam.py showanal [FRAGNUM]
   OR     : zam.py showrest [FRAGNUM]
   OR     : zam.py shownote [FRAGNUM]
   OR     : zam.py decidenext [NFRAG] [byfrag]
   OR     : zam.py clearglobalrest
   OR     : zam.py addglobalrest RESTLIST
   OR     : zam.py delglobalrest RESTLIST
   OR     : zam.py addglobalss RESLIST SS
   OR     : zam.py clearglobalss
   OR     : zam.py edit VARNAME VALUE
   OR     : zam.py newfrag STARTRES STOPRES [PARENTS [REXTIME]]
   OR     : zam.py growfrag FRAGNUM [NUMRES | {NNUMRES CNUMRES}]
   OR     : zam.py duplicatefrag FRAGNUM
   OR     : zam.py makeseedfrags SEEDLEN SEEDGAP [STARTRES STOPRES]
   OR     : zam.py makeconcensusfrags STARTRES STOPRES PARENT [MAXCHILD|RESTLIST]
   OR     : zam.py combineconcensusfrags PARENTS [REXTIME]
   OR     : zam.py fragparents FRAGNUM PARENTS
   OR     : zam.py delfrag FRAGNUM
   OR     : zam.py delchildren FRAGNUM
   OR     : zam.py estimatereplicas FRAGNUM
   OR     : zam.py fragreplicas FRAGNUM NREPLICAS
   OR     : zam.py fragrextime FRAGNUM REXTIME 
   OR     : zam.py fraganaltime FRAGNUM ANALTIME [SKIPTIME]
   OR     : zam.py fraganalopt FRAGNUM {punchon|punchoff}
   OR     : zam.py restartfrag FRAGNUM [REXTIME]
   OR     : zam.py reinitfrag FRAGNUM [REXTIME]
   OR     : zam.py fragstage FRAGNUM STAGE
   OR     : zam.py clearrest FRAGNUM
   OR     : zam.py addrest FRAGNUM RESTLIST
   OR     : zam.py delrest FRAGNUM RESTLIST
   OR     : zam.py copyrest DESTFRAGNUM SOURCEFRAGNUM
   OR     : zam.py addss FRAGNUM RESLIST SS
   OR     : zam.py clearss FRAGNUM
   OR     : zam.py rigidifyfrag FRAGNUM RESLIST
   OR     : zam.py unrigidifyfrag FRAGNUM
   OR     : zam.py anchorfrag FRAGNUM BBRESLIST AARESLIST
   OR     : zam.py unanchorfrag FRAGNUM
   OR     : zam.py addswaprest FRAGNUM SWAPTIME SWAPRESTLIST
   OR     : zam.py repelsaltbridges FRAGNUM {on|off}
   OR     : zam.py scalerest FRAGNUM {on|off}
   OR     : zam.py addumbrellarest FRAGNUM RESTLIST
   OR     : zam.py umbrelladistances FRAGNUM DIST1 DIST2 ...
   OR     : zam.py clearumbrellarest FRAGNUM
   OR     : zam.py fragpairs FRAGNUM {all | hydrophobic}
   OR     : zam.py fragcaps FRAGNUM {none | both | ncap | ccap}
   OR     : zam.py fragnote FRAGNUM [NOTETEXT]
   OR     : zam.py editfrag FRAGNUM VARNAME VALUE
   OR     : zam.py compressdata FRAGNUM FILENAME [anal|prep] [database]
   OR     : zam.py summarize FRAGNUM KEYNAME
   OR     : zam.py deldata FRAGNUM [workonly|dataonly]
   OR     : zam.py savescores FRAGNUM FILEPREFIX

SEQ       : protein sequence
FRAGNUM   : frag number, frag range ("1-2"), or frag list ("[3,5,8]") (base 0)
STARTRES  : number of starting residue (base 1)
STOPRES   : number of stopping residue (base 1)
PARENTS   : list of parent numbers, i.e., "2" or "3,5" 
            optionally add ":#" to indicate which snapshot ("5:1").
RESLIST   : list of residue numbers, "1", "1-4", "1,2", or "1-3,5,7-9" 
NRES      : number of residues to grow on each end
MAXCHILD  : maximum number of children to grow from a parent fragment
STAGE     : fragment stage, i.e. "grow" or "assem"
REXTIME   : rex simulation time in ps, or "auto" to do automatic calculation
ANALTIME  : time in ps used at the end of rex for analysis, or "auto"
SWAPTIME  : initial time in ps to spread restraints across replicas
RESTLIST  : list of restraints; either like "1 8 2 9" or "[(1,8),(2,9)]"
VARNAME   : name of the variable to modify
VALUE     : value of the variable to set
SEEDLEN   : number of residues in seeds
SEEDGAP   : number of residues spacing seeds
"""

Help = """ZAM COMMAND REFERENCE
Capitalized tokens are parameters.  Brackets indicate optional.

zam.py initializerun SEQ
  Starts a zam run by building a database file, restart.dat.
  SEQ = a string of one-letter amino acid codes, UNCAPPED

zam.py restartrun
  Prepares the zam simulation directory for job submission.
  
zam.py show [FRAGNUM]
  Presents a data table summarizing fragments.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
            Default for no FRAGNUM is all fragments.
  
zam.py showglobal
  Shows global variables in the zam run that apply to all fragments.
  
zam.py showall [FRAGNUM]
  Shows all variables in detailed form for the fragments indicated.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
            Default for no FRAGNUM is all fragments.
            
zam.py showanal [FRAGNUM]
  Shows select analysis variables for the fragments indicated.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
            Default for no FRAGNUM is all fragments.
            
zam.py showrest [FRAGNUM]
  Shows the lists of prior (during-REX), initial (optional initial stages
  in REX), proposed (bestn), and new (applied proposed) restraints.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
            Default for no FRAGNUM is all fragments.
            
zam.py shownote [FRAGNUM]
  Shows any text note given by the user.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
            Default for no FRAGNUM is all fragments.

zam.py decidenext [NFRAG] [byfrag|bysnap]
  Makes next set of fragments for simulation.
  NFRAG = number of new fragments to create
  byfrag = use this token to assemble by fragments
  bysnap = use this token to assemble by individual fragment snapshots

zam.py clearglobalrest
  Clears the list of global restraints, which are automatically applied
  to newly created fragments (but not old).
  
zam.py addglobalrest RESTLIST
  Adds restraints to the list of global restraints, which are automatically
  applied to newly created fragments (but not old).
  RESTLIST = list of restraints; either like "1 8 2 9" or "[(1,8),(2,9)]"

zam.py delglobalrest RESTLIST
  Removes restraints from the list of global restraints, which are
  automatically applied to newly created fragments (but not old).

zam.py addglobalss RESLIST SS
  Adds global secondary structure restraints, automatically applied to newly
  created fragments.
  RERESLIST = list of residue numbers, "1", "1-4", "[1,2]", or "[1-3,5,7-9]"
  SS = secondary structure tag, such as "H" or "E" (or "-" for none)

zam.py clearglobalss
  Removes all global secondary structure restraints.
  
zam.py edit VARNAME VALUE
  Generic command for editing any zam global variable.
  VARNAME = string name of the variable
  VALUE = value to give variable
  
zam.py newfrag STARTRES STOPRES [PARENTS [REXTIME]]
  Creates a new fragment.
  STARTRES = number of starting residue
  STOPRES = number of stopping residue
  PARENTS = single parent number ("5"), range ("5-8"), or list ("[5,8,9]").
            For number or list, add number to indicate which snapshot ("5:1").
            Use "[]" for no parents.
  REXTIME = override REX simulation time in ps
  
zam.py growfrag FRAGNUM [NUMRES | NNUMRES CNUMRES]
  Creates a new fragment grown from an existing one by specified number
  of residues on each end.  Automatically checks N and C terminal bounds.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  NUMRES = number of residues to grow on each end
  NNUMRES = number of residues to grow on the n-terminal
  CNUMRES = number of residues to grow on the c-terminal

zam.py duplicatefrag FRAGNUM
  Creates a new fragment with the same residue range, parents,
  restraints, and initial seeds as FRAGNUM.  Stage is reset to
  the default.  Does NOT copy rex time or analysis time.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  
zam.py makeseedfrags SEEDLEN SEEDGAP [STARTRES STOPRES]
  Makes seed fragments evenly distributed across the sequence.
  SEEDLEN = length of seeds in number of residues
  SEEDGAP = gap in number of residues between seeds
  STARTRES = optional starting residue number
  STOPRES = optional stopping residue number
  
zam.py makeconcensusfrags STARTRES STOPRES PARENT [MAXCHILD|RESTLIST]
  Makes a set of indentical fragments from a common parent which each have
  a different restraint.
  STARTRES = number of starting residue
  STOPRES = number of stopping residue
  PARENT = single parent number (cannot be multiple);
           add number to indicate which snapshot ("5:1")  
  MAXCHILD = maximum number of child fragments to generate using restraints
             contained in parent proposed restraints
  RESTLIST = a specific list of restraints to use, one per child fragment;
             list can be either like "1 8 2 9" or "[(1,8),(2,9)]"
             
zam.py combineconcensusfrags PARENTS [REXTIME]
  Makes a new fragment that is a concensus run of a number of parents.
  Only restraints that are common to all parents will be applied.
  PARENTS = single parent number ("5"), range ("5-8"), or list ("[5,8,9]").
            For number or list, add number to indicate which snapshot ("5:1").
  REXTIME = override REX simulation time in ps
  
zam.py fragparents FRAGNUM PARENTS
  Sets the parents of fragments, which are used for initial config making.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  PARENTS = single parent number ("5"), range ("5-8"), or list ("[5,8,9]")
            For number or list, add number to indicate which snapshot ("5:1").
  
zam.py delfrag FRAGNUM
  Deletes fragments, including all associated data files.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  
zam.py delchildren FRAGNUM
  Deletes the child fragments of specified fragments, including all
  associated datafiles.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")

zam.py estimatereplicas FRAGNUM 
  Estimates the total number of replicas in fragments.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")

zam.py fragreplicas FRAGNUM NREPLICA
  Overrides the REX time for fragments. Do not change once a run has started!
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  NREPLICA = number of replicas or "auto" to have automatically determined
  
zam.py fragrextime FRAGNUM REXTIME
  Overrides the REX time for fragments.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  REXTIME = REX time in ps or "auto" to have automatically determined
  
zam.py fraganaltime FRAGNUM ANALTIME [SKIPTIME]
  Overrides the analysis time for fragments, i.e. the ultimate number of ps 
  used in the analysis routines.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  ANALTIME = analysis time in ps or "auto" to have automatically determined
  SKIPTIME = if supplied, will skip a preset number of ps, rather than just
             use those at the end of the trajectory

zam.py fraganalopt FRAGNUM {punchon|punchoff}
  Sets analysis routines to run.  The tokens "punchon" and "punchoff" will
  determine whether time-consuming punch analysis is performed.  If punch
  is not run, the bestn contacts will be taken from contact probability.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  
zam.py restartfrag FRAGNUM [REXTIME]
  Restarts fragments at the REX/analysis stages.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  REXTIME = new REX time in ps
  
zam.py reinitfrag FRAGNUM [REXTIME]
  Restarts fragments at configuration making and deletes all current REX
  and analysis data, i.e., starts a fragment anew.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  REXTIME = REX time in ps
  
zam.py fragstage FRAGNUM STAGE
  Sets the stage of fragments.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  STAGE = string like "grow" or "done" or "old"
  
zam.py clearrest FRAGNUM
  Removes all restraints from fragments.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  
zam.py addrest FRAGNUM RESTLIST
  Adds restraints to fragments.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  RESTLIST = list of restraints; either like "1 8 2 9" or "[(1,8),(2,9)]"
  
zam.py delrest FRAGNUM RESTLIST
  Removes restraints from fragments.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  RESTLIST = list of restraints; either like "1 8 2 9" or "[(1,8),(2,9)]"
  
zam.py copyrest DESTFRAGNUM SOURCEFRAGNUM
  Copies all restraints from one fragment and ADDS them to another.
  DESTFRAGNUM = single fragment ("5"), range ("5-8"), or list ("[5,8,9]")
  SOURCEFRAGNUM = single fragment number ("6")

zam.py addss FRAGNUM RESLIST SS
  Adds secondary structure restraints to a fragment.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  RESLIST = list of residue numbers, "1", "1-4", "[1,2]", or "[1-3,5,7-9]"
  SS = secondary structure tag, such as "H" or "E" (or "-" for none)

zam.py clearss FRAGNUM
  Clears secondary structure restraints in a fragment.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")

zam.py rigidifyfrag FRAGNUM RESLIST
  Adds backbone-rigidifying restraints to segments in fragments.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  RESLIST = list of residue numbers, "1", "1-4", "[1,2]", or "[1-3,5,7-9]"
  
zam.py unrigidifyfrag FRAGNUM
  Removes all rigidifying restraints from fragments.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")

zam.py anchorfrag FRAGNUM BBRESLIST AARESLIST
  Adds anchoring restraints to segments in fragments. Put ref.pdb in conf/.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  BBRESLIST = list of backbone residue numbers, "1", "1-4", "[1,2]", or "[1-3,5,7-9]"
  AARESLIST = list of all-atom residue numbers, "1", "1-4", "[1,2]", or "[1-3,5,7-9]"
  
zam.py unanchorfrag FRAGNUM
  Removes all anchoring restraints from fragments.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")  
  
zam.py addswaprest FRAGNUM SWAPTIME SWAPRESTLIST
  Adds initial restraints that are distributed across replicas for a short
  initial time during REX and then turned off.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  SWAPTIME = initial time to run these restraints, in ps; specify a value
             of 0 to bypass swapped restraints
  SWAPRESTLIST = list of restraints; either like "1 8 2 9" or "[(1,8),(2,9)]"

zam.py clearswaprest FRAGNUM 
  Clears initial swapped restraints that are distributed across replicas .
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")

zam.py repelsaltbridges FRAGNUM {on|off}
  Turns salt bridge repulsion on or off.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")  

zam.py scalerest FRAGNUM {on|off}
  Turns restraint scaling on or off.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")  

zam.py addumbrellarest FRAGNUM RESLIST DIST1 DIST2 ...
  Adds umbrella sampling to a fragment.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  RESTLIST = list of residue pairs; either like "1 8 2 9" or "[(1,8),(2,9)]"

zam.py umbrelladistances FRAGNUM DIST1 DIST2 ...
  Adds umbrella sampling to a fragment.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  DIST* = distances for different umbrellas (must be 2 or more)
  
zam.py clearumbrellarest FRAGNUM
  Removes all umbrella sampling in a fragment.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  
zam.py fragpairs FRAGNUM {all | hydrophobic}
  Specifies whether analysis routines should consider all or just
  hydrophobic residue-residue pairs.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")

zam.py fragcaps FRAGNUM {none | both | ncap | ccap}
  Sets capping options.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  
zam.py fragnote FRAGNUM [NOTETEXT]
  Adds a user-specified note to fragments.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  NOTETEXT = text of the note; the note is erased if none is specified
  
zam.py editfrag FRAGNUM VARNAME VALUE
  Generic command for editing any variable of fragments.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  VARNAME = name of the fragment variable
  VALUE = value to give the variable
  
zam.py compressdata FRAGNUM FILENAME [anal|prep] [database]
  Compresses fragment data into a gzipped tarball.  The optional token "anal"
  will compress just the analysis data.  The optional token "prep" will
  compress just data needed to run the specified fragments, including parent
  configurations and the fragment database (restart.dat).  The optional
  token "database" will include the fragment database (restart.dat).
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  
zam.py summarize FRAGNUM KEYNAME
  Runs the summarize.py html generator on fragments.  
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  KEYNAME = directory name for the results, e.g., "T0283/8mers"
  
zam.py deldata FRAGNUM [workonly|dataonly]
  Deletes the data and work directories for a fragment.  The optional tokens
  "workonly"/"dataonly" will only remove the work/data directories.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")

zam.py savescores FRAGNUM FILEPREFIX
  Saves current decision-making score database to a text file.
  FRAGNUM = single fragment number ("5"), range ("5-8"), or list ("[5,8,9]")
  FILEPREFIX = prefix of files: FILEPREFIX.txt and FILEPREFIX.html
"""


FragShowDat = ["RexTime", "Stage", "Note"]
FragAnalDat = ["MesoEntropy", "MesoGroundPop", "StructScore", "Productive",
               "PunchStabTot", "PunchCoopTot", "PriorRest", "PropRest"]
FragHideDat = ["Seq","Parents","PriorRest","PropRest","NewRest","SwapRest",
               "UmbrRest","ResRigid","ResAnchorBB","ResAnchorAA"]
ZDHideDat = ["Frags", "Rest", "Seq"]
MaxFragSeqLen = 20
StageDflt = "grow"
RexRestartInd = 1
ParentsByIndex = False

import sys
#check for instructions
if len(sys.argv) == 1:
  print Usage
  sys.exit()
else:
  print ""

import os, sys, cPickle, shutil, copy, glob
import zamdata, zamchoice, zamcontact
import sequence, scripttools


def SortUnique(l):
  l = sorted(l)
  return [x for (i,x) in enumerate(l) if not x in l[i+1:]]

def GetFragList(Arg):
  if Arg.lower() == "none":
    return []
  elif Arg == "*" or Arg.lower() == "all":
    return [f for f in zd.Frags]
  else:
    FragList = []
    Arg = Arg.replace("[","").replace("]","")
    for l in Arg.split(","):
      if "-" in l:
        a, b = [int(x) for x in l.split("-")]
        FragList.extend(range(a, b+1))
      else:
        a = int(l)
        FragList.append(a)
    FragList = [i for i in FragList if i >= 0 and i < len(zd.Frags)]
    FragList = SortUnique(FragList)
    FragList = [zd.Frags[i] for i in FragList]
    return FragList
  
def GetParentList(Arg):
  def GetParent(val):
    if ":" in val:
      a,b = val.split(":")[:2]
      if b == "*":
        return int(a), -1
      else:
        return int(a), int(b)
    else:
      return int(val), -1
  if Arg.lower() == "none" or Arg == "[]":
    return []
  else:
    ParentList = []
    Arg = Arg.replace("[","").replace("]","")
    for l in Arg.split(","):
      if "-" in l:
        (a, b), (c, d) = [GetParent(x) for x in l.split("-")]
        ParentList.extend([(i,b) for i in range(a,c)])
        ParentList.append((c,d))
      else:
        a, b = GetParent(l)
        ParentList.append((a,b))
    ParentList = [(a,b) for (a,b) in ParentList if a >= 0 and a < len(zd.Frags)]
    ParentList = SortUnique(ParentList)
    ParentList = [(zd.Frags[a], b) for (a,b) in ParentList]
    return ParentList

def GetRestList(Args):
  l = []
  a,b = None,None
  for itm in Args:
    v = eval(itm)
    if type(v) is list:
      l.extend(v)
    else:
      if a is None:
        a = v
      else:
        b = v
    if not b is None:
      l.append((a,b))
      a,b = None,None
  if not a is None and b is None:
    raise ValueError, "Could not match all restraint pairs."
  l = [(a-1,b-1) for (a,b) in l]
  l = [(min(a,b),max(a,b)) for (a,b) in l]
  if len(l) == 0:
    print "No restraints specified."
    sys.exit()
  return l

def GetResList(Arg):
  global zd
  if Arg.lower() == "none":
    return []
  elif Arg == "*" or Arg.lower() == "all":
    return range(len(zd.Seq))
  else:
    ResList = []
    Arg = Arg.replace("[","").replace("]","")
    for l in Arg.split(","):
      if "-" in l:
        a, b = [int(x) for x in l.split("-")]
        ResList.extend(range(a-1, b))
      else:
        a = int(l)
        ResList.append(a - 1)
  return ResList

def SortUnique(l):
  l = sorted(l)
  return [x for (i,x) in enumerate(l) if not x in l[i+1:]]

def RestListStr(RestList):
  global zd
  return ", ".join(["%s%d-%s%d" % (zd.Seq[a], a+1, zd.Seq[b], b+1)
        for (a,b) in RestList])

def ResListStr(ResList):
  global zd
  return ", ".join(["%s%d" % (zd.Seq[a], a+1) for a in ResList])

def ResListStrShort(ResList):
  global zd
  ResList = SortUnique(ResList)
  l = []
  for i, a in enumerate(ResList):
    if a == ResList[0] or a == ResList[-1]:
      l.append("%s%d" % (zd.Seq[a], a+1))
    elif a == ResList[i-1] + 1 and a == ResList[i+1] - 1:
      if not l[-1] == "-": l.append("-")
    else:
      l.append("%s%d" % (zd.Seq[a], a+1))
  l = ", ".join(l)
  l = l.replace(", -, ", "-")
  return l

def ParentStr(Parent):
  global zd
  Frag, Ind = Parent
  if ParentsByIndex:
    if Ind < 0:
      return "%d:*" % zd.Frags.index(Frag)
    else:
      return "%d:%d" % (zd.Frags.index(Frag), Ind)
  else:
    if Ind < 0:
      return "%s:*" % Frag.Name()
    else:
      return "%s:%d" % (Frag.Name(), Ind)

def ParentsStr(Parents):
  global zd
  return "[" + ",".join([ParentStr(x) for x in Parents]) + "]"

def ShowList(List, BlockLen = 33, NBlock = 2, Indent1 = "  ", Indent2 = "  "):
  #make a new list that fits the block sizes
  l = []
  n = BlockLen - len(Indent2)
  for itm in List:
    l.append(itm[:BlockLen].ljust(BlockLen))
    itm = itm[BlockLen:]
    while len(itm) > 0:
      l.append(Indent2 + itm[:n].ljust(n))
      itm = itm[n:]
  #find the number of list items per column
  n = (len(l) + NBlock - 1) / NBlock
  #make the number of list items an even multiple of NBlock
  l.extend([""]*(n*NBlock - len(l)))
  for i in range(n):
    print Indent1 + " ".join(l[i::n])

def ShowTable(Table, BlockLen = None, ColBorder = " "):
  NCol = len(Table[0])
  NRow = len(Table)
  if BlockLen is None:
    ColW = [max([len(Table[i][j]) for i in range(0,NRow)]) for j in range(0,NCol)]
  else:
    ColW = [BlockLen]*NCol
  for Row in Table:
    for i in range(0,NCol):
      Cell = Row[i].replace("\n","")
      if Cell == "-":
        Cell = "-"*(ColW[i] + len(ColBorder))
      else:
        Cell = Cell.ljust(ColW[i])[:ColW[i]] + ColBorder
      sys.stdout.write(Cell)
    sys.stdout.write("\n")

def PriorRestToList(Frag):
  if len(Frag.PriorRest) > 0:
    return ["%s%d-%s%d" % (Frag.Seq[a-Frag.StartRes], a+1, Frag.Seq[b-Frag.StartRes], b+1)
        for (a,b) in Frag.PriorRest]
  else:
    return []

def PropRestToList(Frag):
  if len(Frag.PropRest) > 0:  
    return ["%s%d-%s%d(%.2f)" % (Frag.Seq[a-Frag.StartRes], a+1, Frag.Seq[b-Frag.StartRes], b+1, s)
        for (s,a,b) in sorted(Frag.PropRest, reverse=True)]
  else:
    return []

def SwapRestToList(Frag):
  if len(Frag.SwapRest) > 0:
    return ["%s%d-%s%d" % (Frag.Seq[a-Frag.StartRes], a+1, Frag.Seq[b-Frag.StartRes], b+1)
        for (a,b) in Frag.SwapRest]
  else:
    return []

def UmbrRestToList(Frag):
  if len(Frag.UmbrRest) > 0:
    return ["%s%d-%s%d" % (Frag.Seq[a-Frag.StartRes], a+1, Frag.Seq[b-Frag.StartRes], b+1)
        for (a,b) in Frag.UmbrRest]
  else:
    return []   


def ShowRest():
  #show the master restraints
  if len(zd.Rest) > 0:
    RestLabels = ["%s%d-%s%d" % (zd.Seq[a], a+1, zd.Seq[b], b+1) for (a,b) in sorted(zd.Rest)]
    print "Master restraints:"
    ShowList(RestLabels, BlockLen = 14, NBlock = 4)
  else:
    print "Master restraints: None"

def ShowZDData(DataKeys = None):
  #show all the properties
  if DataKeys is None:
    DataKeys = zd.__dict__.keys()
  print "Summary of ZAM data:"
  print "--------------------"
  print "Seq: %s" % " ".join(sequence.SeqToList(zd.Seq))
  print "Seq: %s" % sequence.SeqToAA1(zd.Seq)
  print "Len: %d residues" % len(zd.Seq)
  print "Number of frags: %d" % len(zd.Frags)
  ShowRest()
  for k in DataKeys:
    if k in ZDHideDat: continue
    v = zd.__dict__[k]
    print str(k) + ": ", v
  print "" 

def FragSeq(Frag):
  #trim the sequence
  Seq = sequence.SeqToAA1(Frag.Seq)
  if len(Seq) > MaxFragSeqLen:
    Seq = Seq[:MaxFragSeqLen-1] + "+"
  return Seq

def ShowFragData(DataKeys = None, FragList = None, BlockLen = 33, NBlock = 2):
  #show all the properties
  if FragList is None:
    FragList = zd.Frags
  for f in FragList:
    List = []
    print "Fragment %d (%s):" % (zd.Frags.index(f), f.Name())
    List.append("Seq: %s" % FragSeq(f)) 
    List.append("Parents: " +  ParentsStr(f.Parents))
    List.append("ResRigid: " + ResListStrShort(f.ResRigid))
    List.append("ResAnchorBB: " + ResListStrShort(f.ResAnchorBB))
    List.append("ResAnchorAA: " + ResListStrShort(f.ResAnchorAA))
    if DataKeys is None:
      Keys = f.__dict__.keys()
      Keys.sort()
    else:
      Keys = DataKeys
    for k in Keys:
      if k in FragHideDat: continue
      if k in f.__dict__:
        v = f.__dict__[k]
        List.append("%s: " % k + str(v))
      else:
        List.append("")
    ShowList(List, BlockLen = BlockLen, NBlock = NBlock)
    print ""
  print "" 

def ShowFragDataTable(DataKeys = None, FragList = None, ShowSeq = True):
  if FragList is None: FragList = zd.Frags
  if len(zd.Frags) == 0: return
  if DataKeys is None:
    DataKeys = [k for k in zd.Frags[0].__dict__.keys() if not k in FragHideDat]
    DataKeys.sort()
  Table = [["Frag","Name", "Seq", "Parents"] + DataKeys]
  for f in FragList:
    ind = zd.Frags.index(f)
    List = [str(ind), f.Name()]
    if ShowSeq: List.append(FragSeq(f))
    List.append(ParentsStr(f.Parents))
    for k in DataKeys:
      if k == "PriorRest":
        List.append(",".join(PriorRestToList(f)))
      elif k == "PropRest":
        List.append(",".join(PropRestToList(f)))
      elif k == "Note":
        if "Note" in f.__dict__:
          List.append(f.Note[:20])
        else:
          List.append("")
      else:
        if k in f.__dict__:
          v = f.__dict__[k]
          if type(v) is int:
            List.append("%d" % v)
          elif type(v) is str:
            List.append(v)
          else:
            try:
              List.append("%.2f" % v)
            except:
              List.append(str(v))
        else:
          List.append("")
    Table.append(List)
  ShowTable(Table)
  print ""   

def ShowFragRest(FragList = None):
  #show all the properties
  if FragList is None:
    FragList = zd.Frags
  Table = [["Frag", "Name", "PriorRest", "SwapRest", "PropRest", "UmbrRest"]]  
  RowBor = ["-"]*6
  for f in FragList:
    Cols = [[]]*6
    Cols[0] = ["%d" % zd.Frags.index(f)]
    Cols[1] = [f.Name()]
    Cols[2] = PriorRestToList(f)
    Cols[3] = SwapRestToList(f)
    Cols[4] = PropRestToList(f)
    Cols[5] = UmbrRestToList(f)
    MaxLen = max([len(x) for x in Cols])
    Cols = [c + [""]*(MaxLen - len(c)) for c in Cols]
    Table.append(RowBor)
    for i in range(0, MaxLen):
      Table.append([c[i] for c in Cols])
  ShowTable(Table)
  print ""

def ShowFragPropRest(FragList = None):
  #show all the properties
  if FragList is None:
    FragList = zd.Frags
  Table = [["Frag", "Name", "PropRest", "", ""]]
  RowBor = ["-"]*5
  for f in FragList:
    Cols = [[]]*5
    Cols[0] = ["%d" % zd.Frags.index(f)]
    Cols[1] = [f.Name()]
    l = PropRestToList(f)
    n = len(l)
    dn = n/3
    if  n%3 > 0: dn += 1
    Cols[2] = l[:dn]
    Cols[3] = l[dn:2*dn]
    Cols[4] = l[2*dn:]
    MaxLen = max([len(x) for x in Cols])
    Cols = [c + [""]*(MaxLen - len(c)) for c in Cols]
    Table.append(RowBor)
    for i in range(0, MaxLen):
      Table.append([c[i] for c in Cols])
  ShowTable(Table)
  print ""  

def ShowNote(FragList = None):
  #show the notes
  if FragList is None: FragList = zd.Frags
  for f in FragList:
    s = "Fragment %d (%s): " % (zd.Frags.index(f), f.Name())
    if "Note" in f.__dict__ and len(f.Note) > 0:
      s += f.Note 
    else:
      s += "none"
    print s + "\n"

def ShowWarning(FragList):
  print "WARNING: This will delete data from fragments", \
        [zd.Frags.index(f) for f in FragList]
  ret = raw_input("Are you sure?  Enter 'y': ")
  if ret.lower() == 'y':
    return True
  else:
    print "Action cancelled."
    return False


def Load():
  global dat, zd
  #load file
  dat = cPickle.load(file("restart.dat","r"))
  zd = dat[0][1]
  #fix for old data
  zd.UpdateVer()
  #check for zamsetup
  if os.path.isfile("zamsetup.py"): execfile("zamsetup.py")

def Save():
  global dat
  cPickle.dump(dat, file("restart.dat","w"))




Cmd = sys.argv[1].lower()

#check for help file
if Cmd == "help":
  print Help
  sys.exit()
#next check if this is an initialize
elif Cmd == "initializerun":
  #check for zamsetup
  zd = zamdata.ZamDataClass()
  if os.path.isfile("zamsetup.py"): execfile("zamsetup.py")
  #start a new zam run by making a restart file
  Seq = " ".join(sys.argv[2:])
  Seq = sequence.SeqToList(Seq)
  zd = zamdata.ZamDataClass(Seq=Seq)
  dat = [('zd', zd)]
  Save()
  print "Initialized run with sequence: %s" % sequence.SeqToAA1(Seq)
  sys.exit()


#check for the file
if not os.path.isfile("restart.dat"):
  print "Could not find restart.dat file"
  sys.exit()

Load()

if Cmd == "showglobal":
  ShowZDData()

elif Cmd == "show":
  #show the properties
  if len(sys.argv) > 2:
    FragList = GetFragList(sys.argv[2])
    ShowFragDataTable(DataKeys = FragShowDat, FragList = FragList)
  else:
    ShowZDData()
    ShowFragDataTable(DataKeys = FragShowDat)

elif Cmd == "show2":
  #show the properties
  if len(sys.argv) > 2:
    FragList = GetFragList(sys.argv[2])
    ShowFragData(DataKeys = FragShowDat, FragList = FragList)
  else:
    ShowZDData()
    ShowFragData(DataKeys = FragShowDat)

elif Cmd == "showall":
  #show all the properties
  if len(sys.argv) > 2:
    FragList = GetFragList(sys.argv[2])
    ShowFragData(FragList = FragList, BlockLen = 75, NBlock = 1)
  else:
    ShowZDData()
    ShowFragData(BlockLen = 75, NBlock = 1)

elif Cmd == "showanal2":
  #show the analysis results
  if len(sys.argv) > 2:
    FragList = GetFragList(sys.argv[2])
    ShowFragData(DataKeys = FragAnalDat, FragList = FragList)
  else:
    ShowFragData(DataKeys = FragAnalDat)

elif Cmd == "showanal":
  #show the analysis results
  if len(sys.argv) > 2:
    FragList = GetFragList(sys.argv[2])
    ShowFragDataTable(DataKeys = FragAnalDat, FragList = FragList)
  else:
    ShowFragDataTable(DataKeys = FragAnalDat)

elif Cmd == "showrest":
  #show fragment restraints, minus proposed ones
  if len(sys.argv) > 2:
    FragList = GetFragList(sys.argv[2])
    ShowFragRest(FragList = FragList)
  else:
    ShowRest()
    ShowFragRest()

elif Cmd == "showproprest":
  #show proposed fragment restraints
  if len(sys.argv) > 2:
    FragList = GetFragList(sys.argv[2])
    ShowFragPropRest(FragList = FragList)
  else:
    ShowRest()
    ShowFragPropRest()    

elif Cmd == "shownote":
  #show proposed fragment restraints
  if len(sys.argv) > 2:
    FragList = GetFragList(sys.argv[2])
  else:
    FragList = zd.Frags
  ShowNote(FragList)  


elif Cmd == "decidenext":
  #decide what to do next
  AssembleByFrag = None
  if len(sys.argv) > 2:
    s = " ".join(sys.argv[2:]).lower()
    if "byfrag" in s:
      AssembleByFrag = True
      s = s.replace("byfrag","")
    if "bysnap" in s:
      AssembleByFrag = False
      s = s.replace("bysnap","")
    s = s.strip()
    if len(s) > 0:
      zamchoice.MaxFrags = int(s)
  if len(zd.GetRunFrags()) > 0:
    print "Running fragments already exist."
  else:
    print "Loading the contact database..."
    cd = zamcontact.ContactDatabaseClass(zd)
    print "Loading the structure database..."
    sd = zamdata.SnapDatabaseClass(zd, cd.Scores)
    print "Deciding next fragments..."
    Stat = zamchoice.GetNextActions(zd, cd, sd, AssembleByFrag = AssembleByFrag)
    if Stat == "stuck":
      print "ZAM simulation is stuck."
    elif Stat == "done":
      print "ZAM simulation is complete." 
    elif Stat == "stop":
      print "ZAM simulation is stopped."
    elif Stat == "run":
      for f in zd.GetRunFrags():
        f.MakePath(AllPaths = True)
        print "Created new fragment %s as fragment %d" % (f.Name(), zd.Frags.index(f))
        if len(f.Parents) > 0: print "Parent fragments are %s." % ParentsStr(f.Parents)
        if len(f.PriorRest) > 0:
          print "Current restraints copied from parents are:"
          ShowList(PriorRestToList(f), BlockLen = 20, NBlock = 3)
        print ""
    Save()

elif Cmd == "edit":
  #edit a zd property
  zd.__dict__[sys.argv[2]] = eval(" ".join(sys.argv[3:]))
  Save()

elif Cmd == "editfrag":
  #edit a fragment property
  FragList = GetFragList(sys.argv[2])
  Key = sys.argv[3]
  Val = eval(" ".join(sys.argv[4:]))
  for f in FragList:
    f.__dict__[Key] = Val
  Save()

elif Cmd == "delfrag":
  #delete a fragment
  FragList = GetFragList(sys.argv[2])
  if ShowWarning(FragList):
    for Frag in FragList:
      #get rid of the folder
      Frag.DelPath()
      #get rid of the references from children & parents
      Frag.Unlink()
      #delete the fragment
      ind = zd.Frags.index(Frag)
      zd.DelFrag(ind)
      print "Deleted fragment %s" % Frag.Name()
    Save()

elif Cmd == "delchildren":
  #delete all the children of a fragment
  FragList = GetFragList(sys.argv[2])
  ChildFrags = []
  for f in FragList: ChildFrags.extend(Frag.GetAllChildren())
  if ShowWarning(ChildFrags):
    for Frag in FragList:
      for Child in Frag.GetAllChildren():
        Child.Unlink()
        Child.DelPath()
        if Child in zd.Frags:
          ind = zd.Frags.index(Child)
          zd.DelFrag(ind)
          print "Deleted fragment child %s of %s" % (Child.Name(), Frag.Name())
    Save()  

elif Cmd == "newfrag":
  #make a new fragment
  StartRes = max(int(sys.argv[2]), 1) - 1
  StopRes = min(int(sys.argv[3]), len(zd.Seq)) - 1
  Stage = StageDflt
  if len(sys.argv) > 4:
    Parents = GetParentList(sys.argv[4])
  else:
    Parents = []
  f = zd.AddFrag(StartRes, StopRes, Parents, Stage)
  if StartRes == 0 and StopRes == len(zd.Seq) - 1:
    ret = raw_input("Fragment is full sequence. Add caps to existing sequence (y/n)? ")
    ret = ret[0].lower()
    if ret == 'n': f.CapN, f.CapC = False, False
  f.MakePath(AllPaths = True)
  if len(sys.argv) > 5:
    f.CalcRexTime = False
    f.RexTime = float(sys.argv[5])
  print "Created new fragment %s as fragment %d" % (f.Name(), zd.Frags.index(f))
  if len(Parents) > 0:
    print "Parent fragments are %s." % ParentsStr(Parents)
  if len(f.PriorRest) > 0:
    print "Current restraints copied from parents are:"
    ShowList(PriorRestToList(f), BlockLen = 20, NBlock = 3)
  print ""
  Save()

elif Cmd == "growfrag":
  #grow into a new fragment
  FragList = GetFragList(sys.argv[2])
  if len(sys.argv) > 4:
    NRes, CRes = int(sys.argv[3]), int(sys.argv[4])
  else:
    NRes = int(sys.argv[3])
    CRes = NRes
  for Frag in FragList:
    StartRes = max(Frag.StartRes - NRes, 0)
    StopRes = min(Frag.StopRes + CRes, len(zd.Seq) - 1)
    Stage = StageDflt
    Parents = [(Frag, -1)]
    f = zd.AddFrag(StartRes, StopRes, Parents, Stage)
    if StartRes == 0 and StopRes == len(zd.Seq) - 1:
      ret = raw_input("Fragment grown from %s is full sequence. Decap (y/n)? " % Frag.Name())
      ret = ret[0].lower()
      if ret == 'y': f.CapN, f.CapC = False, False
    f.MakePath(AllPaths = True)
    print "Created new fragment %s as fragment %d" % (f.Name(), zd.Frags.index(f))
    if len(Parents) > 0:
      print "Parent fragments are %s." % ParentsStr(Parents)
    if len(f.PriorRest) > 0:
      print "Current restraints copied from parents are:"
      ShowList(PriorRestToList(f), BlockLen = 20, NBlock = 3)
    print ""
  Save()  

elif Cmd == "duplicatefrag":
  #duplicates fragments
  FragList = GetFragList(sys.argv[2])
  for Frag in FragList:
    StartRes = Frag.StartRes
    StopRes = Frag.StopRes
    Stage = StageDflt
    Parents = copy.copy(Frag.Parents)
    f = zd.AddFrag(StartRes, StopRes, Parents, Stage)
    f.MakePath(AllPaths = True)
    #copy restraint params
    f.PriorRest = copy.deepcopy(Frag.PriorRest)
    f.SwapRest = copy.deepcopy(Frag.SwapRest)
    f.RexTimeSwapRest = Frag.RexTimeSwapRest
    #copy conf
    SeedFiles = glob.glob(os.path.join(Frag.BasePath, "conf/seed*.pdb"))
    for Seed in SeedFiles:
      NewSeed = os.path.join(f.BasePath, "conf/%s" % os.path.basename(Seed))
      shutil.copy(Seed, NewSeed)
    #copy nrep
    f.NRep = Frag.NRep
    print "Created new fragment %s as fragment %d" % (f.Name(), zd.Frags.index(f))
    if len(Parents) > 0:
      print "Parent fragments are %s." % ParentsStr(Parents)
    if len(f.PriorRest) > 0:
      print "Restraints are:"
      ShowList(PriorRestToList(f), BlockLen = 20, NBlock = 3)
    print ""
  Save()

elif Cmd == "fragparents":
  #change the parents of a fragment
  FragList = GetFragList(sys.argv[2])
  Parents = GetParentList(sys.argv[3])
  for Frag in FragList:
    Frag.Parents = Parents
    print "Set fragment %s to have parents %s" % (Frag.Name(), ParentsStr(Parents))
  Save()

elif Cmd == "makeseedfrags":
  #makes new seeds
  StartRes = 0
  StopRes = len(zd.Seq) - 1
  NRes = int(sys.argv[2])
  ResGap = int(sys.argv[3])
  if len(sys.argv) > 4:
    StartRes = max(int(sys.argv[4]) - 1, StartRes)
    StopRes = min(int(sys.argv[5]) - 1, StopRes)
  SegList = [(i, i+NRes-1) for i in range(StartRes, StopRes-NRes+2, ResGap)]
  LastSeg = (StopRes-NRes+1, StopRes)
  if not LastSeg in SegList: SegList.append(LastSeg)
  for (a,b) in SegList:
    f = zd.AddFrag(a, b, Stage = StageDflt)
    f.MakePath(AllPaths = True)
  print "Made %d seeds." % len(SegList)
  Save()

elif Cmd == "makeconcensusfrags":
  #grows a new set of fragments from a common parent, each with a different restraint
  StartRes = max(int(sys.argv[2]), 1) - 1
  StopRes = min(int(sys.argv[3]), len(zd.Seq)) - 1
  Stage = StageDflt
  n = int(sys.argv[4])
  Frag = zd.Frags[n]
  Parents = [(Frag,-1)]
  NChild = None
  Rest = None
  #check the argument and get the restraints
  if len(sys.argv) == 6:
    val = eval(sys.argv[5])
    if type(val) is int:
      NChild = val
    else:
      Rest = GetRestList(sys.argv[5])
  elif len(sys.argv) > 6:
    Rest = GetRestList(sys.argv[5:])
  else:
    Rest = [(a,b) for (s,a,b) in Frag.PropRest]
  if NChild is None:
    NChild = len(Rest)
  #make the new fragments
  for i in range(0, NChild):
    f = zd.AddFrag(StartRes, StopRes, Parents, Stage)
    f.MakePath(AllPaths = True)
    f.AddRest([Rest[i]])
    a,b = Rest[i]
    print "Made new fragment %s from parent %s with restraint %s%d-%s%d" % \
          (f.Name(), Frag.Name(), zd.Seq[a], a+1, zd.Seq[b], b+1)
    print "Current restraints are:"
    ShowList(PriorRestToList(f), BlockLen = 20, NBlock = 3)
  Save()

elif Cmd == "combineconcensusfrags":
  #combine all previous fragments into a concensus
  Parents = GetParentList(sys.argv[2])
  if len(Parents) == 0:
    print "ERROR: No parents specified."
    sys.exit()
  #find the common restraints
  RestList = []
  for f1 in Parents:
    for r in f1.PriorRest:
      InOthers = ([r in f2.PriorRest for f2 in Parents if not f1==f2].count(False) == 0)
      if InOthers and not r in RestList:
        RestList.append(r)
  #get the start and stop res
  StartRes = min([f.StartRes for f in Parents])
  StopRes = max([f.StopRes for f in Parents])
  #setup the new fragment
  Stage = StageDflt
  f = zd.AddFrag(StartRes, StopRes, Parents, Stage)
  f.PriorRest = RestList
  f.MakePath(AllPaths = True)
  if len(sys.argv) > 3:
    f.CalcRexTime = False
    f.RexTime = float(sys.argv[3])
  print "Created new concensus fragment %s as fragment %d" % (f.Name(), zd.Frags.index(f))
  print "Parent fragments are %s." % ParentsStr(Parents)
  if len(f.PriorRest) > 0:
    print "Current restraints copied from parents are:"
    ShowList(PriorRestToList(f), BlockLen = 20, NBlock = 3)
  Save()

elif Cmd == "estimatereplicas":
  #set the number of replicas
  import zamrex
  FragList = GetFragList(sys.argv[2])
  NRep = 0
  for Frag in FragList:
    if Frag.NRep == 0:
      NRep += zamrex.CalcNRep(Frag)
    else:
      NRep += abs(Frag.NRep)
  print "Estimated total number of replicas is %d" % NRep
  
elif Cmd == "fragreplicas":
  #set the number of replicas
  FragList = GetFragList(sys.argv[2])
  for Frag in FragList:
    v = sys.argv[3].lower()
    #issue a warning
    if Frag.NRep > 0:
      print "WARNING for fragment %s: if already simulated with a previous" % Frag.Name()
      print "number of replicas, changing this value will have destructive effects!"
      ret = raw_input("Continue (y/n)? ")
      ret = ret[0].lower()
      if ret == "n": continue      
    if v == "auto":
      Frag.NRep = 0
      print "Set number of replicas for %s to auto" % (Frag.Name(),)
    else:
      Frag.NRep = int(sys.argv[3])
      print "Set number of replicas for %s to %d" % (Frag.Name(), Frag.NRep)
  Save()  

elif Cmd == "fragrextime":
  #set the rex time
  FragList = GetFragList(sys.argv[2])
  for Frag in FragList:
    v = sys.argv[3].lower()
    if v == "auto":
      Frag.CalcRexTime = True
      print "Set rex time for %s to auto" % Frag.Name()
    else:
      Frag.CalcRexTime = False
      Frag.RexTime = float(sys.argv[3])
      print "Set rex time for %s to %.1f ps" % (Frag.Name(), Frag.RexTime)
  Save()

elif Cmd == "fraganaltime":
  #set the anal time
  FragList = GetFragList(sys.argv[2])
  for Frag in FragList:
    v = sys.argv[3].lower()
    if v == "auto":
      Frag.CalcAnalTime = True
      Frag.SkipTime = -1.
      print "Set anal time for %s to auto" % Frag.Name()
    else:
      Frag.CalcAnalTime = False
      Frag.AnalTime = float(sys.argv[3])
      if len(sys.argv) > 4:
        Frag.SkipTime = float(sys.argv[4])
        print "Set anal time for %s to %.1f ps, skipping %.1ps" \
              % (Frag.Name(), Frag.AnalTime, Frag.SkipTime)
      else:
        print "Set anal time for %s to %.1f ps" % (Frag.Name(), Frag.AnalTime)
  Save()

elif Cmd == "fraganalopt":
  #set the anal options
  FragList = GetFragList(sys.argv[2])
  for Frag in FragList:
    Frag.RunPunch = (Frag.RunPunch or ("punchon" in sys.argv)) \
                     and not ("punchoff" in sys.argv)
    print "Set punch flag for fragment %s to %s" % \
          (Frag.Name(), str(Frag.RunPunch))
  Save()

elif Cmd == "restartfrag":
  #restart a fragment
  FragList = GetFragList(sys.argv[2])
  for Frag in FragList:
    #reset the restart index and stage to that for REX
    Frag.RestartInd = RexRestartInd
    Frag.Stage = StageDflt
    #delete the analysis data
    Frag.DelAnalysis()  
    #change the rextime if desired
    if len(sys.argv) > 3:
      Frag.CalcRexTime = False
      Frag.RexTime = float(sys.argv[3])
    print "Restarted fragment %s" % Frag.Name()
  Save()

elif Cmd == "reinitfrag":
  #erases all previous rex data and reinitializes a fragment
  FragList = GetFragList(sys.argv[2])
  if ShowWarning(FragList):
    for Frag in FragList:
      #reset the restart index and stage
      Frag.RestartInd = 0
      if Frag.Stage == 'done': Frag.Stage = StageDflt
      #delete the analysis data
      Frag.DelAnalysis()
      Frag.DelData()
      #reset number of replicas
      Frag.NRep = 0
      #change the rextime if desired
      if len(sys.argv) > 3:
        Frag.RexTime = float(sys.argv[3])
      print "Reinitialized fragment %s" % Frag.Name()
    Save()

elif Cmd == "fragstage":
  #change the stage of a fragment
  FragList = GetFragList(sys.argv[2])
  for Frag in FragList:
    Frag.Stage = sys.argv[3]
    print "Set fragment stage for %s" % Frag.Name()
  Save()

elif Cmd == "clearglobalrest":
  #clear the restraints
  zd.Rest = []
  print "Cleared global restraints."
  Save()

elif Cmd == "addglobalrest":
  #add a new restraint
  RestList = GetRestList(sys.argv[2:])
  for (a,b) in RestList:
    zd.AddRest(a, b, Force = True)
  print "Added global restraints:  %s" % RestListStr(RestList)
  Save()

elif Cmd == "delglobalrest":
  #delete a restraint
  RestList = GetRestList(sys.argv[2:])
  RestMissing = [x for x in RestList if not x in zd.Rest]
  RestFound = [x for x in RestList if x in zd.Rest]
  if len(RestList) > 0:
    zd.Rest = [x for x in zd.Rest if not x in RestFound]
    print "Deleted global restraints: %s" % RestListStr(RestFound)
    Save()
  if len(RestMissing) > 0:
    print "Could not find global restraints: %s" % RestListStr(RestMissing)

elif Cmd == "addglobalss":
  #adds global secondary structure restraints
  ResList = GetResList(sys.argv[2])
  SS = sys.argv[3]
  zd.AddSS(ResList, SS)
  print "Added global secondary structure restraints for %s: %s" \
          % (SS[0], ResListStrShort(ResList))
  Save()

elif Cmd == "clearglobalss":
  #clears global secondary structure restraints
  zd.SS = ""
  print "Cleared global secondary structure restraints." 
  Save()

elif Cmd == "clearfragrest" or Cmd == "clearrest":
  #clear the restraints for a frag
  FragList = GetFragList(sys.argv[2])
  for Frag in FragList:
    Frag.PriorRest = []
    print "Cleared restraints for fragment %s" % Frag.Name()
  Save()

elif Cmd == "addfragrest" or Cmd == "addrest":
  #add a new restraint for a fragment
  FragList = GetFragList(sys.argv[2])
  RestList = GetRestList(sys.argv[3:])
  for Frag in FragList:
    Frag.AddRest(RestList)
    print "Added restraints to fragment %s: %s" \
          % (Frag.Name(), RestListStr(RestList))
  Save()

elif Cmd == "delfragrest" or Cmd == "delrest":
  #delete a restraint from a fragment
  FragList = GetFragList(sys.argv[2])
  RestList = GetRestList(sys.argv[3:])
  for Frag in FragList:
    RestFound = [x for x in RestList if x in Frag.PriorRest]
    RestMissing = [x for x in RestList if not x in Frag.PriorRest]
    Frag.PriorRest = [x for x in Frag.PriorRest if not x in RestFound]
    if len(RestFound) > 0:
      print "Deleted restraints in fragment %s: %s" \
              % (Frag.Name(), RestListStr(RestFound))
    if len(RestMissing) > 0:
      print "Restraints not found in frag %s: %s" \
              % (Frag.Name(), RestListStr(RestMissing)) 
  Save()

elif Cmd == "copyfragrest" or Cmd == "copyrest":
  #copies restraints from a previous fragment
  FragList = GetFragList(sys.argv[2])
  SrcFrag = zd.Frags[int(sys.argv[3])]
  for Frag in FragList:
    Frag.AddRest(SrcFrag.PriorRest)
    print "Copied restraints to fragment %s from fragment %s: %s" \
          % (Frag.Name(), SrcFrag.Name(), RestListStr(SrcFrag.PriorRest))
  Save()

elif Cmd == "addfragss" or Cmd == "addss":
  #adds secondary structure restraints
  FragList = GetFragList(sys.argv[2])
  ResList = GetResList(sys.argv[3])
  SS = sys.argv[4]
  for Frag in FragList:
    Frag.AddSS(ResList, SS)
    print "Added secondary structure restraints for %s to fragment %s: %s" \
          % (SS[0], Frag.Name(), ResListStrShort(ResList))
  Save()

elif Cmd == "clearfragss" or Cmd == "clearss":
  #clears secondary structure restraints
  FragList = GetFragList(sys.argv[2])
  for Frag in FragList:
    Frag.SS = ""
    print "Cleared secondary structure restraints in fragment %s." % Frag.Name()
  Save()

elif Cmd == "rigidifyfrag":
  #rigidify a segment in a fragment
  FragList = GetFragList(sys.argv[2])
  ResList = GetResList(sys.argv[3])
  for Frag in FragList:
    Frag.AutoRigid = False
    r = Frag.GetOverlapRes(ResList)
    Frag.ResRigid = SortUnique(Frag.ResRigid + r)
    print "Rigidified residues in fragment %s: %s" % (Frag.Name(), ResListStrShort(r))
  Save()

elif Cmd == "unrigidifyfrag":
  #rigidify a segment in a fragment
  FragList = GetFragList(sys.argv[2])
  for Frag in FragList:
    Frag.AutoRigid = True
    Frag.ResRigid = []
    print "Unrigidified fragment %s" % Frag.Name()
  Save()

elif Cmd == "anchorfrag":
  #rigidify a segment in a fragment
  FragList = GetFragList(sys.argv[2])
  ResListBB = GetResList(sys.argv[3])
  ResListAA = GetResList(sys.argv[4])
  for Frag in FragList:
    r = Frag.GetOverlapRes(ResListBB)
    Frag.ResAnchorBB = SortUnique(Frag.ResAnchorBB + r)
    if len(r) > 0: print "Anchored backbone res in fragment %s: %s" % (Frag.Name(), ResListStrShort(r))
    r = Frag.GetOverlapRes(ResListAA)
    Frag.ResAnchorAA = SortUnique(Frag.ResAnchorAA + r)
    if len(r) > 0: print "Anchored backbone res in fragment %s: %s" % (Frag.Name(), ResListStrShort(r))
  Save()

elif Cmd == "unanchorfrag":
  #rigidify a segment in a fragment
  FragList = GetFragList(sys.argv[2])
  for Frag in FragList:
    Frag.ResAnchorBB = []
    Frag.ResAnchorAA = []
    print "Unanchored fragment %s" % Frag.Name()
  Save()

elif Cmd == "addswaprest":
  #setup initial restraints
  FragList = GetFragList(sys.argv[2])
  if sys.argv[3] == 'auto':
    CalcSwapRestTime = True
    SwapRestTime = 0.
  else:
    CalcSwapRestTime = False
    SwapRestTime = float(sys.argv[3])
  if sys.argv[4] == 'auto':
    CalcSwapRest = True
    SwapRestList = []
  else:
    CalcSwapRest = False
    SwapRestList = GetRestList(sys.argv[4:])
  for Frag in FragList:
    Frag.CalcSwapRest = CalcSwapRest
    Frag.CalcSwapRestTime = CalcSwapRestTime
    Frag.RexTimeSwapRest = SwapRestTime
    Frag.SwapRest = [(a,b) for (a,b) in SwapRestList
                     if Frag.RestInRange(a,b)]
    print "Set init restraint time for %s to %.1f using restraints: %s" \
          % (Frag.Name(), Frag.RexTimeSwapRest, RestListStr(Frag.SwapRest))
  Save()

elif Cmd == "clearswaprest":
  #clears swapped restraints
  FragList = GetFragList(sys.argv[2])
  for Frag in FragList:
    Frag.CalcSwapRest = False
    Frag.CalcSwapRestTime = False
    Frag.RexTimeSwapRest = 0.
    Frag.SwapRest = []
    print "Cleared swapped restraints for %s" % Frag.Name()
  Save()  

elif Cmd == "repelsaltbridges":
  #turn salt bridge repulsion on or off
  FragList = GetFragList(sys.argv[2])
  if sys.argv[3].lower() in ["on", "true"]:
    RepelSalt = True
    Text = "on"
  elif sys.argv[3].lower() in ["off", "false"]:
    RepelSalt = False
    Text = "off"
  else:
    print "Error: this command requires 'on' or 'off'"
    sys.exit()
  for Frag in FragList:
    Frag.RepelSalt = RepelSalt
    print "Turned salt bridge repulsion to %s in fragment %s" % (Text, Frag.Name())
  Save()

elif Cmd == "scalerest":
  #turn restraint scaling on or off
  FragList = GetFragList(sys.argv[2])
  if sys.argv[3].lower() in ["on", "true"]:
    ScaleRest = True
    Text = "on"
  elif sys.argv[3].lower() in ["off", "false"]:
    ScaleRest = False
    Text = "off"
  else:
    print "Error: this command requires 'on' or 'off'"
    sys.exit()
  for Frag in FragList:
    Frag.ScaleRest = ScaleRest
    print "Turned restraint scaling to %s in fragment %s" % (Text, Frag.Name())
  Save()  

elif Cmd == "addumbrellarest":
  #add a new restraint for a fragment
  FragList = GetFragList(sys.argv[2])
  RestList = GetRestList(sys.argv[3:])
  for Frag in FragList:
    if len(Frag.UmbrDist) == 0: Frag.UmbrDist = [10., 20.]
    for Rest in RestList:
      Rest = sorted(Rest)
      if not Rest in Frag.UmbrRest: Frag.UmbrRest.append(Rest)
    print "Added umbrella restraints to fragment %s: %s" \
          % (Frag.Name(), RestListStr(RestList))
  Save()

elif Cmd == "umbrelladistances":
  #set the umbrella restraint distances
  FragList = GetFragList(sys.argv[2])
  Dist = " ".join(sys.argv[3:])
  Dist = Dist.replace(","," ")
  Dist = [float(x) for x in Dist.split()]
  for Frag in FragList:
    Frag.UmbrDist = Dist
    print "Set umbrella restraint distances in fragment %s" % Frag.Name()
  Save()

elif Cmd == "clearumbrellarest":
  #clear the restraints for a frag
  FragList = GetFragList(sys.argv[2])
  for Frag in FragList:
    Frag.UmbrRest = []
    Frag.UmbrDist = []
    print "Cleared umbrella restraints for fragment %s" % Frag.Name()
  Save()  
  
elif Cmd == "fragpairs":
  #set residue pairs for analysis
  FragList = GetFragList(sys.argv[2])
  AllPairs = (sys.argv[3].lower() == "all")
  for Frag in FragList:
    Frag.AllPairs = AllPairs
    if AllPairs:
      print "Set fragment %s to use all contact pairs." % Frag.Name()
    else:
      print "Set fragment %s to use hydrophobic contact pairs." % Frag.Name()
  Save()

elif Cmd == "fragcaps":
  #set capping options
  FragList = GetFragList(sys.argv[2])
  s = " ".join(sys.argv[3:]).lower()
  CapN, CapC = True, True
  if "none" in s: CapN, CapC = False, False
  if "both" in s: CapN, CapC = True, True
  if "ncap" in s: CapN = True
  if "ccap" in s: CapC = True
  for Frag in FragList:
    Frag.CapN, Frag.CapC = CapN, CapC
    print "Set fragment %s (N-cap, C-cap) to (%s,%s)." % (Frag.Name(), CapN, CapC)
  Save()
  
elif Cmd == "fragnote":
  #change the note for a fragment
  FragList = GetFragList(sys.argv[2])
  if len(sys.argv) > 3:
    Note = " ".join(sys.argv[3:])
    if Note.lower() == "none": Note = ""
  else:
    Note = ""
  for Frag in FragList:
    Frag.Note = Note
    print "Set note for fragment %s." % Frag.Name()
  Save()  

elif Cmd == "compressdata":
  #save to a gzipped file
  FragList = GetFragList(sys.argv[2])
  FileName = sys.argv[3]
  Args = [x.lower() for x in sys.argv[4:]]
  AddPrep = "prep" in Args
  AddAnal = "anal" in Args
  AddAll = not AddPrep and not AddAnal
  AddDB = "database" in Args or AddPrep
  if AddDB:
    #add the restart.dat file
    Dirs = ["restart.dat"]
  else:
    Dirs = []
  if AddAll:
    #add the base paths
    for Frag in FragList:
      Dirs.append(Frag.BasePath)
  elif AddPrep:
    for Frag in FragList:
      #add conf folder of self and parents
      Dirs.extend([os.path.join(f.BasePath, "conf/")
                   for (f,x) in [(Frag,-1)] + Frag.Parents])
  elif AddAnal:
    for Frag in FragList:
      #add just analysis and conf folders
      Dirs.append(os.path.join(Frag.BasePath, "anal/"))
      Dirs.append(os.path.join(Frag.BasePath, "conf/"))
  Dirs = " ".join(Dirs)
  os.system("tar -cvzf %s %s" % (FileName, Dirs))

elif Cmd == "restartrun":
  #get ready for a restart
  if os.path.isfile("server.txt"):
    os.remove("server.txt")
  print "Ready for new job submission."

elif Cmd == "summarize":
  FragList = GetFragList(sys.argv[2])
  KeyName = sys.argv[3]
  if len(FragList) == 0: sys.exit()
  PathStr = " ".join([os.path.basename(os.path.dirname(f.BasePath))
                      for f in FragList])
  SumCmd1 = "summarize.py . %s -bestn 2 -f %s" % (KeyName, PathStr)
  SumCmd2 = "summarize.py . %s -index -f %s" % (KeyName, PathStr)
  os.system(SumCmd1)
  os.system(SumCmd2)

elif Cmd == "deldata":
  #delete fragment data
  FragList = GetFragList(sys.argv[2])
  Work = not "dataonly" in sys.argv
  Data = not "workonly" in sys.argv
  if ShowWarning(FragList):
    for Frag in FragList:
      #get rid of the data
      Frag.DelData(Work = Work, Data = Data, Conf = False)
      print "Removed data for fragment %s" % Frag.Name()

elif Cmd == "savescores":
  #save scores to a text file
  FragList = GetFragList(sys.argv[2])
  fn = sys.argv[3]
  print "Loading the contact database..."
  cd = zamcontact.ContactDatabaseClass(zd, Update = False)
  cd.UpdateFrags(zd, FragList)
  cd.UpdateScores()
  cd.WriteScores(fn + ".txt")
  cd.WriteHtml(fn + ".html")
  n = len(cd.Scores)
  Fmt = "%-3d %-3d %-8.3f %-5d %-8.3f"
  s = "\n".join([Fmt % (i, j, cd.LogitSum[i,j], cd.Counts[i,j], cd.Scores[i,j])
                 for i in range(n) for j in range(i+3,n) if cd.Counts[i,j] > 0])
  s = "i   j   LogitSum Count Score   \n" + s
  file(fn + ".dat", "w").write(s)
      
elif Cmd == "mss":
  #calculate information content
  FragList = GetFragList(sys.argv[2])
  fn = sys.argv[3]
  from numpy import *
  import protein
  CM = protein.ProteinClass(Pdb = fn).ResContactMap()
  print "Loading the contact database..."
  cd = zamcontact.ContactDatabaseClass(zd, Update = False)
  cd.UpdateFrags(zd, FragList)
  cd.UpdateScores()
  n = len(zd.Seq)
  W1, W2 = 0., 0.
  lnPnat0, lnPnon0 = log(0.5), log(0.5)
  for i in range(n):
    for j in range(i+3,n):
      if cd.Counts[i,j] == 0: continue
      if CM[i,j] == 0:
        W1 += zamcontact.LogitToLogP(-cd.Scores[i,j]) - lnPnon0
      else:
        W1 += zamcontact.LogitToLogP(cd.Scores[i,j]) - lnPnat0
        W2 += zamcontact.LogitToLogP(cd.Scores[i,j]) - lnPnat0
  print "Information metrics: %.4f  %.4f" % (W1, W2)
    
elif Cmd == "mss2":
  #print fragments
  FragList = GetFragList(sys.argv[2])
  fn = sys.argv[3]
  from numpy import *
  import protein, rmsd
  pnat = protein.ProteinClass(Pdb = fn)
  print "Loading the contact database..."
  cd = zamcontact.ContactDatabaseClass(zd, Update = False)
  cd.UpdateFrags(zd, zd.Frags)
  cd.UpdateScores()
  print "Loading the structure database..."
  sd = zamdata.SnapDatabaseClass(zd, cd.Scores, RemoveRedundant = False)
  N = len(zd.Seq)
  CO = fromfunction(lambda i, j: abs(i-j), (N,N), dtype = int)
  P = zamcontact.LogitToP(cd.Scores)
  s = ""
  for Frag in FragList:
    a, b = Frag.StartRes, Frag.StopRes + 1
    ScoreCO = sum(P[a:b,a:b] * CO[a:b,a:b]) / sum(P[a:b,a:b])
    Snaps = sd.GetSnaps(Frag)
    for Snap in Snaps:
      p = protein.ProteinClass(Pdb = Snap.Pdb)
      ThisCO = sum([abs(x-y) for (x,y) in Snap.Contacts]) / float(len(Snap.Contacts))
      r = rmsd.RMSDProteinClass(pnat, p, AlignSeq = True)[0]
      s += "%s %.2f %.2f %.2f %.2f\n" % (Snap.Name, r, Snap.Score, ThisCO, ScoreCO)
  file("co.dat", "w").write(s)

else:
  print "Command not recognized"
  
