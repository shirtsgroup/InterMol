#!/usr/bin/env python

#last modified: 05-01-07

#CONVENTION: Any command that takes Pdb as an argument can use
#either the pdb filename or the actual pdb data 

#CONVENTION: Residue and atom numbers start as 1 as in the Pdb file



Usage = """Modifies pdb files.

Usage     : pdbtools.py cap INPUTPDB [OUTPUTPDB]
   OR     : pdbtools.py decap INPUTPDB [OUTPUTPDB]
   OR     : pdbtools.py dehydrogen INPUTPDB [OUTPUTPDB]
   OR     : pdbtools.py extendc INPUTPDB SEQ [OUTPUTPDB]
   OR     : pdbtools.py extendn INPUTPDB SEQ [OUTPUTPDB]
   OR     : pdbtools.py extend INPUTPDB SEQBEFORE SEQAFTER [OUTPUTPDB]
   OR     : pdbtools.py generate SEQ [OUTPUTPDB]
   OR     : pdbtools.py renumber INPUTPDB OUTPUTPDB
   OR     : pdbtools.py rotate INPUTPDB RESNUM PHI PSI [OUTPUTPDB]
   OR     : pdbtools.py showseq INPUTPDB
   OR     : pdbtools.py showoverlaps INPUTPDB [OVERLAPDIST]
   OR     : pdbtools.py splice INPUTPDB1 INPUTPDB2 OUTPUTPDB
   OR     : pdbtools.py spliceopt INPUTPDB1 INPUTPDB2 OUTPUTPDB
   OR     : pdbtools.py trim INPUTPDB STARTRES STOPRES [OUTPUTPDB]
   OR     : pdbtools.py standardize INPUTPDB [OUTPUTPDB]
   OR     : pdbtools.py alignseq INPUTPDB1 INPUTPDB2

INPUTPDB  : input pdb file
OUTPUTPDB : output pdb file
STARTRES  : starting residue number
STOPRES   : stopping residue number
RESNUM    : residue number
SEQ       : one-letter sequence codes
PHI,PSI   : desired angles after rotation, in degrees
"""

#check for instructions
import sys
if __name__ == "__main__" and len(sys.argv) == 1:
  print Usage
  sys.exit()
  

import re, string, urllib, os, tempfile, copy, random, gzip
import sequence, protein


#global variables
Caps = ["ACE", "NHE", "NME"]
TerminiAtoms = ["OXT", "H2", "H3"]
OverlapDist = 0.65    #cutoff for detecting atomic overlap, in angstrom


def IsFileName(PdbString):
  """Returns True if PdbString is a filename, False if it is the data."""
  return not "\n" in PdbString

def IsData(PdbString):
  """Returns True if PdbString is the Pdb data, False if it is the filename."""
  return "\n" in PdbString

def ReturnPdbData(PdbString):
  """PdbString can be a pdb filename or the actual file contents."""
  if IsData(PdbString):
    return PdbString
  else:
    if os.path.isfile(PdbString):
      #check for gzip
      if PdbString.endswith(".gz"):
        return gzip.GzipFile(PdbString,"r").read()
      else:
        return file(PdbString,"rU").read()
    else:
      raise IOError, "Pdb file %s not found." % PdbString

def SavePdbData(PdbString, PdbFile):
  "Saves a pdb string to a file."
  if PdbFile is not None:
    file(PdbFile, "w").write(PdbString)
    

#-------- INFORMATION --------

def Atoms(Pdb):
  "Returns an array of atom names."
  Pdb = ReturnPdbData(Pdb)
  return [line[12:16] for line in Pdb.split("\n") if line.startswith("ATOM")]

def AtomRes(Pdb):
  "Returns an array of residue names for each atom."
  Pdb = ReturnPdbData(Pdb)
  return [line[22:26] for line in Pdb.split("\n") if line.startswith("ATOM")]

def AtomNum(Pdb):
  "Returns an array of atom numbers for each atom."
  Pdb = ReturnPdbData(Pdb)
  return [int(line[6:11]) for line in Pdb.split("\n") if line.startswith("ATOM")]

def AtomPos(Pdb):
  "Returns an array of atom coordinates."
  Pdb = ReturnPdbData(Pdb)
  return [[float(l[30:38]), float(l[38:46]), float(l[46:54])]
          for l in Pdb.split("\n") if l.startswith("ATOM")]

def AtomResNums(Pdb):
  "Returns the numbers of the residues for each atom the pdb file data."
  Pdb = ReturnPdbData(Pdb)
  r = []
  ThisResNum = ""
  for l in Pdb.split("\n"):
    if (not l.startswith("ATOM")):
      continue
    ResNum = l[22:29]
    if ThisResNum != ResNum:
      r.append(ResNum)
  return r 

def Seq(Pdb):
  "Returns an array of residue names for the sequence."
  Pdb = ReturnPdbData(Pdb)
  s = []
  ThisResNum = ""
  for l in Pdb.split("\n"):
    if (not l.startswith("ATOM")):
      continue
    ResNum = l[22:29]
    if ThisResNum != ResNum:
      s.append(l[17:20])
      ThisResNum = ResNum
  return s

def ResLen(Pdb):
  "Returns the number of residues."
  return len(Seq(Pdb))

def AtomLen(Pdb):
  "Returns the number of atoms."
  return len(Atoms(Pdb))

def ResNums(Pdb):
  "Returns the numbers of the residues in the pdb file data."
  Pdb = ReturnPdbData(Pdb)
  r = []
  ThisResNum = ""
  for l in Pdb.split("\n"):
    if (not l.startswith("ATOM")):
      continue
    ResNum = l[22:29]
    if ThisResNum != ResNum:
      r.append(ResNum)
  return r  

def GetCoords(Pdb):
  "Returns the Pdb coordinates."
  Pdb = ReturnPdbData(Pdb)
  Coords = [[float(s[30:38]), float(s[38:46]), float(s[46:54])]
            for s in Pdb.split("\n") if s.startswith("ATOM")]
  return Coords

def GetOverlaps(Pdb, MinDist = None, FirstOnly = False):
  """Returns a list of tuples (a,b) of overlaps between atoms a and b.
a and b start at zero and increment consecutively and are not taken
from the pdb data."""
  if MinDist is None: MinDist = OverlapDist
  #load atom positions
  Coords = GetCoords(Pdb)
  N = len(Coords)
  Overlaps = []
  MinDistSq = MinDist*MinDist
  #check atom distances
  for i in range(0,N):
    for j in range(i+1,N):
      DistSq = sum([(Coords[i][k]-Coords[j][k])**2 for k in range(0,3)])
      if DistSq < MinDistSq:
        Overlaps.append((i,j))
        if FirstOnly: return Overlaps
  return Overlaps

def HasOverlap(Pdb, MinDist = OverlapDist):
  return len(GetOverlaps(Pdb, MinDist = MinDist, FirstOnly = True)) > 0

def HasHydrogens(Pdb):
  "Indicates whether or not there are hydrogens."
  Pdb = ReturnPdbData(Pdb)
  OutPdb = [s for s in Pdb.split("\n") if s[0:4]=="ATOM" and s[13]=="H"]
  return len(OutPdb) > 0

def HasCaps(Pdb):
  "Indicates if there are cap residues in a pdb file."
  Pdb = ReturnPdbData(Pdb)
  OutPdb = [s for s in Pdb.split("\n") if s[0:4]=="ATOM" and s[17:20] in Caps]
  return len(OutPdb) > 0


#--------PDB MODIFIERS--------

def Renumber(Pdb, StartAtom = 1, StartRes = 1, StartChain = " ",
             OutPdbFile = None):
  """Renumbers pdb residues and atoms starting at StartRes and StartAtoms,
  and starting the chain at no specification."""
  ChainLetters = [" "] + list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
  Pdb = ReturnPdbData(Pdb)
  OutPdb = []
  r = StartRes - 1
  a = StartAtom - 1
  c = ChainLetters.index(StartChain)
  ThisResNum = ""
  ThisChain = ""
  #first get all the chains
  Chains = []
  for l in Pdb.split("\n"):
    if l.startswith("ATOM"):
      if not l[21] in Chains:
        Chains.append(l[21])
  Chains.sort()
  ChainDict = {}
  for i in range(0, len(Chains)):
    ind = (c + i) % len(ChainLetters)
    ChainDict[Chains[i]] = ChainLetters[ind]
  #now renumber
  for l in Pdb.split("\n"):
    if l.startswith("ATOM"):
      a += 1
      ResNum = l[22:29]
      if ThisResNum != ResNum:
        r += 1
        ThisResNum = ResNum
      ChainNum = ChainDict[l[21]]
      s = l[:6] + str(a).rjust(5) + l[11:21] + ChainNum \
          + str(r).rjust(4) + "   " + l[29:]
      OutPdb.append(s)
    else:
      OutPdb.append(l)
  OutPdb = "\n".join(OutPdb)
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb 

def Dehydrogen(Pdb, OutPdbFile = None):
  "Removes hydrogens from a pdb file."
  Pdb = ReturnPdbData(Pdb)
  OutPdb = [s for s in Pdb.split("\n") if not (s[0:4]=="ATOM" and s[13]=="H")]
  OutPdb = "\n".join(OutPdb)
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Decap(Pdb, OutPdbFile = None):
  "Removes caps and terminal atoms from a pdb file."
  Pdb = ReturnPdbData(Pdb)
  OutPdb = []
  for l in Pdb.split("\n"):
    if l.startswith("ATOM"):
      if l[12:16].strip() == "H1":
        OutPdb.append(l[:12] + " H  " + l[16:])
      elif not l[12:16].strip() in TerminiAtoms and not l[17:20] in Caps:
        OutPdb.append(l)
    else:
      OutPdb.append(l)
  OutPdb = "\n".join(OutPdb)
  OutPdb = Renumber(OutPdb)
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Cap(Pdb, OutPdbFile = None):
  "Adds caps ACE and NME to a pdb file."
  #this uses proteinclass
  Pdb = ReturnPdbData(Pdb)
  p = protein.ProteinClass(Pdb = Pdb)
  CapNRes = protein.ProteinClass(Seq = ">")
  CapCRes = protein.ProteinClass(Seq = "<")
  p = CapNRes + p + CapCRes
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def CapN(Pdb, OutPdbFile = None):
  "Adds n-terminal ACE cap to a pdb file."
  #this uses proteinclass
  Pdb = ReturnPdbData(Pdb)
  p = protein.ProteinClass(Pdb = Pdb)
  CapNRes = protein.ProteinClass(Seq = ">")
  p = CapNRes + p
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def CapC(Pdb, OutPdbFile = None):
  "Adds c-terminal NME cap to a pdb file."
  #this uses proteinclass
  Pdb = ReturnPdbData(Pdb)
  p = protein.ProteinClass(Pdb = Pdb)
  CapCRes = protein.ProteinClass(Seq = "<")
  p = p + CapCRes
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Splice(Pdb1, Pdb2, OutPdbFile = None):
  "Joins two pdb files together."
  #this uses proteinclass
  Pdb1 = ReturnPdbData(Pdb1)
  Pdb2 = ReturnPdbData(Pdb2)
  p1 = protein.ProteinClass(Pdb = Pdb1)
  p2 = protein.ProteinClass(Pdb = Pdb2)
  p3 = p1 + p2
  OutPdb = p3.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Generate(Seq, OutPdbFile = None):
  "Generates a pdb file of extended sequence Seq."
  #uses proteinclass
  p = protein.ProteinClass(Seq = Seq)
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Extend(Pdb, SeqBefore, SeqAfter, OutPdbFile = None, Opt = True, N = 5):
  "Extends a pdbfile in either or both directions with given sequences."
  #uses proteinclass
  OutPdb = ReturnPdbData(Pdb)
  p = protein.ProteinClass(Pdb = Pdb)
  n = len(p)
  if len(SeqBefore) > 0:
    p = protein.ProteinClass(Seq = SeqBefore) + p
    if Opt and n > 0: p.OptimizePep(len(SeqBefore), N)
  if len(SeqAfter) > 0:
    p = p + protein.ProteinClass(Seq = SeqAfter)
    if Opt and n > 0: p.OptimizePep(n, N)
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Extract(Pdb, StartRes, StopRes, OutPdbFile = None):
  "Extracts a portion of a pdb file between specified residue nums."
  Pdb = ReturnPdbData(Pdb)
  OutPdb = []
  ThisResNum = ""
  ResIndex = 0
  for l in Pdb.split("\n"):
    if not l.startswith("ATOM"): continue
    ResNum = l[22:29]
    if ThisResNum != ResNum:
      ThisResNum = ResNum
      ResIndex += 1
    if ResIndex >= StartRes and ResIndex <= StopRes:
      OutPdb.append(l)
  OutPdb = "\n".join(OutPdb)
  OutPdb = Renumber(OutPdb)
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb
    
def Trim(Pdb, Chains = [], ResNums = [], OutPdbFile = None):
  """Selects out specific chains or residue numbers from a pdb file.
Here, residue nums are those in the file, and not necessarily in
sequential order starting at 1."""
  Pdb = ReturnPdbData(Pdb)
  pdblines = Pdb.split("\n")
  OutPdb = []
  for line in pdblines:
    if (not line.startswith("ATOM")):
      if line.startswith("TER"):
          line = line[:17] + lastres + line[20:23] + lastresnum + line[26:]
      OutPdb.append(line)
      continue
    chain = line[21]
    resnum = int(line[23:26])
    if (len(Chains) and chain not in Chains):
      continue
    if (len(ResNums) and resnum not in ResNums):
      continue
    lastres = line[17:20]
    lastresnum = line[23:26]
    OutPdb.append(line)
  OutPdb = "\n".join(OutPdb)
  OutPdb = Renumber(OutPdb)
  SavePdbData(OutPdb, OutPdbFile)    
  return OutPdb

def Rotate(Pdb, ResNum, Phi = None, Psi = None, OutPdbFile = None):
  "Rotates a residue to specified phi, psi."
  #this uses proteinclass
  Pdb = ReturnPdbData(Pdb)
  p = protein.ProteinClass(Pdb = Pdb)
  if not ResNum in range(0, len(p)):
    raise IndexError, "Residue number %d not found" % ResNum
  try:
    p.RotateToPhiPsi(ResNum, Phi, Psi)
  except StandardError:
    print "Could not perform rotation."
    return Pdb
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def SpliceOpt(Pdb1, Pdb2, OutPdbFile = None, N = 5):
  "Splices two pdbs together and rotates to optimize for non-overlap."
  #this uses proteinclass
  Pdb1 = ReturnPdbData(Pdb1)
  Pdb2 = ReturnPdbData(Pdb2)
  p1 = protein.ProteinClass(Pdb = Pdb1)
  p2 = protein.ProteinClass(Pdb = Pdb2)
  p3 = p1 + p2
  if len(p1) > 0 and len(p2) > 0:
    p3.OptimizePep(len(p1), N)
  OutPdb = p3.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def RandDihedrals(Pdb, OutPdbFile = None, DeltaAng = 5.):
  "Adds a random angle to the torsions along the backbone."
  #this uses proteinclass
  Pdb = ReturnPdbData(Pdb)
  p = protein.ProteinClass(Pdb = Pdb)
  for i in range(0,len(p)):
    Phi, Psi = p.PhiPsi(i)
    if not Phi is None and not Psi is None:
      Phi += DeltaAng * (2.*random.random() - 1.)
      Psi += DeltaAng * (2.*random.random() - 1.)
      p.RotateToPhiPsi(i, Phi, Psi)
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Standardize(Pdb, OutPdbFile = None):
  "Formats to standard residue and atom names, etc."
  Pdb = ReturnPdbData(Pdb)
  OutPdb = []
  for l in Pdb.split("\n"):
    if l.startswith("ATOM"):
      #fix histidines
      if l[17:20] in ["HID","HIE","HIP"]:
        l = l[:17] + "HIS" + l[20:]
      elif l[17:20] in ["CYX", "CYM"]:
        l = l[:17] + "CYS" + l[20:]
      #fix amide hydrogens
      if l[12:16] == " HN ":
        l = l[:12] + " H  " + l[16:]
    OutPdb.append(l)
  OutPdb = "\n".join(OutPdb)
  OutPdb = Renumber(OutPdb)
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb  
     

#--------PDB DOWNLOAD--------

def Download(Id, OutPdbFile = None):
  "Downloads a pdb file."
  pdbmirror = "ftp://ftp.rcsb.org/pub/pdb"
  if (not re.match("^[1-9][1-9a-z]{3}", Id)):
    print "warning: " + Id + " not a valid PDB code"
    return ""
  Id = string.lower(Id)
  url = pdbmirror + "/data/structures/all/pdb/pdb" + Id + ".ent.Z"
  urlhandle = urllib.urlopen(url)
  pdbz = urlhandle.read()
  urlhandle.close()
  OutPdb = __Uncompress(pdbz)

  SavePdbData(OutPdb, OutPdbFile)  
  return OutPdb

def __Uncompress(Data):
  "Accessory for uncompressing downloaded pdb."
  #uin, uout = os.popen2("uncompress")
  #
  #uin.write(data)
  #uin.close()
  
  compfd, comppath = tempfile.mkstemp(".Z")
  compfile = os.fdopen(compfd, "w")
  compfile.write(Data)
  compfile.close()
  
  uout = os.popen("uncompress -c " + comppath)
  
  uncompdata = uout.read()
  uout.close()
  
  return uncompdata



def GetPdbFilenames(argv, NArg):
  if len(argv) > 3 + NArg:
    return sys.argv[2], sys.argv[-1]
  else:
    return sys.argv[2], sys.argv[2]

#command-line running
if __name__ == "__main__":
  Cmd = sys.argv[1].lower()

  if Cmd == "cap":
    InPdb, OutPdb = GetPdbFilenames(sys.argv, 0)
    pdb = Cap(InPdb, OutPdbFile = OutPdb)
    print "Capped %s as %s." % (InPdb, OutPdb)

  elif Cmd == "decap":
    InPdb, OutPdb = GetPdbFilenames(sys.argv, 0)
    pdb = Decap(InPdb, OutPdbFile = OutPdb)
    print "Decapped %s as %s." % (InPdb, OutPdb)

  elif Cmd == "dehydrogen":
    InPdb, OutPdb = GetPdbFilenames(sys.argv, 0)
    pdb = Dehydrogen(InPdb, OutPdbFile = OutPdb)
    print "Dehydrogenated %s as %s." % (InPdb, OutPdb)

  elif Cmd == "extendc":
    InPdb, OutPdb = GetPdbFilenames(sys.argv, 1)
    Seq = sys.argv[3]
    pdb = Extend(InPdb, [], Seq, OutPdbFile = OutPdb, Opt = True)
    print "Extended and optimized for overlaps %s as %s." % (InPdb, OutPdb)

  elif Cmd == "extendn":
    InPdb, OutPdb = GetPdbFilenames(sys.argv, 1)
    Seq = sys.argv[3]
    pdb = Extend(InPdb, Seq, [], OutPdbFile = OutPdb, Opt = True)
    print "Extended and optimized for overlaps %s as %s." % (InPdb, OutPdb)

  elif Cmd == "extend":
    InPdb, OutPdb = GetPdbFilenames(sys.argv, 2)
    SeqBefore, SeqAfter = sys.argv[3:5]
    pdb = Extend(InPdb, SeqBefore = SeqBefore, SeqAfter = SeqAfter,
                 OutPdbFile = OutPdb, Opt = True)
    print "Extended and optimized for overlaps %s as %s." % (InPdb, OutPdb)

  elif Cmd == "generate":
    Seq, OutPdb = sys.argv[2:4]
    pdb = Generate(Seq, OutPdbFile = OutPdb)
    print "Generated as %s." % OutPdb

  elif Cmd == "renumber":
    InPdb, OutPdb = GetPdbFilenames(sys.argv, 0)
    pdb = Renumber(InPdb, OutPdbFile = OutPdb)
    print "Renumbered %s as %s." % (InPdb, OutPdb)

  elif Cmd == "rotate":
    InPdb, OutPdb = GetPdbFilenames(sys.argv, 3)
    ResNum = int(sys.argv[3]) - 1
    Phi, Psi = [float(x) for x in sys.argv[4:6]]
    pdb = Rotate(InPdb, ResNum, Phi, Psi, OutPdbFile = OutPdb)
    print "Rotated %s as %s." % (InPdb, OutPdb)    

  elif Cmd == "showseq":
    InPdb = sys.argv[2]
    Seq = sequence.SeqToList(Seq(InPdb))
    Nums = range(1, len(Seq)+1)
    n = 15
    i = 0
    print "\n" + sequence.SeqToAA1(Seq) + "\n"
    while i < len(Seq):
      print " ".join(["%-3d" % x for x in Nums[i:i+n]])
      print " ".join(["%-3s" % x for x in Seq[i:i+n]])
      print ""
      i += n

  elif Cmd == "showoverlaps":
    InPdb = sys.argv[2]
    if len(sys.argv) > 3:
      OverlapDist = float(sys.argv[3])
    aname = [x.strip() for x in Atoms(InPdb)]
    anum = AtomNum(InPdb)
    Overlaps = GetOverlaps(InPdb)
    if len(Overlaps) == 0:
      print "No overlaps detected."
    else:
      print "Overlaps (less than %.2f A):" % (OverlapDist)
      for (a,b) in Overlaps:
        print "  %s%d-%s%d" % (aname[a], anum[a], aname[b], anum[b])

  elif Cmd == "splice":
    InPdb1, InPdb2, OutPdb = sys.argv[2:5]
    pdb = Splice(InPdb1, InPdb2, OutPdbFile = OutPdb)
    print "Spliced %s and %s as %s." % (InPdb1, InPdb2, OutPdb)

  elif Cmd == "spliceopt":
    InPdb1, InPdb2, OutPdb = sys.argv[2:5]
    pdb = SpliceOpt(InPdb1, InPdb2, OutPdbFile = OutPdb)
    print "Spliced and optimized for overlaps %s and %s as %s." \
          % (InPdb1, InPdb2, OutPdb)

  elif Cmd == "trim":
    InPdb, OutPdb = GetPdbFilenames(sys.argv, 2)
    StartRes = int(sys.argv[3])
    StopRes = int(sys.argv[4])
    pdb = Extract(InPdb, StartRes, StopRes, OutPdbFile = OutPdb)
    print "Extracted residues from %s as %s." % (InPdb, OutPdb)

  elif Cmd == "standardize":
    InPdb, OutPdb = GetPdbFilenames(sys.argv, 2)
    pdb = Standardize(InPdb, OutPdbFile = OutPdb)
    print "Standardized pdb file from %s to %s" % (InPdb, OutPdb)

  elif Cmd == "alignseq":
    InPdb1, InPdb2 = sys.argv[2:]
    Seq1, Seq2 = Seq(InPdb1), Seq(InPdb2)
    Map = sequence.SeqMapClass(Seq1, Seq2)
    if len(Map) <= 1:
      print "Pdb files do not align."
    else:
      print "%s: residues %d to %d (%d of %d)" % (InPdb1, Map.a + 1, Map.b,
                                                  len(Map), len(Seq1))
      print "%s: residues %d to %d (%d of %d)" % (InPdb2, Map.c + 1, Map.d,
                                                  len(Map), len(Seq2))
      Seq = Seq1[Map.a:Map.b]
      Nums1, Nums2 = range(Map.a+1, Map.b+1), range(Map.c+1, Map.d+1)
      n = 15
      i = 0
      print "\nOverlap sequence:"
      print sequence.SeqToAA1(Seq) + "\n"
      print sequence.SeqToAA1(Seq2[Map.c:Map.d]) + "\n"
      while i < len(Seq):
        print " ".join(["%-3d" % x for x in Nums1[i:i+n]])
        print " ".join(["%-3d" % x for x in Nums2[i:i+n]])
        print " ".join(["%-3s" % x for x in Seq[i:i+n]])
        print ""
        i += n
  else:
    print "Command not recognized."
    