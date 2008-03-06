#!/usr/bin/env python

#LAST MODIFIED: 05-04-07

#DESCRIPTION: Provides reading and writing routines for popular
#coordinate file formats.  Currently: pdb and amber trajectory.

from numpy import *
import copy, os, gzip, pdbtools

#Masks for backbone atoms
NoMask = []
PdbBackboneMask = ["CA", "C", "N"]
CrdBackboneMask = ["CA", "C", "N"]
BackboneMask = ["CA", "C", "N"]
AlphaCarbonMask = ["CA"]
Caps = ["NHE","NME","ACE"]


def GetPdbSeq(PdbFile):
  "Gets the sequence from a PdbFile."
  LineList = [l for l in open(PdbFile, "rU").readlines() if l[0:4] == 'ATOM']
  Seq = []
  ThisResNum = 0
  for l in LineList:
    ResNum = int(l[22:29])
    if ThisResNum != ResNum:
      Seq.append(l[17:20].upper())
      ThisResNum = ResNum
  return Seq    

def GetPdbCoords(PdbFile, Mask = NoMask):
  "Gets coordinates from a Pdb file."
  if os.path.isfile(PdbFile):
    f = open(PdbFile, "rU")
    Coords = []
    while True:
      s = f.readline()
      if len(s) == 0: break
      if s[0:4] == "ATOM" and (s[13:15].strip() in Mask or Mask == NoMask):
        Coords.append([float(s[30:38]),float(s[38:46]), float(s[46:54])])
    f.close()
    return array(Coords,float)
  else:
    raise IOError, "Pdb file does not exist."

def SavePdbCoords(Coords, PdbFile, PdbTmplFile, Mask = NoMask):
  "Saves a coordinate array to a pdb file using a template pdb file."
  if os.path.isfile(PdbTmplFile):
    ft = open(PdbTmplFile, "rU")
    f = open(PdbFile, "w")
    AtomNum = 0
    while True:
      s = ft.readline()
      if len(s) == 0: break
      if s.startswith("TER") or s.startswith("MODEL") or s.startswith("ENDMDL"):
        f.write(s)
      elif s.startswith("ATOM") and (s[13:15].strip() in Mask or Mask == NoMask):
        s = s[0:30] + "%8.3f%8.3f%8.3f" % (Coords[AtomNum,0],
            Coords[AtomNum,1], Coords[AtomNum,2]) + s[54:]
        AtomNum += 1
        f.write(s)
    ft.close()
    f.close()
    if not AtomNum == len(Coords):
      print "Warning: pdb template file %s has %d coordinates, but I have %d." \
            % (PdbTmplFile, AtomNum, len(Coords))
  else:
    raise IOError, "Warning: could not find template pdb file."

def ParseCrdString(s, AtomNames = [], Mask = NoMask):
  """Takes a text string of Crd coordinates and parses into an array."""
  #parse into a n by 3 array
  try:
    vals = s.replace("\n","")
    vals = [float(vals[i:i+8]) for i in range(0, len(vals), 8)]
    Coords = array(vals, float)
  except ValueError:
    return
  if mod(len(Coords), 3) == 0:
    Coords = reshape(Coords, (-1,3))
    #if using mask remove the extraneous coordinates
    if not Mask == NoMask and len(AtomNames) == len(Coords):
      Coords = compress([a.strip() in Mask for a in AtomNames], Coords, 0)
    return Coords
  else:
    raise IOError, "Improper number of coordinates found in Crd file."


def ParseRstString(s, AtomNames = [], Mask = NoMask):
  """Takes a text string of Rst coordinates and parses into an array;
Note that this will return velocities as well."""
  #parse into a n by 3 array
  try:
    vals = s.replace("\n","")
    vals = [float(vals[i:i+12]) for i in range(0, len(vals), 12)]
    Coords = array(vals, float)
  except ValueError:
    return
  if mod(len(Coords), 3) == 0:
    Coords = reshape(Coords, (-1,3))
    #if using mask remove the extraneous coordinates
    if not Mask == NoMask and len(AtomNames) == len(Coords):
      Coords = compress([a.strip() in Mask for a in AtomNames], Coords, 0)
    return Coords
  else:
    raise IOError, "Improper number of coordinates found in Rst file."


def GetPrmtopBlock(PrmtopFile, Flag, BlockLen = None):
  "Gets Prmtop data for a specified Flag."
  if os.path.isfile(PrmtopFile):
    f = open(PrmtopFile, "rU")
    Dat = []
    line = f.readline()
    while len(line) > 0 and not line.startswith(Flag):
      line = f.readline()
    f.readline()
    line = f.readline()
    while len(line) > 0 and not line[0:1] == "%":
      if BlockLen is None:
        Dat.extend(line.split())
      else:
        Dat.extend([line[i:i+BlockLen] for i in xrange(0,len(line)-1,BlockLen)])
      line = f.readline()
    f.close()
    return Dat
  else:
    raise IOError, "Could not find Prmtop file."
    return None

def GetPrmtopAtomNames(PrmtopFile):
  "Gets the names of atoms from a Prmtop file."
  Names = GetPrmtopBlock(PrmtopFile, "%FLAG ATOM_NAME", BlockLen = 4)
  return Names

def AmbToPdbAtomNames(Names):
  "Converts Amber-style atom names to Pdb-style."
  return [an[-1] + an[:-1] for an in Names]

def GetPrmtopSeq(PrmtopFile):
  "Gets the names of residues from a Prmtop file."
  Seq = GetPrmtopBlock(PrmtopFile, "%FLAG RESIDUE_LABEL", BlockLen = 4)
  Seq = [s.strip() for s in Seq]
  return Seq

def GetPrmtopAtomRes(PrmtopFile):
  "Gets the atom's residue numbers from a Prmtop file."
  ResPtr = GetPrmtopBlock(PrmtopFile, "%FLAG RESIDUE_POINTER", BlockLen = 8)
  if ResPtr is None: return None
  AtomNames = GetPrmtopAtomNames(PrmtopFile)
  ResPtr = [int(x) for x in ResPtr] + [len(AtomNames) + 1]
  AtomRes = []
  for i in range(0, len(ResPtr)-1):
    AtomRes.extend([i] * (ResPtr[i+1] - ResPtr[i]))
  return AtomRes


def GetCrdCoords(CrdFile, PrmtopFile = "", Mask = NoMask):
  "Gets coordinates from a Crd and Prmtop set of files."
  if len(Mask) > 0:
    AtomNames = GetPrmtopAtomNames(PrmtopFile)
    if AtomNames is None: return
  else:
    AtomNames = []
  if os.path.isfile(CrdFile):
    f = open(CrdFile, "rU")
    f.readline()  #bypass the header lines
    f.readline()
    s = f.read()
    f.close()
    return ParseCrdString(s, AtomNames, Mask)
  else:
    raise IOError, "Crd file not found."


def SaveCrdCoords(Coords, CrdFile):
  "Saves coordinates to a Crd file."
  NPerLine = 10
  Fmt = "%8.3f"
  f = open(CrdFile, "w")
  f.write("ACE".ljust(80) + "\n")
  p = [Fmt % x for x in Coords.flat]
  r = mod(len(p),NPerLine)
  if r > 0: p.extend([""]*(NPerLine-r))
  for i in xrange(0, len(p), NPerLine):
    f.write("".join(p[i:i+NPerLine]) + "\n")
  f.close()


def GetRstCoords(RstFile, PrmtopFile = "", Mask = NoMask):
  "Gets coordinates from an amber restart file."
  if len(Mask) > 0:
    AtomNames = GetPrmtopAtomNames(PrmtopFile)
    if AtomNames is None: return
  else:
    AtomNames = []
  if os.path.isfile(RstFile):
    f = open(RstFile, "rU")
    f.readline()  #bypass the header lines
    f.readline()
    s = f.read()
    f.close()
    return ParseRstString(s, AtomNames, Mask)
  else:
    raise IOError, "Rst file not found."


def SaveRstCoords(Coords, RstFile):
  "Saves coordinates to an amber restart file."
  NPerLine = 6
  Fmt = "%12.7f"
  f = open(RstFile, "w")
  f.write("ACE".ljust(80) + "\n")
  f.write("%5d  0.0000000E+00\n" % len(Coords))
  p = [Fmt % x for x in Coords.flat]
  r = mod(len(p),NPerLine)
  if r > 0: p.extend([""]*(NPerLine-r))
  for i in xrange(0, len(p), NPerLine):
    f.write("".join(p[i:i+NPerLine]) + "\n")
  f.close()  
  

def AmbToPdb(RstFile, PrmtopFile, PdbFile, AAtm = False, BRes = False):
  """Covert parmtop and coordinate files into a PDB string
* RstFile: file path to amber restart file file
* PrmtopFile: file path to parmtop file
* PdbFile: pdb file to write
* AAtm: boolean, use Amber atom names in the PDB file
* BRes: boolean, use PDB-standard residue names"""
  ambcmd = "ambpdb -p " + PrmtopFile
  if (AAtm):
    ambcmd += " -aatm"
  if (BRes):
    ambcmd += " -bres"
  ambpdb = os.popen(ambcmd + " < " + RstFile + " 2> /dev/null", "r")
  pdb = ambpdb.read()
  ambpdb.close()
  file(PdbFile,"w").write(pdb)


def SavePdbCoordsAmb(Coords, PrmtopFile, PdbFile, Standardize = True):
  """Saves a coordinate array to a pdb file using an amber prmtop file."""
  AtomNames = GetPrmtopAtomNames(PrmtopFile)
  AtomNames = AmbToPdbAtomNames(AtomNames)
  AtomRes = GetPrmtopAtomRes(PrmtopFile)
  Seq = GetPrmtopSeq(PrmtopFile)
  if Standardize:
    for i in range(0, len(Seq)):
      if Seq[i] in ["HIE","HID","HIP"]:
        Seq[i] = "HIS"
  Pdb = "MODEL\n"
  for i in range(0, len(AtomNames)):
    (x, y, z) = Coords[i,:]
    an, ar = AtomNames[i], Seq[AtomRes[i]]
    Pdb += "ATOM  %5d %4s %3s  %4d    %8.3f%8.3f%8.3f \n" \
      % (i+1, an, ar, AtomRes[i] + 1, x, y, z)
  Pdb += "TER\nENDMDL\n"
  file(PdbFile, "w").write(Pdb)


def SavePdbCoordsAmb2(Coords, PrmtopFile, PdbFile, AAtm = False, BRes = False):
  """Saves a coordinate array to a pdb file using an amber prmtop file."""
  SaveRstCoords(Coords, "rst.tmp")
  AmbToPdb("rst.tmp", PrmtopFile, PdbFile, AAtm, BRes)
  os.remove("rst.tmp")


def GetTrjLen(TrjFile, PrmtopFile):
  "Gets the number of frames in a trajectory."
  #get the number of atoms using the prmtop file
  Names = GetPrmtopAtomNames(PrmtopFile)
  if Names is None: return 0
  NAtom = len(Names)
  #calculate the number of lines to read per coordinate set
  BlockLen = (3 * NAtom) / 10
  if (3 * NAtom) % 10 > 0: BlockLen += 1
  #use gzip?
  if TrjFile.split(".")[-1].strip().lower() == "gz":
    FileMthd = gzip.GzipFile
  else:
    FileMthd = file
  #open the file and count the lines
  NLine = -1
  f = FileMthd(TrjFile, "r")
  while not f.readline() == "":
    NLine += 1
  f.close()
  return NLine/BlockLen


class TrjClass:
  "Provides a class for reading successive sets of coordinates from Amber trajectory files."
  
  def __init__(self, TrjFile, PrmtopFile, Mask = NoMask,
               NSkip = 0, NRead = None, NStride = 1,
               LinkPos = None):
    """Initializes the class and opens the trajectory file for reading.
* TrjFile: string name of trj file
* PrmtopFile: string name of prmtop file
* Mask: list of strings; filter for atom names (default is no mask/empty list)
* NSkip: number of configurations to skip
* NRead: maximum number of configurations to read (default is all)
* NStride: stride between configuration frames (default is 1)
* LinkPos: an outside array that is updated automatically as coords are read"""
    IsFile1, IsFile2 = os.path.isfile(PrmtopFile), os.path.isfile(TrjFile)
    if IsFile1 and IsFile2:
      #set the filenames
      self.TrjFile = TrjFile
      self.PrmtopFile = PrmtopFile
      #set the frames to skip, read, and stride
      self.NSkip = NSkip
      self.NRead = NRead
      self.NStride = max(NStride, 1)
      if self.NRead < 0: self.NRead = None
      #set the byte fields
      self.BytesHead = 0
      self.BytesCoords = 0
      self.BytesTot = 0
      self.NCoords = 0
      self.SliceNCoords = 0
      #set the counters
      self.Count = 0
      self.Index = -1
      self.SliceIndex = -1
      #get the atom names, seq, and residue nums
      self.AtomNames = GetPrmtopAtomNames(self.PrmtopFile)
      self.AtomNames = AmbToPdbAtomNames(self.AtomNames)
      self.NAtom = len(self.AtomNames)
      if self.NAtom == 0:
        raise IOError, "There was an error reading the Prmtop file."
        return
      self.AtomRes = GetPrmtopAtomRes(self.PrmtopFile)
      self.Seq = GetPrmtopSeq(self.PrmtopFile)
      #set the mask option
      self.Mask = Mask
      #set the linked pos
      self.LinkPos = LinkPos
      #set the file method
      if self.TrjFile.split(".")[-1].strip().lower() == "gz":
        self.__FileMthd = gzip.GzipFile
      else:
        self.__FileMthd = open
      self.__Trj = None
      #count and reset
      self.__CountBytes()
      self.Reset()

    elif IsFile1:
      raise IOError, "Could not find %s." % TrjFile
    else:
      raise IOError, "Could not find %s." % PrmtopFile


  def __CountBytes(self):
    """Counts the number of bytes per configuration by reading the first set."""
    #calculate the number of lines to read per coordinate set
    r = mod(3 * self.NAtom, 10)
    self.__NLine = (3 * self.NAtom) / 10
    if r > 0: self.__NLine += 1
    #close the trj file if it's open
    self.Close()
    #open the trj file
    try:
      self.__Trj = self.__FileMthd(self.TrjFile, "r")
    except IOError:
      raise IOError, "There was an error opening the trajectory file."
      return
    #read the header
    self.BytesHead = len(self.__Trj.readline())
    #read the first config
    self.BytesCoords = 0
    for i in range(0, self.__NLine):
      self.BytesCoords += len(self.__Trj.readline())
    #unfortunately for gzip files, we can't seek to end, so we need to
    #do this bit to get the total size:
    self.__Trj.seek(0)
    self.BytesTot = 0
    #try to catch errors in the gzipped file; stride 1024 bytes
    try:
      while len(self.__Trj.read(1024)) > 0: pass
    except IOError:
      #find the exact byte number
      try:
        while len(self.__Trj.read(1)) > 0: pass
      except IOError:
        pass
    self.BytesTot = max(self.__Trj.tell(), 0)
    #set the total number of configs
    self.NCoords = (self.BytesTot - self.BytesHead) / self.BytesCoords
    self.NSkip = min(self.NSkip, self.NCoords)
    a = max(0, self.BytesTot - self.BytesHead - self.BytesCoords*self.NSkip)
    b = self.NStride * self.BytesCoords
    self.SliceNCoords = a / b
    if a % b > 0: self.SliceNCoords += 1
    #set the limit of how many to read in
    if self.NRead is None:
      self.NRead = self.SliceNCoords
    else:
      self.SliceNCoords = min(self.NRead, self.SliceNCoords)


  def __getitem__(self, ind):
    #this index is relative to the sliced version
    if self.__Trj is None:
      raise ValueError
    if type(ind) is int:
      #check for reverse notation
      if ind < 0:
        ind += self.SliceNCoords
      #check bounds
      if ind < 0 or ind >= self.SliceNCoords:
        raise ValueError, "Index out of bounds for trj class."
      else:
        #calculate the absolute index
        self.SliceIndex = ind
        self.Index = self.NSkip + self.NStride * ind
        #seek to the right position
        self.__Trj.seek(self.BytesHead + self.BytesCoords*self.Index)
        #read the data
        s = self.__Trj.read(self.BytesCoords)
        #check to see if we ran out of coordinates
        if len(s) < self.BytesCoords:
          raise IOError
          return
        #parse the crd string
        Pos = ParseCrdString(s, self.AtomNames, self.Mask)
        #update a linked coord array
        if not self.LinkPos is None: self.LinkPos[:,:] = Pos
        return Pos
    else:
      raise ValueError

  def __len__(self):
    return self.SliceNCoords

  def __iter__(self):
    self.Reset()
    return self

  def next(self):
    ind = self.SliceIndex + 1
    if ind < self.SliceNCoords:
      self.Count += 1
      return self[ind]
    else:
      self.Reset()
      raise StopIteration
  

  def Reset(self):
    "Resets current configuration to trajectory start."
    #skip to the right configuration
    self.__Trj.seek(self.BytesHead + self.NSkip * self.BytesCoords)
    #reset the indices
    self.Index, self.SliceIndex = -1, -1
    #set the coordinate set number, corresponding to last set read in
    self.Count = 0

  def GetNextCoords(self, Mask = None):
    "Returns the next set of coordinates in the Trj file, or None if the end is reached."
    OldMask = Mask
    if not Mask is None:
      self.Mask = Mask
    ind = self.SliceIndex + 1
    if ind < self.SliceNCoords:
      self.Count += 1
      result = self[ind]
    else:
      result = None
    Mask = OldMask
    return result

  def GetPastIndices(self):
    "Returns the configuration indices of all read in so far."
    return range(self.NSkip, self.Index+1, self.NStride)

  def Close(self):
    "Closes the open trajectory file."
    if self.__Trj is not None:
      self.__Trj.close()
      self.__Trj = None
    

class PdbListClass:
  "Provides a class for reading successive sets of coordinates from pdb files."
  
  def __init__(self, PdbFileList, Mask = NoMask,
               LinkPos = None):
    """Initializes the class and checks for pdb file existence.
* PdbFileList: list of string names of pdb files
* Mask: list of strings; filter for atom names (default is no mask/empty list)
* LinkPos: an outside array that is updated automatically as coords are read"""
    #check for file existence
    self.PdbFileList = [f for f in PdbFileList if os.path.isfile(f)]
    Diff = len(PdbFileList) - len(self.PdbFileList)
    if Diff > 0:
      raise IOError, "Could not find %d files; removing from list" % Diff
    #set the number of positions
    self.LastLen = -1
    #set the total count
    self.TotalCount = len(self.PdbFileList)
    #set the mask option
    self.Mask = Mask
    #set the linked pos
    self.LinkPos = LinkPos
    #get the sequence, atom names, and atom residues
    f = self.PdbFileList[0]
    self.AtomNames = pdbtools.Atoms(f)
    self.AtomRes = pdbtools.AtomRes(f)
    self.Seq = pdbtools.Seq(f)
    #reset
    self.Reset()
    
  def Reset(self):
    "Resets current configuration to list start."
    #set the coordinate set number, corresponding to last set read in
    self.Count = 0
    self.LastLen = -1
    #set the index
    self.Index = -1

  def GetNextCoords(self, Mask = None):
    "Returns the next set of pdb coordinates, or None if list end is reached."
    if self.Count < len(self.PdbFileList):
      self.Count += 1
      self.Index += 1
      if Mask is None:
        Mask = self.Mask
      Pos = GetPdbCoords(self.PdbFileList[self.Count-1], Mask)
      if self.LastLen > 0 and not self.LastLen == len(Pos):
        raise IOError, "Configuration read with different number of atoms from last read."
      self.LastLen = len(Pos)
      #update the linked pos
      if not self.LinkPos is None: self.LinkPos[:,:] = Pos
      return Pos

  def GetPastIndices(self):
    "Returns the configuration indices of all read in so far."
    return range(0, len(self.PdbFileList))

    
