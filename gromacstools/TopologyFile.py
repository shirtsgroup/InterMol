# TopologyFile.py
#
# HISTORY
# VAV: 1/15/08 - Just created the object to stopre chunks of text.  Later will add functionality

import sys, os, tempfile, string, copy

# Functions

def print_linelist(linelist):
  """Helper function to print lists of lines as \n'ed blocks"""
  for line in linelist:
    print line.replace('\n','')


# Classes 

  
class TopologyFile(object):
    
  def __init__(self, templatefile, debug=True): 
    """A class to store and modify Gromacs *.top topology files.
    
    INPUT    
        templatefile        a *.top file to read in as a template
    
    OUTPUT
        None
    """

    self.debug = debug 
    # process infile name
    self.templatefile = templatefile
    
    # read infile contents and store as lines
    self.lines = []
    self.includes = []      # raw lines of text, i.e. '#include "ffamber03.itp"'
    self.ifdefs = []        # list of list of raw lines of text for each #ifdef (#ifdef'ed  #includes are in here too)
    self.molecules = []     # a list of TopologyFileMolecule object classes 
    self.systemTitle = ''   # a one-line title for the system
    self.nmols = {}         # a dictionary of {molecule name: #mols}
    
    [self.lines, self.includes, self.ifdefs, self.molecules, self.systemTitle, self.nmols] = self.read(templatefile)


  def read(self, filename, debug=True):
    """Reads the lines of the template *.top, and returns
    includes, molecules, systemTitle, and nmols"""
    
    if debug:
      print 'Reading topology file', filename
    
    # read lines
    fin = open(filename, 'r')
    lines = fin.readlines()
    fin.close()
    
    # keep lines as the full list of top file lines, but parse out pop()s from templines
    templines = copy.copy(lines)

    # parse #ifdefs  (must be done before #includes, because #ifdefs may have includes --VV)
    [ifdefs, templines] = self.parse_ifdefs(templines, debug=debug)
   
    # parse #includes
    [includes, templines] = self.parse_includes(templines, debug=debug)
       
          
    # parse the [ system ] title
    [systemTitle, templines] = self.parse_systemTitle(templines, debug=debug)
    
    # parse the [ molecules ] title to get nmols counts
    [nmols, templines] = self.parse_nmols(templines, debug=debug)
        
    # NOW, the rest of the lines should be [ moleculetypes ]
    molecules = self.buildMolecules(templines, debug=debug)    # a list of TopologyFileMolecules() objects
    
    return [lines, includes, ifdefs, molecules, systemTitle, nmols]

  def write(self, filename):
    """Writes the lines of the topology *.top file"""
    fout = open(filename, 'w')
    for line in self.lines:
      fout.write(line)
    fout.close()
  
  
  def parse_includes(self, templines, debug=True):
    """parses the #include lines in the topology file
    RETURNS  [includes, linesRemaining]:
        include            the #include lines
        linesRemaining     the lines left after the #includes are parsed out """

    includes = []
    i = 0
    while i < len(templines):
      if len(templines[i]) >= len('#include '):
        if templines[i][0:len('#include ')] == '#include ':
          includes.append(templines.pop(i))
        else:
          i=i+1
      else:
        i=i+1
    if debug:
      print 'includes found:\n', print_linelist(includes)
      print

    return [includes, templines]
  
  def parse_ifdefs(self, templines, debug=False):
    """parses the #ifdef lines in the topology file
    RETURNS  [ifdefs, linesRemaining]:
        ifdefs            the #ifdef lines
        linesRemaining    the lines left after the #ifdef are parsed out """


    ifdefs = []
    i = 0
    while i < len(templines):
      print 'templines[',i,']',templines
      if len(templines[i]) >= len('#ifdef '):
        if templines[i][0:len('#ifdef ')] == '#ifdef ':
          ifdefs.append([])
          while templines[i][0:len('#endif')] != '#endif': 
            ifdefs[-1].append(templines.pop(i))
          ifdefs[-1].append(templines.pop(i))  # this should be the #endif line
        else:
          i=i+1
      else:
        i=i+1
    
    return [ifdefs, templines]


  def parse_systemTitle(self, templines, debug=False):
    """parses the systemTitle in the topology file
    RETURNS [systemTitle, linesRemaining]:
        systemTitle       the system title string
        linesRemaining    the lines left after the #ifdef are parsed out """
   
    i = 0
    systemTitle = ''
    while i < len(templines):
      if len(templines[i]) >= len('[ system ]'):
        if templines[i][0:len('[ system ]')] == '[ system ]':
          templines.pop(i)
          # skip comments
          while templines[i][0] == ';':
            templines.pop(i)
          systemTitle = templines.pop(i)
        else:
          i=i+1         
      else:
        i=i+1
        
    if debug:
      print 'systemTitle found:\n', systemTitle 
      print

    return [systemTitle, templines]


  def parse_nmols(self, templines, debug=False):
    """parses the number of mols lines in the topology file
    RETURNS [nmols, linesRemaining]:
        nmols             dictionary of {name:#mols}
        linesRemaining    the lines left after the #ifdef are parsed out """
   
    i = 0
    nmols = {}
    while i < len(templines):
      if len(templines[i]) >= len('[ molecules ]'):
        if templines[i][0:len('[ molecules ]')] == '[ molecules ]' :
          templines.pop(i)
          # skip comments
          while templines[i][0] == ';':
            templines.pop(i)
          fields = templines.pop(i).split()
          while len(fields) >= 2:
            nmols[fields[0]] = fields[1]
            if i < len(templines):
              fields = templines.pop(i).split()
            else:
              fields = []
      i=i+1
    if debug:
      print 'nmols found:\n', nmols 
      print

    return [nmols, templines]


  def buildMolecules(self, templines, debug=False):
    """parse (remaining) lines, and return a list of TopologyFileMolecule
    RETURNS:
        molecules        a list of TopologyFileMolecule objects
    """
    
    # compile a list of line chunks
    
    ### get rid of leading lines
    foundone = False
    while not foundone:
      if len(templines[0]) >= len('[ moleculetype ]'):
        if templines[0][0:len('[ moleculetype ]')] == '[ moleculetype ]':
          foundone = True
        else:
          templines.pop(0)
      else:
          templines.pop(0)

    chunks = []
    while len(templines) > 0:
      if len(templines[0]) >= len('[ moleculetype ]'):
        if templines[0][0:len('[ moleculetype ]')] == '[ moleculetype ]':
          chunks.append([])
          chunks[-1].append( templines.pop(0) )
        else:
          chunks[-1].append( templines.pop(0) )
      else:
        chunks[-1].append( templines.pop(0) )

    if debug:
      print 'chunks found:\n'
      for chunk in chunks:
         print print_linelist(chunk)
         print

    # convert the text chunks into TopologyFileMolecule objects
    molecules = []
    for chunk in chunks:
      molecules.append( TopologyFileMolecule(chunk) )

    if debug:
      print 'moleculetypes found:\n', molecules 
      print


    return molecules

  
  
  def rebuild_lines(self):
    """Rebuild the self.lines from the data"""

    self.lines = []
    for line in self.includes:
      self.lines.append(line)
    self.lines.append('\n')
    
    for ifdef in self.ifdefs:
      for line in ifdef:
        self.lines.append(line)
      self.lines.append('\n')

    for molecule in self.molecules:
      self.lines.append(repr(molecule))
    self.lines.append('\n')

    self.lines.append('[ system ]\n')
    self.lines.append('; Name\n')
    self.lines.append(self.systemTitle)
    self.lines.append('\n')
    
    self.lines.append('[ molecules ]\n')
    self.lines.append('; Compound        #mols\n')
    for key,value in self.nmols.iteritems():
      self.lines.append('%-16s %s\n'%(key, value))
    self.lines.append('\n')
     

  
  def __repr__(self):
    """returns a (rebuilt) string representing the contents of the topology file"""
    
    self.rebuild_lines()
    return string.joinfields(self.lines,'')
  
  

  
  
#####################

class TopologyFileMolecule(object):
  """A data structure object to store TopologyFile molecules:
      [ moleculetype ]
      [ atoms ]
      [ bonds ]
      [ pairs ]
      [ angles ]
      [ dihedrals ]      
  """
  
  
  def __init__(self, lines):
    """Initializes the TopologyFileMolecule object.

    Here's a template to remind ourselves:
    
[ moleculetype ]
; Name            nrexcl
Protein             3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
1 amber99_34      1    ALA      N      1  -0.404773      14.01   ; qtot -0.4048
2 amber99_17      1    ALA      H      2   0.294276      1.008   ; qtot -0.1105
3 amber99_11      1    ALA     CA      3  -0.027733      12.01   ; qtot -0.1382
4 amber99_19      1    ALA     HA      4   0.120802      1.008   ; qtot -0.01743
5 amber99_11      1    ALA     CB      5  -0.229951      12.01   ; qtot -0.2474
6 amber99_18      1    ALA    HB1      6   0.077428      1.008   ; qtot -0.17
7 amber99_18      1    ALA    HB2      7   0.077428      1.008   ; qtot -0.09252
8 amber99_18      1    ALA    HB3      8   0.077428      1.008   ; qtot -0.01509
9  amber99_2      1    ALA      C      9   0.570224      12.01   ; qtot 0.5551
10 amber99_41      1    ALA      O     10  -0.555129         16   ; qtot 0

[ bonds ]
;  ai    aj funct            c0            c1            c2            c3
1     2     1
1     3     1
3     4     1
3     5     1
3     9     1
5     6     1
5     7     1
5     8     1
9    10     1

[ pairs ]
;  ai    aj funct            c0            c1            c2            c3
1     6     1
1     7     1
1     8     1
1    10     1
2     4     1
2     5     1
2     9     1
4     6     1
4     7     1
4     8     1
4    10     1
5    10     1
6     9     1
7     9     1
8     9     1

[ angles ]
;  ai    aj    ak funct            c0            c1            c2            c3
2     1     3     1
1     3     4     1
1     3     5     1
1     3     9     1
4     3     5     1
4     3     9     1
5     3     9     1
3     5     6     1
3     5     7     1
3     5     8     1
6     5     7     1
6     5     8     1
7     5     8     1
3     9    10     1

[ dihedrals ]
;  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5
2     1     3     4     3
2     1     3     5     3
2     1     3     9     3
1     3     5     6     3
1     3     5     7     3
1     3     5     8     3
4     3     5     6     3
4     3     5     7     3
4     3     5     8     3
9     3     5     6     3
9     3     5     7     3
9     3     5     8     3
1     3     9    10     3
4     3     9    10     3
5     3     9    10     3

    """
        
    # raw lines for everything in the molecule   
    self.lines = lines

    self.moleculetype_headers = []
    self.moleculetype = {}      # {Name, nrexcl}

    # dictionaries for the data 
    # NOTE that all the data is stored as strings, but get_ and set_ as values
    self.atoms_headers = []
    self.atoms = {}           # a list of atom dictionaries {atomnum: [list]}     
    self.bonds_headers = []
    self.bonds = {}           # a list of bonds dictionaries {(ai,aj): [list]}     
    self.pairs_headers = []
    self.pairs = {}           # a list of pairs dictionaries {(ai,aj): [list]}     
    self.angles_headers = []
    self.angles = {}          # a list of angles dictionaries {(ai,aj,ak): [list]}     
    self.dihedrals_headers = []
    self.dihedrals = {}          # a list of dihedrals dictionaries {(ai,aj,ak,al): [list]}     
    
    self.parse(self.lines)

    return
  

  def parse(self, lines):
    """parse to the lines of the topology files to the substituent objects"""

    i = 0
    while i < len(lines):
      
      if lines[i].strip() == '[ moleculetype ]':
          i+=1
          while len(lines[i].strip()) > 0:
              # header comment
              if lines[i][0] == ';':
                self.moleculetype_headers = lines[i].replace(';','').split()
              else:
                fields = lines[i].split()
                name = fields[0]
                self.moleculetype[name] = fields  # keep the atom number in the fields (so they match up with the headers)
              i+=1

      
      if lines[i].strip() == '[ atoms ]':
          i+=1
          while len(lines[i].strip()) > 0:
              # header comment
              if lines[i][0] == ';':
                self.atoms_headers = lines[i].replace(';','').split()
              else:
                fields = lines[i].split()
                atomnum = int(fields[0])
                self.atoms[atomnum] = fields  # keep the atom number in the fields (so they match up with the headers)
              i+=1
           
      
      if lines[i].strip() == '[ bonds ]':
          i+=1
          while len(lines[i].strip()) > 0:
              # header comment
              if lines[i][0] == ';':
                self.bonds_headers = lines[i].replace(';','').split()
              else:
                fields = lines[i].split()
                bondnums = (int(fields[0]),int(fields[1]))
                self.bonds[bondnums] = fields  # keep the atom number in the fields (so they match up with the headers)
              i+=1
      
      if lines[i].strip() == '[ pairs ]':
          i+=1
          while len(lines[i].strip()) > 0:
              # header comment
              if lines[i][0] == ';':
                self.pairs_headers = lines[i].replace(';','').split()
              else:
                fields = lines[i].split()
                pairnums = (int(fields[0]),int(fields[1]))
                self.pairs[pairnums] = fields  # keep the atom number in the fields (so they match up with the headers)
              i+=1

      
      if lines[i].strip() == '[ angles ]':
          i+=1
          while len(lines[i].strip()) > 0:
              # header comment
              if lines[i][0] == ';':
                self.angles_headers = lines[i].replace(';','').split()
              else:
                fields = lines[i].split()
                anglenums = (int(fields[0]),int(fields[1]),int(fields[2]))
                self.angles[anglenums] = fields  # keep the atom number in the fields (so they match up with the headers)
              i+=1

      
      if lines[i].strip() == '[ dihedrals ]':
          i+=1
          while len(lines[i].strip()) > 0:
              # header comment
              if lines[i][0] == ';':
                self.dihedrals_headers = lines[i].replace(';','').split()
              else:
                fields = lines[i].split()
                dihedralnums = (int(fields[0]),int(fields[1]),int(fields[2]),int(fields[3]))
                self.dihedrals[dihedralnums] = fields  # keep the atom number in the fields (so they match up with the headers)
              i+=1

              
      i+=1  # advance if no key words found...
      

  # GETTING AND SETTING
  
  def get_atomtype(self, atomnum):
    """Returns the forcefield atom type of the specified atom number.
    Returns None if that atom number does not exist."""
    
    if self.atoms.has_key(atomnum):
      return self.atoms[atomnum][1]
    else:
      return None
    
  def set_atomtype(self, atomnum, atomtype):
    """Set the forcefield atom type of the specified atom number to the given atomtype string."""
    if self.atoms.has_key(atomnum):
        self.atoms[atomnum][1] = atomtype    

  def get_atomresnum(self, atomnum):
    """Returns the forcefield atom type of the specified atom number.
    Returns None if that atom number does not exist."""
    
    if self.atoms.has_key(atomnum):
      return int(self.atoms[atomnum][2])
    else:
      return None

  def set_atomresnum(self, atomnum, atomresnum):
    """Set the atomresnum (int) of the specified atom number."""
    if self.atoms.has_key(atomnum):
        self.atoms[atomnum][2] = str(atomresnum)

        
  def get_atomresname(self, atomnum):
    """Returns the forcefield atom type of the specified atom number.
    Returns None if that atom number does not exist."""
    
    if self.atoms.has_key(atomnum):
      return self.atoms[atomnum][3]
    else:
      return None

  def set_atomresname(self, atomnum, atomresname):
    """Set the atomresnum (int) of the specified atom number."""
    if self.atoms.has_key(atomnum):
        self.atoms[atomnum][3] = atomresname
        
  def get_atomname(self, atomnum):
    """Returns the forcefield atom type of the specified atom number.
    Returns None if that atom number does not exist."""
    
    if self.atoms.has_key(atomnum):
      return self.atoms[atomnum][4]
    else:
      return None

  def set_atomname(self, atomnum, atomname):
    """Set the atomresnum (int) of the specified atom number."""
    if self.atoms.has_key(atomnum):
        self.atoms[atomnum][4] = atomname
        
  def get_cgnr(self, atomnum):
    """Returns the forcefield atom type of the specified atom number.
    Returns None if that atom number does not exist."""
    
    if self.atoms.has_key(atomnum):
      return int(self.atoms[atomnum][5])
    else:
      return None

  def set_cgnr(self, atomnum, cgnr):
    """Set the atomresnum (int) of the specified atom number."""
    if self.atoms.has_key(atomnum):
        self.atoms[atomnum][5] = str(cgnr)
        
  def get_charge(self, atomnum):
    """Returns the forcefield atom type of the specified atom number.
    Returns None if that atom number does not exist."""
    
    if self.atoms.has_key(atomnum):
      return float(self.atoms[atomnum][6])
    else:
      return None

  def set_charge(self, atomnum, charge):
    """Set the atomresnum (int) of the specified atom number."""
    if self.atoms.has_key(atomnum):
        self.atoms[atomnum][6] = '%8.6f'%charge
        
  def get_mass(self, atomnum):
    """Returns the forcefield atom type of the specified atom number.
    Returns None if that atom number does not exist."""
    
    if self.atoms.has_key(atomnum):
      return float(self.atoms[atomnum][7])
    else:
      return None

  def set_mass(self, atomnum, mass):
    """Set the atomresnum (int) of the specified atom number."""
    if self.atoms.has_key(atomnum):
        self.atoms[atomnum][7] = '%6.4f'%mass

        
  # ADDING AND REMOVING
  
  def remove_atom(self, atomnum, Renumber=False):
    """Remove the specified atom number and all refernces to it in the bonded parameters.
    If Renumber=True, renumber the remaning atoms."""
    
    self.atoms.pop(atomnum)

    for key in self.bonds.keys():
      if (key[0] == atomnum) | (key[1] == atomnum):
        self.bonds.pop(key)
        
    for key in self.pairs.keys():
      if (key[0] == atomnum) | (key[1] == atomnum):
        self.pairs.pop(key)
        
    for key in self.angles.keys():
      if (key[0] == atomnum) | (key[1] == atomnum) | (key[2] == atomnum):
        self.angles.pop(key)
        
    for key in self.dihedrals.keys():
      if (key[0] == atomnum) | (key[1] == atomnum) | (key[2] == atomnum) | (key[3] == atomnum):
        self.dihedrals.pop(key)
        
    if Renumber:
      self.renumber_atoms()
      
  def renumber_atoms(self):
    """Renumber the atoms in consecutive order from 1 to N"""
    pass
  
        
  # FORMATTING        
        
  def moleculetypes_lines(self):
    """Build the [ moleculetype ] lines from the data"""
    
    lines = []
    lines.append('[ moleculetype ]\n')
    lines.append('; %-16s %8s '%tuple(self.moleculetype_headers[0:2]) +'\n')
    for key, value in self.moleculetype.iteritems():  # there should only be one key, the molecule name
        lines.append('%-16s %8s '%tuple(value[0:2])+'\n')
    lines.append('\n')
    return lines

  
  def atoms_lines(self):
    """Build the [ atoms ] lines from the data"""
    
    lines = []
    lines.append('[ atoms ]\n')
    lines.append('; %6s %-10s %6s %6s %6s %6s %-10s %-10s '%tuple(self.atoms_headers[0:8])  \
                       + string.joinfields(self.atoms_headers[8:],' ')+'\n')
    sorted_keys = self.atoms.keys()
    sorted_keys.sort()
    for key in sorted_keys:  
        lines.append('%8s %-10s %6s %7s %6s %6s %-10s %-10s '%tuple(self.atoms[key][0:8])  \
                       + string.joinfields(self.atoms[key][8:],' ') +'\n')
    lines.append('\n')
    return lines

  
  def bonds_lines(self):
    """Build the [ bonds ] lines from the data"""
    
    lines = []
    lines.append('[ bonds ]\n')
    lines.append('; %6s %6s %6s '%tuple(self.bonds_headers[0:3]) + string.joinfields(self.bonds_headers[3:],' ')+'\n')
    sorted_keys = self.bonds.keys()
    sorted_keys.sort()
    for key in sorted_keys:  
        lines.append('%8s %6s %6s '%tuple(self.bonds[key][0:3]) + string.joinfields(self.bonds[key][3:],' ') +'\n')
    lines.append('\n')
    return lines

  def pairs_lines(self):
    """Build the [ pairs ] lines from the data"""
    
    lines = []
    lines.append('[ pairs ]\n')
    lines.append('; %6s %6s %6s '%tuple(self.pairs_headers[0:3]) + string.joinfields(self.pairs_headers[3:],' ')+'\n')
    sorted_keys = self.pairs.keys()
    sorted_keys.sort()
    for key in sorted_keys:  
        lines.append('%8s %6s %6s '%tuple(self.pairs[key][0:3]) + string.joinfields(self.pairs[key][3:],' ') +'\n')
    lines.append('\n')
    return lines

  def angles_lines(self):
    """Build the [ angles ] lines from the data"""
    
    lines = []
    lines.append('[ angles ]\n')
    lines.append('; %6s %6s %6s %6s '%tuple(self.angles_headers[0:4]) + string.joinfields(self.angles_headers[4:],' ')+'\n')
    sorted_keys = self.angles.keys()
    sorted_keys.sort()
    for key in sorted_keys:  
         lines.append('%8s %6s %6s %6s '%tuple(self.angles[key][0:4]) + string.joinfields(self.angles[key][4:],' ') +'\n')
    lines.append('\n')
    return lines

  def dihedrals_lines(self):
    """Build the [ angles ] lines from the data"""
    
    lines = []
    lines.append('[ dihedrals ]\n')
    lines.append('; %6s %6s %6s %6s %6s '%tuple(self.dihedrals_headers[0:5]) + string.joinfields(self.dihedrals_headers[5:],' ')+'\n')
    sorted_keys = self.dihedrals.keys()
    sorted_keys.sort()
    for key in sorted_keys:  
         lines.append('%8s %6s %6s %6s %6s '%tuple(self.dihedrals[key][0:5]) + string.joinfields(self.dihedrals[key][5:],' ') +'\n')
    lines.append('\n')
    return lines
  
  def rebuild_lines(self):
    """Rebuild the self.lines from the data"""
    
    self.lines = []
    self.lines = self.lines + self.moleculetypes_lines()
    self.lines = self.lines + self.atoms_lines()
    self.lines = self.lines + self.bonds_lines()
    self.lines = self.lines + self.pairs_lines()
    self.lines = self.lines + self.angles_lines()
    self.lines = self.lines + self.dihedrals_lines()


    
    
    
  def __repr__(self):
    """Return representation of the TopologyFileMolecule object, in this case,
    rebuilding the lines that would be in the *.top file"""
        
    self.rebuild_lines()
    return string.joinfields(self.lines,'')
    

  