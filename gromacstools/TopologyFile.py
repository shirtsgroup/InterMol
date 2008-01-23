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
    
  def __init__(self, templatefile, debug=False): 
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
    
    [self.lines, self.includes, self.ifdef, self.molecules, self.systemTitle, self.nmols] = self.read(templatefile)


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
    
    # parse #includes
    [includes, templines] = self.parse_includes(templines, debug=debug)
       
    # parse #ifdefs
    [ifdefs, templines] = self.parse_ifdefs(templines, debug=debug)
          
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
  
  
  def parse_includes(self, templines, debug=False):
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
    if debug:
      print 'ifdefs found:\n'
      for ifdef in ifdefs:
        print_linelist(ifdef)
        print
    
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
    (Doesn't do anything right now execpt store all the lines corresponding to the [ moleculetype ]
    """
    
    # raw lines for everything in the molecule   
    self.lines = lines
    
    return
  