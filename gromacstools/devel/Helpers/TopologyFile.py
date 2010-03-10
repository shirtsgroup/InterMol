# TopologyFile.py
#
# HISTORY
# VAV: 1/15/08 - Just created the object to stopre chunks of text.  Later will add functionality

import sys, os, tempfile, string, copy

# Classes 

  
class TopologyFile(object):
    
    def __init__(self, topfile=None, debug=True): 
        """A class to store and modify Gromacs *.top topology files.
        
        INPUT    
            templatefile        a *.top file to read in as a template
        
        OUTPUT
            None
        """
    
        self.debug = debug 
        
        # read infile contents and store as lines
        self.lines = []
        self.includes = []      # raw lines of text, i.e. '#include "ffamber03.itp"'
        self.ifdefs = []        # list of list of raw lines of text for each #ifdef (#ifdef'ed  #includes are in here too)
        self.molecules = []     # a list of TopologyFileMolecule object classes 
        self.systemTitle = ''   # a one-line title for the system
        self.nmols = {}         # a dictionary of {molecule name: #mols}
        
        self.nmolsOrder = []    # It's very important to preserve the *ordering* of the defined [ moleculetype ]
                                # at the end of the *.top file.  This list specifies the ordering of the self.nmols.keys()
        
        if not topfile == None:
            [self.lines, self.includes, self.ifdefs, \
             self.molecules, self.systemTitle, self.nmols, self.nmolsOrder] = self.read(topfile)
    
    
    def read(self, filename, debug=False):
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
        [nmols, nmolsOrder, templines] = self.parse_nmols(templines, debug=debug)
            
        # NOW, the rest of the lines should be [ moleculetype ] blocks
        molecules = self.buildMolecules(templines, debug=debug)    # a list of TopologyFileMolecules() objects
        
        return [lines, includes, ifdefs, molecules, systemTitle, nmols, nmolsOrder]
    
    def write(self, filename):
        """Writes the lines of the topology *.top file"""
        
        self.rebuild_lines()
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
          else:
            i=i+1
    
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
            nmolsOrder        a list of ordered keys
            linesRemaining    the lines left after the #ifdef are parsed out """
       
        i = 0
        nmols = {}
        nmolsOrder = []
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
                nmolsOrder.append(fields[0])
                if i < len(templines):
                  fields = templines.pop(i).split()
                else:
                  fields = []
          i=i+1
        if debug:
          print 'nmols found:\n', nmols 
          print
    
        return [nmols, nmolsOrder, templines]
    
  
    def buildMolecules(self, templines, debug=False):
        """parse (remaining) lines, and return a list of TopologyFileMolecule
        RETURNS:
            molecules        a list of TopologyFileMolecule objects
        """
        
        # If there are no temporary lines to parse, return an empty list of TopologyFileMolecules
        if len(templines) == 0:
            return []
        
        # Strip out any line before '[ moleculetype ]'
        while templines[0].count( '[ moleculetype ]' ) == 0:
            templines.pop(0)
          
        # Go through the lines, looking for the start of each Block, splitting at this point
        Blocks = []       
        for line in templines:
            if line.count('[ moleculetype ]') > 0:
                Blocks.append([])
            Blocks[-1].append(line)
        
        # Debug
        # print '[ moleculetype ] Blocks found:\n'
        # for Block in Blocks:
        #     print string.joinfields(Block)
             
    
        # convert the text chunks into TopologyFileMolecule objects
        molecules = []
        for Block in Blocks:
          molecules.append( TopologyFileMolecule(Block) )
        
        return molecules
    
    def AddMolecule(self, MoleculeObject, NumberOfMolecules):
        """Add a TopologyFileMolecule object to the topology file."""
        
	self.molecules.append(MoleculeObject)
	
	# Add the nmols lines for each molecule
	self.nmols[ MoleculeObject.getMoleculeName() ] = NumberOfMolecules
	self.nmolsOrder.append( MoleculeObject.getMoleculeName() )
	
	
   
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
        for key in self.nmolsOrder:
          self.lines.append('%-16s %s\n'%(key, self.nmols[key]))
    
    
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
    
        
        # Formatting Specification
        #
        # NOTE: The [ moleculetype ] block contain a series of directives, each that start with "[ moleculetype ]".
        #
        # Within the [ moleculetype ] block are a series of **sub-blocks**:
        #
        #     [ atoms ], [ bonds ], [ pairs ], [ dihedrals ], ...
        #
        # Each sub-block consists of sets of lines demarcated by blank lines (i.e. \newlines).  In some cases,
        # there may be multiple sub-blocks that need to collected together.  An example of this is the [ dihedrals ]
        # sub-block -- there is often two subblocks, the first specifying proper dihedrals and the second specifying
        # improper dihedrals.
        #
        #;;;;;;
        # [ moleculetype ]
        # ; Name            nrexcl
        # Protein             3
        # 
        # [ atoms ]
        # ;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
        #      1 amber99_39      1   NPRO      N      1     -0.202      14.01   ; qtot -0.202
        #      2 amber99_17      1   NPRO     H1      2      0.312      1.008   ; qtot 0.11
        #      3 amber99_17      1   NPRO     H2      2      0.312      1.008   ; qtot 0.422
        # ...
        # 
        # [ atoms ]
        # ;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
        # 
        # [ bonds ]
        # ;  ai    aj funct            c0            c1            c2            c3
        #     1     2     1
        #     1     3     1
        # ...
        #
        # [ dihedrals ]
        # ;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
        #
        #####################


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
            
        # Store raw lines for everything in the molecule   
        self.lines = lines
        
        # Create SubBlock objects for each [ moleculetype ], [ atoms ], [ bonds ], [ pairs ], [ angles ], [ dihedrals ]
        # NOTE that each SubBlock object will parse the full set of lines to find its data 
        
        self.subBlocks = {}
        
        self.subBlocks['moleculetype'] = TopologyFileMoleculeSubBlock('moleculetype', inlines=lines,
                 headerFormat='; %-16s %8s ',
                 dataFormat='%-16s %8s ',
                 dataIndexedBy = 'fields[0]' )
        
        self.subBlocks['atoms'] = TopologyFileMoleculeSubBlock('atoms', inlines=lines,
                 headerFormat='; %6s %-10s %6s %6s %6s %6s %-10s %-10s ',
                 dataFormat='%8s %-10s %6s %6s %6s %6s %-10s %-10s ',
                 dataIndexedBy = 'int(fields[0])' )
        
        self.subBlocks['bonds'] = TopologyFileMoleculeSubBlock('bonds', inlines=lines,
                 headerFormat='; %6s %6s %6s ',
                 dataFormat='%8s %6s %6s ',
                 dataIndexedBy = '(int(fields[0]), int(fields[1]))' )
        
        self.subBlocks['pairs'] = TopologyFileMoleculeSubBlock('pairs', inlines=lines,
                 headerFormat='; %6s %6s %6s ',
                 dataFormat='%8s %6s %6s ',
                 dataIndexedBy = '(int(fields[0]), int(fields[1]))' )
        
        self.subBlocks['angles'] = TopologyFileMoleculeSubBlock('angles', inlines=lines,
                 headerFormat='; %6s %6s %6s %6s ',
                 dataFormat='%8s %6s %6s %6s ',
                 dataIndexedBy = '(int(fields[0]), int(fields[1]), int(fields[2]))' )
        
        self.subBlocks['dihedrals'] = TopologyFileMoleculeSubBlock('dihedrals', inlines=lines,
                 headerFormat='; %6s %6s %6s %6s %6s ',
                 dataFormat='%8s %6s %6s %6s %6s ',
                 dataIndexedBy = '(int(fields[0]), int(fields[1]), int(fields[2]), int(fields[3]))' )
        
        
        return
      
  
    # GETTING AND SETTING
    
    def getMoleculeName(self):
	"""Return the name of the molecule, as stored in the [ moleculetype ] data."""
	return self.subBlocks['moleculetype'].data.keys()[0]   # there should be only one index token -- the molecule name
	
    
    def get_atomtype(self, atomnum):
        """Returns the forcefield atom type of the specified atom number.
        Returns None if that atom number does not exist."""
        
        if self.subBlocks['atoms'].data.has_key(atomnum):
          return self.subBlocks['atoms'].data[atomnum][1]
        else:
          return None
      
    def set_atomtype(self, atomnum, atomtype):
        """Set the forcefield atom type of the specified atom number to the given atomtype string."""
        if self.subBlocks['atoms'].data.has_key(atomnum):
            self.subBlocks['atoms'].data[atomnum][1] = atomtype    
  
    def get_atomresnum(self, atomnum):
        """Returns the forcefield atom type of the specified atom number.
        Returns None if that atom number does not exist."""
        
        if self.subBlocks['atoms'].data.has_key(atomnum):
          return int(self.subBlocks['atoms'].data[atomnum][2])
        else:
          return None
  
    def set_atomresnum(self, atomnum, atomresnum):
        """Set the atomresnum (int) of the specified atom number."""
        if self.subBlocks['atoms'].data.has_key(atomnum):
            self.subBlocks['atoms'].data[atomnum][2] = str(atomresnum)
  
          
    def get_atomresname(self, atomnum):
        """Returns the forcefield atom type of the specified atom number.
        Returns None if that atom number does not exist."""
        
        if self.subBlocks['atoms'].data.has_key(atomnum):
          return self.subBlocks['atoms'].data[atomnum][3]
        else:
          return None
  
    def set_atomresname(self, atomnum, atomresname):
        """Set the atomresnum (int) of the specified atom number."""
        if self.subBlocks['atoms'].data.has_key(atomnum):
            self.subBlocks['atoms'].data[atomnum][3] = atomresname
          
    def get_atomname(self, atomnum):
        """Returns the forcefield atom type of the specified atom number.
        Returns None if that atom number does not exist."""
        
        if self.subBlocks['atoms'].data.has_key(atomnum):
          return self.subBlocks['atoms'].data[atomnum][4]
        else:
          return None
    
    def set_atomname(self, atomnum, atomname):
        """Set the atomresnum (int) of the specified atom number."""
        if self.subBlocks['atoms'].data.has_key(atomnum):
            self.subBlocks['atoms'].data[atomnum][4] = atomname
          
    def get_cgnr(self, atomnum):
        """Returns the forcefield atom type of the specified atom number.
        Returns None if that atom number does not exist."""
        
        if self.subBlocks['atoms'].data.has_key(atomnum):
          return int(self.subBlocks['atoms'].data[atomnum][5])
        else:
          return None
  
    def set_cgnr(self, atomnum, cgnr):
        """Set the atomresnum (int) of the specified atom number."""
        if self.subBlocks['atoms'].data.has_key(atomnum):
            self.subBlocks['atoms'].data[atomnum][5] = str(cgnr)
          
    def get_charge(self, atomnum):
        """Returns the forcefield atom type of the specified atom number.
        Returns None if that atom number does not exist."""
        
        if self.subBlocks['atoms'].data.has_key(atomnum):
          return float(self.subBlocks['atoms'].data[atomnum][6])
        else:
          return None
  
    def set_charge(self, atomnum, charge):
        """Set the atomresnum (int) of the specified atom number."""
        if self.subBlocks['atoms'].data.has_key(atomnum):
            self.subBlocks['atoms'].data[atomnum][6] = '%8.6f'%charge
          
    def get_mass(self, atomnum):
        """Returns the forcefield atom type of the specified atom number.
        Returns None if that atom number does not exist."""
        
        if self.subBlocks['atoms'].data.has_key(atomnum):
          return float(self.subBlocks['atoms'].data[atomnum][7])
        else:
          return None
  
    def set_mass(self, atomnum, mass):
        """Set the atomresnum (int) of the specified atom number."""
        if self.subBlocks['atoms'].data.has_key(atomnum):
            self.subBlocks['atoms'].data[atomnum][7] = '%6.4f'%mass
  
    def set_atomnum(self, atomnum, newatomnum):
        """Set the atomnum (int) of the specified atom number."""
        if self.subBlocks['atoms'].data.has_key(atomnum):
            self.subBlocks['atoms'].data[atomnum][0] = str(newatomnum)
    
          
    # ADDING AND REMOVING
    
    def remove_atom(self, atomnum, Renumber=False):
        """Remove the specified atom number and all refernces to it in the bonded parameters.
        If Renumber=True, renumber the remaning atoms."""
        
        # Since the data is indexed different for each SubBlock, let's go through it indiviudually
        for key in self.subBlocks['atoms'].data.keys():
          if (key[0] == atomnum):
            self.subBlocks['atoms'].data.pop(key)
        
        for key in self.subBlocks['bonds'].data.keys():
          if (key[0] == atomnum) | (key[1] == atomnum):
            self.subBlocks['bonds'].data.pop(key)
            
        for key in self.subBlocks['pairs'].data.keys():
          if (key[0] == atomnum) | (key[1] == atomnum):
            self.subBlocks['pairs'].data.pop(key)
            
        for key in self.subBlocks['angles'].data.keys():
          if (key[0] == atomnum) | (key[1] == atomnum) | (key[2] == atomnum):
            self.subBlocks['angles'].data.pop(key)
            
        for key in self.subBlocks['dihedrals'].data.keys():
          if (key[0] == atomnum) | (key[1] == atomnum) | (key[2] == atomnum) | (key[3] == atomnum):
            self.subBlocks['dihedrals'].data.pop(key)

        # Renumber the atoms so they're in order    
        if Renumber:
            self.renumber_atoms()
          
          

    def rebuild_lines(self):
      """Rebuild the self.lines from the data"""
      
      # Rebuild the lines of the constituent SubBlocks
      self.lines = []
      for s in [ 'moleculetype', 'atoms', 'bonds', 'pairs', 'angles', 'dihedrals' ]:
          self.lines += self.subBlocks[s].rebuild_lines()
          
      return self.lines
  
      
    def __repr__(self):
      """Return representation of the TopologyFileMolecule object, in this case,
      rebuilding the lines that would be in the *.top file"""
          
      self.rebuild_lines()
      return string.joinfields(self.lines,'')


class TopologyFileMoleculeSubBlock(object):
    """A data structure object to store the SubBlocks of TopologyFile molecules."""
        
    def __init__(self, name, inlines=None,
                 headerFormat='; %-16s %8s ',
                 dataFormat='%-16s %8s ',
                 dataIndexedBy = '(int(fields[0]), int(fields[1]))'):
        
        """Initializes the TopologyFileMoleculeSubBlock object."""
        
        self.name = name
        self.alllines = inlines   # a list of all the lines that came in for parsind
        self.lines = []           # a list of all lines relevant to this subBlock
        
        self.headerFormat = headerFormat
        self.headerFields = ['' for i in range(headerFormat.count('%'))]
        self.numberHeaderFields = headerFormat.count('%')

        self.dataFormat = dataFormat
        self.numberDataFields = dataFormat.count('%')
        self.dataIndexedBy = dataIndexedBy
        self.data = {}   # A dict of data fields (strings) indexed by atomnum
        
        # If a list of lines were provided at initialization, parse them
        if not inlines == None:
            self.parseLines(self.alllines)
        
    def parseLines(self, lines):
        """Parse the list of lines and store the resulting object data."""
      
        # Go through the lines, looking for a starting lines containing self.name
        # and ending lines that are blank.
        
        StoreData = False
        for line in lines:
            
            # Determine if we should start or stop storing data
            if line.count(self.name) > 0:
                StoreData = True
            if len(line.strip()) == 0:
                StoreData = False
                
            if StoreData:
                
                # Add this line to the relevant list of lines
                self.lines.append(line) 
                
                # Parse data from the lines
                if line.count(self.name) > 0:
                    pass
                elif line.strip()[0] == ';':
                    self.headerFields = line.replace(';','').split()
                else:   # split the line into text fields
                    fields = line.split()   
                    self.data[ eval(self.dataIndexedBy) ] = fields


    def rebuild_lines(self):
        """Rebuild a list of lines representing the SubBlock data from the stored values."""
        
        # If there's no data in this SubBlock, then we don't need to return any lines!
        # if len(self.data) == 0:
        #     return []
        
        # Otherwise, let's rebuild the lines
        self.lines = []
        
        # Append the SubBlock title line
        self.lines.append('[ %s ]\n'%self.name)
        
        # Append the header comment line
        if len(self.headerFields) >= self.numberHeaderFields:
            self.lines.append(self.headerFormat%tuple(self.headerFields[0:self.numberHeaderFields])  \
                           + string.joinfields(self.headerFields[self.numberHeaderFields:],' ')+'\n')
            
        # Append the data lines
        sorted_keys = self.data.keys()
        sorted_keys.sort()
        for key in sorted_keys:  
            if len(self.data[key]) >= self.numberDataFields:
                 self.lines.append(self.dataFormat%tuple(self.data[key][0:self.numberDataFields])  \
                           + string.joinfields(self.data[key][self.numberDataFields:],' ') +'\n')
                 
        # A blank line always         
        self.lines.append('\n')
        return self.lines
        

    def __repr__(self):
        """Return representation of the object, in this case,
        rebuilding the lines that would be in the *.top file"""
            
        self.rebuild_lines()
        return string.joinfields(self.lines,'')
        
        

# TopologyFile.py exceptions
class TopologyFileFormatError(Exception): pass


  
