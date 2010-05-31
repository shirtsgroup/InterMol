#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""

"""

import sys, os, tempfile, string, copy

import simtk.unit as units


from Topology import * 

# Classes 

  
class GromacsTopology(Topology):
    
    """
    DESCRIPTION

    A derived class of Topology.py but with:
    
	1) extra attributes for gromacs topology parameters, #ifdefs, #defines
	2) extra methods for reading and writing to GROMACS *.top files
	
	
    REQUIREMENTS 	

    * Needs an installation of the simtk.units package.  See: https://simtk.org/home/python_units

    * The environment variable $GMXLIB needs to be defined, listing the directory where Gromacs forcefield and *.itp files
      are stored.  This can be added to .bash_profile, and sourced, for example:
      
      export GMXLIB="/home/my-account/local/gromacs-3.1.4/share/top"
        
    TO DO
    
    * Some parameters will be specified using the #included *.itp, while others are specified in the file.
    * Be able to write specified sets of parameters to *.itp files at will. 
    
    NOTES    
    
    """
    
    def __init__(self, topology=None, topfile=None, name=None):
        """
        Create a new GromacsTopology object (a derived Topology class)

        If a GromacsTopology object is specified, it will be queried to construct the class.
        If a *.top filename is specified, it will be read in.

        """
	
	# name of GromacsTopology (will also be used as name of [ system ] in the *.top file)
	self.name = name
	if self.name == None: 
            self.name = "Untitled"         
        
        self.parameters     = list()   # These are GromacsParameter[i] objects, akin to Force objects,
				       # but holding Force parameters for various atomtype combinations 
        self.molecules      = list()   # molecules[i] is the ith TopologySystem object
	
	if (topfile):
	    self.GromacsTopologyFileObject = GromacsTopologyFile(topfile=topfile)
	else:
	    self.GromacsTopologyFileObject = None
	    
	return
    
    def setName(self, name):
	self.name = name
	return
    
		
	
    def readTopologyFile(self, topfile):
	"""Read in the contents of the topology file.
	
	NOTES
	"""
	# Read in the *.top file, creating a new GromacsTopologyFile object
	self.GromacsTopologyFileObject = GromacsTopologyFile(topfile=topfile)
	
	
	# Translate all the parameter directives to GromacsParameter objects
	for parm in self.GromacsTopologyFileObject.ParameterDirectives:
	    self.parameters.append( ConvertParameterDirectiveToGromacsParameter(parm) )

	
	# Translate each set of molecule directives to TopologySystem objects
	for mol in range(len(self.GromacsTopologyFileObject.MoleculeDefinitionDirectives)):
	    self.molecules.append( self.ConvertMoleculeDirectivesToTopologySystem(self.GromacsTopologyFileObject.MoleculeDefinitionDirectives[mol]))	    
    
	return
    
    

    def ConvertMoleculeDirectivesToTopologySystem(self, MoleculeDirectives):
	"""Takes a list of MoleculeDirectives (Directive objects), and uses it
	to create a TopologySystem object withe correct forces."""
	pass # VAV (May 30, 2010):  Needs to be implemented
	
    def ConvertParameterDirectiveToGromacsParameter(self, parmDirective):
	"""Takes a parmDirective (Directive object), and uses it
	to create a GromacsParameter object to store the forcefield parameters."""
	pass # VAV (May 30, 2010):  Needs to be implemented
	
	

    def writeTopologyFile(self, filename, ExpandIncludes=False, RebuildDirectives=False):
	"""Write the Gromacs topology *.top file to file."""
	
	# Rebuild the Directives from the information in the Topology objects.
	if (RebuildDirectives):
	    pass
	
	self.GromacsTopologyFileObject.write(filename, ExpandIncludes=ExpandIncludes)
	return
    
    
	
class GromacsTopologyFile(object):  
    
    def __init__(self, topfile=None, defines=None): 
        """A class to store and modify information in Gromacs *.top topology files.
	
	A Gromacs *.top topology file stores not only the connectivity for each molecule in the system, but (usually as #includes),
	all of the forcefield parameters for all possible (or at least all common) combinations of bonded atom types. 
	
	Some information about default #define parameters, from gmx 3.1.4 manual:
	
	"You can use any defines to control options in your customized topology files. Options that are
	already available by default are:
	
	define = -DFLEX_SPC
	    Will tell grompp to include FLEX SPC in stead of SPC into your topology, this is necessary to make conjugate
	    gradient work and will allow steepest descent to minimize further.
	define = -DPOSRE
	    Will tell grompp to include posre.itp into your topology, used for position restraints.
        """
    	
	if (defines):
	    self.defines = defines
	else:
	    self.defines = list()  # A list of strings that have been #define'd

	# Add default #defines    
	self.addDefine('FLEX_SPC')
	self.addDefine('POSRE')
		
	# *.top file contents will be stored as lines
        self.lines = []
	
        # lines after preprocessing
	self.processed_lines = []
		
	# Store an ordered list of Directive objects, to be created after parsing the lines 
	self.directives = []  

	
	self.ParameterDirectives = list()
	self.MoleculeDefinitionDirectives = list()  # this is a list of lists (for multiple molecules)	
	self.SystemDirectives = list()
	
	# Read in the contents of the topology file, if supplied
	if not topfile == None:
            self.read(topfile)
	    	
	
	"""
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
	"""
    
    

	
	
    def addDefine(self, defineString):
	"""Add a #define string to the list of self.defines."""
	
	if self.defines.count(defineString) == 0:
	    self.defines.append(defineString)
	return 
	
	
	
    def read(self, filename):
        """
	Read in the contents of a Gromacs topology file.
	
	WHAT ARE WE PARSING?
	
	The contents of a *.top file is formatted in of blocks of Directives.  A Directive looks like this: [ my_directive ]	
	From the gmx 3.1.4 manual, here's the description of what to expect:
		
	    * semicolon (;) and newline surround comments
	    * on a line ending with \ the newline character is ignored.
	    * directives are surrounded by [ and ]
	    * the topology consists of three levels:
	      - the parameter level (see Table 5.3)
	      - the molecule level, which should contain one or more molecule definitions (see Table 5.4)
	      - the system level: [ system ], [ molecules ]


	HOW DO WE PARSE THIS?
	
	The basic strategy is straighforward:  For each Directive's block of text, we'll create and ordered list of
	Directive() objects, which are private container classess.  The only wrinkles here are that Gromacs *.top topology 
	file have preprocesser statements:
	
	    #include "myfile.itp"
	
		Inserts a file's contents at a specific position in the *.top file
	
	    #ifdef USE_THIS
	    <statements>
	    #endif
	
		Only reads the statements if the variable USE_THIS is the defined (in the *.mdp, usually)
	    
	To deal with this, we do our *own* preprocessing when we read in the file, in two passes:
	
	    1) We parse all #ifdef statements accoring to our internal list in self.defines (these can be set by the user 
	       before reading the *.top file using methods addDefine() and delDefine()  )
	    2) We consider each #include as a special kind of Directive class called an IncludeDirective, which in turn stores
	       it's own GromacsTopology object (!)  Arbitrary levels of #include-nesting thus be read, and their parameters
	       unpacked into the container classes ParameterInfo or MoleculeDefinitionInfo classes.  However, since the original
	       structure of nested #includes (wel, at least the top level) is still intact, we can write a *.top to file using
	       only the highest-level (i.e. original)( top level.
	    
	There also is a convention for comments that we assume is in place:
	
	    1) Any free-floating "; <comment>" line gets its own Directive (CommentDirective)
	    2) Any "; <comment>"s betwewn a directive statement and directive data is considered to be the "header" for that directive
	
	"""
	
	# Read in the *.top topology file lines
	self.lines = self.readLines(filename)
	
	# Perform pre-processing of the #ifdef and #endif lines 
	self.processedLines = self.preprocess()
	
	# Parse the contents into a list of Directives
	self.parseDirectives()
	
	# Organize the Directives into Parameter, MoleculeDefinition, and System groups
	self.organizeDirectives()
		

    def readLines(self, filename):
	
	fin = open(filename, 'r')
	lines = fin.readlines()
	fin.close()
	return lines


    def preprocess(self):
	
	processedLines = []

        # Parse the #ifdef statements
	includeThisLine = True
        for line in self.lines:
	    
	    if line[0:6] == '#ifdef':
		fields = line.split()
		if len(fields) > 1:
		    defineString = fields[1]
		    if self.defines.count(defineString) > 0:
			includeThisLine = True
		    else: 
			includeThisLine = False
	    elif line[0:6] == '#endif':
		includeThisLine = True

	    else:
		if includeThisLine:
		    processedLines.append(line)
		
        # Concatenate any line continuations 
	i = 0
	while i < len(processedLines)-1:
	    if len(processedLines[i].strip()) > 0:
		if processedLines[i].strip()[-1] == '\\':
		    processedLines[i] = processedLines[i].strip() + processedLines.pop(i+1)
	    i += 1
	     
		
	return processedLines


    def parseDirectives(self):
	"""Parse the lines of the topology file contents into an ordered list of Directive containers.""" 
	
	self.directives = []
	
	index = 0
        while index < len(self.processedLines):
	    
	    line = self.processedLines[index]
	    if len(line) > 0:
		
		# Is the line an include?
		if line.count('#include') > 0 and line[0] == '#':
		    # Then create an IncludeDirective
		    self.directives.append( self.IncludeDirective(line) )

		# Is the line a comment (or a set of comments)?
		elif line[0] == ';':
		    linetxt = ''
		    while self.processedLines[index][0] == ';':
		        linetxt += self.processedLines[index]
			index += 1
			if not (index < len(self.processedLines)):
			    break		    
		    # Then create an CommentDirective
		    self.directives.append( self.CommentDirective(linetxt) )

		# Is the line a define?
		elif line[0] == '#':
		    # Then create an DefineDirective
		    self.directives.append( self.DefineDirective(line) )
		    # Also, add the #define string to our list of defines
		    self.addDefine(line.strip().split()[1])

		# Is the line a blank line?
		elif line.strip() == '':
		    # ignore it
		    pass

		# Is the line the start of a directive?
	        elif line[0] == '[' and line.count(']') > 0:
		    myDirective = self.Directive(line, header='')
		    index += 1
		    # until we hit a blank line (or the end of the file) ...
		    while len(self.processedLines[index].strip()) > 0:
			# comments get concatenated to the header
			if self.processedLines[index][0] == ';':
			    myDirective.header += self.processedLines[index]
			# ... and data lines get appended to Directive.lines.
			else:
			    myDirective.lines.append( self.processedLines[index] ) 
		        index += 1
			
			if not (index < len(self.processedLines)):
			    break
			
		    self.directives.append( myDirective ) 
			

	    # Go to the next line
	    index += 1

    def organizeDirectives(self):
	"""Organize the Directives into Parameter, MoleculeDefinition, and System groups"""
	
	# Make an expanded list of Directives from the nested IncludeDirectives, skipping CommentDirectives and DefineDirectives
	expandedDirectives = self.getAllDirectives([])

	# group the System Directives first
	index = 0
	while index < len(expandedDirectives):
	    if expandedDirectives[index].name.count( '[ system ]' ) > 0:
		self.SystemDirectives.append( expandedDirectives.pop(index) )
	    elif expandedDirectives[index].name.count( '[ molecules ]' ) > 0:
		self.SystemDirectives.append( expandedDirectives.pop(index) )
	    else:
		index += 1
		
	# group the Parameter Directives next
	index = 0
	while index < len(expandedDirectives) and ( expandedDirectives[index].name.count( '[ moleculetype ]' ) == 0 ):
	    self.ParameterDirectives.append( expandedDirectives.pop(index) )

	# finally group lists of lists of molecule directives   
	index = 0
	while index < len(expandedDirectives):
	    if expandedDirectives[index].name.count( '[ moleculetype ]' ) > 0:
		self.MoleculeDefinitionDirectives.append( [] )
	    self.MoleculeDefinitionDirectives[-1].append(expandedDirectives.pop(index))
	    
	return
	    

    def getAllDirectives(self, appendToList, IgnoreComments=True, IgnoreDefines=True):
	"""Get the complete list of directives by traversing all nested IncludeDirectives."""
	
	for d in self.directives:
	    if d.classType == 'Directive':
		appendToList.append(d)
	    elif d.classType == 'IncludeDirective':
		appendToList = d.GromacsTopologyFileObject.getAllDirectives(appendToList)
        return appendToList

    def testDirectives(self):
	"""Print out all the Directives, in order."""
	
	for d in self.directives:
	    print '%%% Directive:', d.classType
	    print d

    def write(self, filename, ExpandIncludes=False):
	"""Write the GromacsTopologyFile to file"""
	
	fout = open(filename, 'w')
	for d in self.directives:
	    if d.classType == 'IncludeDirective':
		if ExpandIncludes:
		    fout.write( repr(d) )
		else:
		    fout.write( d.name + '\n' ) # directives need trailing '\n' for padding
	    else:
		fout.write( repr(d) )
	fout.close()
	
	return
	    
    
    
    #=============================================================================================
    # Containers for *.top file DIRECTIVES, e.g.: [ defaults ], [ atomtypes ]
    #=============================================================================================
    
    class Directive(object):
	"""A base class for for the Directive container"""
	
	def __init__(self, name, header=None):
	
	    self.classType = 'Directive'
	    self.name = name        # Can be 'name', or be '[ name ]\n' (It will be fixed to be the latter).
	    self.header = header    # any "; <comment>\n" text between the name and data (can have multiple newlines).
	    self.lines = []	    	    
	    return
	
	def fixName(self):
	    """if name is 'name', convert it to '[ name ]' """
	    self.name = '[ ' + self.name.replace('[','').replace(']','').strip() + ' ]\n'
	    return

 	def getName(self):
	    return self.name
	
	def setName(self, name):
	    self.name = name
	    return
	    
	def getHeader(self, header):
	    self.header = header
	    return

	def setHeader(self, header):
	    self.header = header
	    return

	def __repr__(self):
	    """Returns a string representation that can be printed or written to file."""
	    
	    outtxt = self.name + self.header
	    for line in self.lines:
		outtxt += line
	    outtxt += '\n'  # each directive needs a trailing newline for coding
	    return outtxt
	    


    #=============================================================================================
    # special DIRECTIVES
    
    class CommentDirective(Directive):
	""" any line that starts with a semicolon """

	def __init__(self, comment):
	
	    self.classType = 'CommentDirective'
	    self.name = comment	    	    
	    return

	def __repr__(self):
	    """Returns a string representation that can be printed or written to file."""
	    return self.name+'\n'  # each directive needs a trailing newline for coding

    class DefineDirective(Directive):
	""" any line that starts with '#define'  """

	def __init__(self, name):
	
	    self.classType = 'DefineDirective'
	    self.name = name	    	    
	    return

	def __repr__(self):
	    """Returns a string representation that can be printed or written to file."""
	    return self.name+'\n'  # each directive needs a trailing newline for coding

	
    class IncludeDirective(Directive):
	"""A special Directive container that points to the filename of the #include file,
	and has as an attribute a self.GromacsTopologyFileObject for that #include file. """
	
	def __init__(self, includeString):
	    self.classType = 'IncludeDirective'
	    self.name  = includeString   # should be entire '#include\n'
	    self.includeFileName = self.name.replace('#include ','').replace('"','').strip()
	    
	    if not os.path.exists(self.includeFileName):
		self.includeFileName = os.path.join(os.environ['GMXLIB'], self.includeFileName)
		
		
	    self.GromacsTopologyFileObject = GromacsTopologyFile(topfile=self.includeFileName)
	    return

	def __repr__(self):
	    """Returns a string representation that can be printed or written to file."""
	    outtxt = ''
	    for d in self.GromacsTopologyFileObject.directives:
		outtxt += repr(d)
	    outtxt += '\n'  # each directive needs a trailing newline for coding
	    return outtxt 



     
