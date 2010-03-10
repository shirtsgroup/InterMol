import sys, os

from mmtools.gromacstools.devel.Structure import *

# Modification Histroy:
# VAV:  Feb 8, 2010:    recasting as an IndexFile.py object with methods
#
# Based on atomSelect.py by Dan Ensign
# Modification history of atomSelect.py:
# VAV:  June 12, 2007:  changed AtomSelection.generateIndex() to return None if error (not 0)
# VAV:  June 11, 2007:  Fixed getSelections() so it returns IndexGrop objects, not just their indiex number 
# VAV:  June 7,  2007:  Added __main__ interface, and getIndexGroupText() to return strings (not print to stdout), fixed formatAtomList() (addded spaces)
# VAV:  June 8,  2007:  Fixed grofile-reading bug where the last line wasn't being read!
# GRB:  June 22, 2007:  Change code allow for atom/res #'s with 6+ digits
#
#
# TO DO:
# * make the 'all' keyword works
# * encapsulate in an Objeect file with methods 


class IndexFile( object ):
    """A class for reading, writing, and creating GROMACS *.ndx index files from various atom selections.
    Methods include a way to create atom groups by way of a simple text script."""
    
    
    def __init__(self, ndxfile=None):
	"""Initialize the Index File object."""
	
	self.ndxfile = ndxfile
	self.indexGroups = []   # a list of IndexGroup objects
	
	if self.ndxfile != None:
	    self.readIndexFile(self.ndxfile)

	    
    def readIndexFile(self, ndxfile):
	"""Read in the information from a GROMACS *.ndx index file, allocating a list of IndexGroup 
	objects self.indexGroups.  If self.indexGroups already exists, the atom groups will be added
	to the existing list, avoiding duplicates."""
	
	fin = open(ndxfile,'r')
	lines = fin.readlines()
	fin.close()
	
	# get rid of blank lines at the top
	while lines[0] == '':
	    lines.pop(0)
	    
	# the first line should be  '[ name ]'
	try:
	    while len(lines) > 0:
		name = lines.pop(0).replace('[','').replace(']','').strip()
		atomlist = []
		while lines[0][0] != '[':
		    atomList += [int(s) for s in lines.pop(0).split()]
		# avoid duplicates     
		addMe = True
		for i in range(len(self.indexGroups)):
		    if self.indexGroups[i].name == name:
			if atomList == self.indexGroups[i].atomList:
			    addMe = False
			    print 'WARNING: Index group %s is a redundant duplicate -- NOT ADDING.'%name
			else:
			    print 'Competing index group %s has same name existing index group, but a different atomList!'%name
			    raise IndexFileFormatError		      
		self.indexGroups.append( IndexGroup(name, atomList) )
	except:
	    raise IndexFileFormatError

	
    def writeIndexFile(self, ndxfile, skipEmptyGroups=True):
	"""Write the atom groups to a GROMACS *.ndx index file."""

	fout = open(ndxfile, 'w')
	for group in self.indexGroups:
	    if skipEmptyGroups and (len(group.atomList)==0):
		continue
	    else:
		fout.write(group.getIndexGroupText())
		fout.write('\n') # blank line separating index groups
	fout.close()
    
    
    def addIndexGroup(self, name, atomList):
	"""Add a new index group to the index file object."""
	
        self.indexGroups.append( IndexGroup(name, atomList) )
	
    def removeIndexGroup(self, name):
	"""Remove the index group(s) corresponding to the name supplied."""
	
	for i in range(len(self.indexGroups)):
	    if self.indexGroups[i].name == name:
		self.indexGroups.pop(i)
	
    
    def getAtomListFromGrofile(self, grofile, ResName=None, ResNum=None, AtomName=None, AtomNum=None):
	"""Returns a list of atom indices that meet all of the criteria of the specified quanities:
    
        REQUIRED
	grofile        The *.gro structure file to select atom indices from
	
        OPTIONS
	ResName        a string 'TYR LEU', or a list of strings ['TYR', LEU'].  Left unspecified, will select all (default)
	ResNum         a number 6, or a list of numbers [4,5,6,7,8,9,45].   Left unspecified, will select all (default)
	AtomName       a string 'N CA CB', or a list of strings ['N', CA', 'CB'].  Left unspecified, will select all (default)
	AtomNum        a number 345, or a list of numbers [4,5,6,7,8,9,45].   Left unspecified, will select all (default)
	"""
	
	print 'Filtering atoms matching:'
	if ResName:   print '    ResName  =', ResName
	else:         print '    ResName  = ALL'
	if ResNum:    print '    ResNum   =', ResNum
	else:         print '    ResNum   = ALL'
	if AtomName:  print '    AtomName =', AtomName
	else:         print '    AtomName = ALL'
	if AtomNum:   print '    AtomNum  =', AtomNum
	else:         print '    AtomNum  = ALL'
	print '...working'

	
	g = structureTools.GromacsStructureFromGrofile(grofile)
	
	# get a list of all the atom indices
	all_indices = range(1,g.natoms+1)
	
	# filter the ResName list for our selections
	if ResName != None:
	    if isinstance(ResName, str): ResName = ResName.split()
	    ResName_indices = []
	    for i in range(len(g.atoms)):
		if sum([(g.atoms[i].resname.strip()==res) for res in ResName]):
		       ResName_indices.append( g.atoms[i].atomnum )
	else:
	    ResName_indices = all_indices
	    
	# filter the ResNum list for our selections
	if ResNum != None:
	    if isinstance(ResNum, int): ResNum = [ResNum]
	    ResNum_indices = []
	    for num in ResNum:
		ResNum_indices += [a.atomnum for a in g.atoms if a.resnum==num]
	else:
	    ResNum_indices = all_indices
	    
	# filter the AtomName list for our selections
	if AtomName != None:
	    if isinstance(AtomName, str): AtomName = AtomName.split()
	    AtomName_indices = []
	    for name in AtomName:
		AtomName_indices += [a.atomnum for a in g.atoms if a.atomname.strip()==name]
	else:
	    AtomName_indices = all_indices
	    
	# filter the ResName list for our selections
	if AtomNum != None:
	    if isinstance(AtomNum, int): AtomNum = [AtomNum]
	    AtomNum_indices = []
	    for num in AtomNum:
		AtomNum_indices += [a.atomnum for a in g.atoms if a.atomnum.strip()==num]
	else:
	    AtomNum_indices = all_indices
	
	# finding the intersection of all lists...
	IntersectionOfAllLists = list(set(all_indices) & set(ResName_indices) & set(ResNum_indices) \
				      & set(AtomName_indices) & set(AtomNum_indices))
	print '...Done'
	
	return IntersectionOfAllLists
	

class IndexGroup( object ):
    """ class for one index group in an index file"""
    
    def __init__( self, name, atomList, width = 80, precision = 6 ):
        self.name = self.formatName( name )
        self.atomList = atomList
        self.width = width
        self.precision = precision

    def formatName( self, name ):
        name = name.strip()
        if name[0] != "[" :
            name = "[ " + name
        if name[-1] != "]" :
            name = name + " ]"
        return name 

    def formatAtomList( self ):
        # self.atom list is just a list of numbers; this will format them 
        # into a string 'width' characters wide (this is just a guideline, 
        # though -- we'll assume that the strings are 5 characters (precision
        # 10000) so we'll print width/5 atoms in a row

	width = self.width              # 80 chars is the default here
        precision = self.precision      # number of characters in the atom number (6 is a safe bet -- up to 999999)

        formattedList=""
	numberInRow = width / precision
	buffer = ''
  
	# sort the list, in case it's not
	self.atomList.sort()
	
        for atom in self.atomList:          
	    formattingString = "%%-%ds " % precision
            atomString = formattingString % str(atom)
            buffer += atomString
            if len(buffer) > width:
                formattedList += buffer + "\n"
		buffer = ''
	# dump the rest of the buffer
	formattedList += buffer + '\n'	
        return formattedList     
            
    def displayIndexGroup( self ):
        print self.name
        print self.formatAtomList()

    def getIndexGroupText( self ):
	return self.name + '\n' + self.formatAtomList() + '\n'
    


# IndexFile.py exceptions
class IndexFileFormatError(Exception): pass


