import sys, os

# Modification history
#
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
            name = "[" + name
        if name[-1] != "]" :
            name = name + "]"
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
  
	# print 'self.atomList',self.atomList
	
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
    

class AtomSelection( object ):
    # parses an atom selection
    def __init__( self, selections, grofile ):
        self.selections = selections
        self.grofile = grofile

        self.atomNames = []
        self.residues = []
        self.includeHydrogen = True
	self.onlyCA = False

        self.lastAtom = -1
        self.lastRes = -1
        self.lastReadRes = -1

        self.parse()
        self.index = self.generateIndex()

    def parse( self ):
        self.name = self.selections[0]
        putativeCommands = self.selections[1:]
        
        for command in putativeCommands :
            if "not hydrogen" in command :
                self.includeHydrogen = False
            if "only CA" in command:
                self.onlyCA = True
            if "residue" in command :
                self.getResNums( command )

    def getResNums( self, command ):
        command = command.replace("residues","")
        command = command.replace("residue","")
        command = command.replace(" ","")

        resNumbers = command.split(",")
        self.residues = []
        for resRange in resNumbers:
            if "-" in resRange:
                resRangeValues = resRange.split("-")
                min = int( resRangeValues[0] )
                max = 1 + int( resRangeValues[1] )
                myRange = range( min, max )
                for resno in myRange: 
                    resno = "%d" % resno
                    self.residues.append( resno.strip() )
            else:
                self.residues.append( resRange.strip() )

    def getAtomNames( self, command ):
        self.atomNames = command.split()[1:]

    def generateIndex( self ):
        FILE = open( self.grofile )
	grofile = FILE.readlines()[2:-1]    # changed from [2:-2] -- it wasn't picking up the last line!!!!!    
        FILE.close()

        atomlist = []        
        name = self.name.strip()

        for line in grofile :
            ( isHydrogen, isCA, resnum, atomNumber ) = self.parseGroLine( line )

            do = False 

            if resnum in self.residues : do = True
            if self.onlyCA == True and isCA == False : do = False
            if self.includeHydrogen == False and isHydrogen == True : do = False

            if do == True : atomlist.append( int( atomNumber ) )

        if atomlist : 
            myNewIndex = IndexGroup( name, atomlist )
        else:
            myNewIndex = None
 
        return myNewIndex

    def parseGroLine( self, line ):
        # gro line looks like:
        #     1NLEU  HD21   17   2.833   2.474   1.648  0.0000  0.0000  0.0000
        # 0123456789012345678901234567890123456789012345678901234567890123456789
        isHydrogen = False
        isCA = False

        # only get res/atom numbers first time
        resnum = ""
        if self.lastRes == -1:
          resnum = line[0:5].strip()
          self.lastReadRes = int(resnum)
          self.lastRes = int(resnum)
        else:
           # if current res# from reading line != lastRes then move to next res
           if self.lastReadRes != int(line[0:5]):
             self.lastReadRes = int(line[0:5])
             resnum = str(self.lastRes + 1)
           else:
             resnum = str(self.lastRes)
        self.lastRes = int(resnum)
        if self.lastAtom == -1:
          atomnum = line[15:20].strip()
        else:
          atomnum = str(self.lastAtom + 1)
        self.lastAtom = int(atomnum)

        atomname = line[9:15].strip()
        if "H" in atomname : isHydrogen = True
        if "CA" in atomname : isCA = True

        return ( isHydrogen, isCA, resnum, atomnum )

def parseSelectionsByGroup( selections ):
    # have a file like:
    # [ nameofgroup1 ]
    # cmd1
    # [ nameofgroup2 ]
    # cmd2
    # cmd3
    # and there might be blank lines
    # just return lists 
    indexGroups = []
    i = -1
    for line in selections:
         if line[0] == "[" :
             # this indicates a new index group
             indexGroups.append( [] )
             i += 1
         if i > -1 :
             indexGroups[ i ].append( line )

    return indexGroups

def getSelections( file, grofile ):
    SELECTIONS = open( file )
    selections = SELECTIONS.readlines()
    SELECTIONS.close()

    # selections file consists of blank lines, lines beginning with these keywords, and      
    # titles. Keywords: 'all', 'not', 'atomname', 'residues'
    # titles begin and end with "["/"]"tmptsYYJC
    # example: 
    # [ helix1-CA ]
    # residues 4-10
    # only CA
    # this little file will produce an index group called 'helix1-CA' which consists of atoms
    # in resiudes 4 through 10 which are called 'CA'

    selections = parseSelectionsByGroup( selections )

    # each element in 'selections' represents an index group. AtomSelection object will parse
    # them.
    
    indexGroups = []

    for selection in selections:
        iGroup = AtomSelection( selection, grofile ).index
        indexGroups.append( iGroup ) 

    return indexGroups

if __name__ == "__main__":
    
    usage = """
    USAGE:  atomSelect.py  selectionsFile  groFile
    
    Displays a set of make_ndx-style index groups for a given *.gro file based on
    the selections made in the selectionsFile  
    
    INPUT
    
    selectionsFile        a file consisting of blank lines, lines beginning with these
		          keywords, and titles. Keywords: 'all', 'not', 'atomname', 'residues'
		          titles begin and end with "[" and "]"
			  
		          example: 
		          [ helix1-CA ]
		          residues 4-10
		          only CA
			  
		          this little file will produce an index group called 'helix1-CA'
			  which consists of atoms in resiudes 4 through 10 which are called 'CA'

    selectionsFile        a file consisting of blank lines, lines beginning with these
		          keywords, and titles. Keywords: 'all', 'not', 'atomname', 'residues'.
    """
    
    if len(sys.argv) < 3:
	print usage
	sys.exit(1)
    
    selectionsFile = sys.argv[1]
    groFile = sys.argv[2]
    
    indexGroups = getSelections( selectionsFile, groFile )
    #print "there are %d index groups" % len( indexGroups )
    for group in indexGroups:
        try:
            group.displayIndexGroup() 
        except:
            pass
 
