from os.path import exists
# from numpy import *	# do I need this?
from os import system,unlink

# version 2-27-2007 3:42 pm (though not changed from DATE variable)
### global variables ###
AUTHOR="Dan Ensign"
VERSION="0.1"
DATE="27 February 2007"

### to do ###
# - right now, getfield.py leaves a shitload of huge .gro files on the disk, probably from
#   xtc2gro. I *think* the unlink I added stopped this, but I need to be sure
# - probably, xtc2gro should be turned back into a function that returns a list of GromacsStructure
#   objects. A second function -- getsnapshot -- should be added which gets the ith snapshot
#   from the xtc file. This is still insanely slow -- a direct (C?) xtc reader would be better. Also
#   add a third function, which just returns a list of gro files in text form. For the last, make
#   certain that these structures can be converted to GromacsStructure objects. Three ways to read
#   xtc files add a nice level of control.
# - add a "trimmer" to the GromacsStructure, useful for getting rid of e.g., water
# - add an index file writer superior to make_ndx (ie, that recognizes "NLEU" is an amino acid)
# - read charges into a GromacsStructure from a .top file
# - find some way to make representation of big structures more efficient!!
# - allow two ways of building a GromacsStructure object from file:
#  	gro=GromacsStructure() ; gro.load("conf.gro")
#		OR
#	gro=GromacsStructure("conf.gro")

# three letter amino acid code to one letter amino acid code
tlc2olc = {
"ALA" : "A" ,
"CYS" : "C" ,
"CYN" : "C" ,
"CYX" : "C" ,
"ASP" : "D" ,
"GLU" : "E" ,
"PHE" : "F" ,
"CPHE" : "F" ,
"GLY" : "G" ,
"HIS" : "H" ,
"HIE" : "H" ,
"HIP" : "H" ,
"ILE" : "I" ,
"NLE" : "J" ,
"LYS" : "K" ,
"LYP" : "K" ,
"LEU" : "L" ,
"NLEU" : "L" ,
"MET" : "M" ,
"ASN" : "N" ,
"PRO" : "P" ,
"GLN" : "Q" ,
"ARG" : "R" ,
"SER" : "S" ,
"THR" : "T" ,
"VAL" : "V" ,
"TRP" : "W" ,
"TYR" : "Y" 
} 

### classes ###
"""
class types: 
------------
- GromacsAtom		- represents one atom. Now just for .gro but will later include topology, trajectory, etc.
- GromacsAtoms		- a fancy list whose elements are all GromacsAtom objects
- GromacsStructure	- represents the information in  a gro file
"""

class GromacsAtom:
	"""A type which contains all the known information about one atom in a gromacs structure,
	topology, or trajectory. (For now just works with the gro structure file.)

	attributes:
	-----------
	- resnum
	- resname
	- atomname
	- atomnum
	- x, y, and z
	- vx, vy, and vz

	methods:
	--------
	- __copy__()	- makes a copy of the atom (TEST!)
	- position()	- returns the position as a 3-d tuple (make this work on the velocity also)
	- printgroatom()- prints the atom as it's formatted in the gro file (used by the GromacsStructure type)
	"""
	
	def __init__(self, resnum=0, resname="", atomname="", atomnum=0, x=0.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0):
		self.resnum=resnum
		self.resname=resname
		self.atomname=atomname
		self.atomnum=atomnum
		self.x=x
		self.y=y
		self.z=z
		#self.r=position()
		self.vx=vx
		self.vy=vy
		self.vz=vz 

	def __copy__(self):
                return GromacsAtom(self.resnum, self.resname, self.atomname, self.atomnum, self.x, self.y, self.z)
	def position(self):
		# is there a way to make this return the velocity vector as well?
		return (self.x, self.y, self.z)

	def printgroatom(self):
		"""prints the atom as it would be formatted in the .gro file"""
		resnum="%5d" % self.resnum
		resname="%-4s" % self.resname # should be left justified
		atomname="%6s" % self.atomname
		atomnum="%5d" % self.atomnum
		pos="% 8.3f% 8.3f% 8.3f" % (self.x, self.y, self.z)
		vel="% 8.4f% 8.4f% 8.4f" % (self.vx, self.vy, self.vz)
		return resnum+resname+atomname+atomnum+pos+vel

class GromacsAtoms( list ):
	"""This is a special list of composed of GromacsAtom objects. There might be certain
	operations that make more sense doing them to all the atoms (say, like printing them).
	Individual atoms are accessed through this object's list elements.

	attributes:
	-----------
	none

	methods:
	--------
	- __copy__()	- makes a copy of the atom list (TEST!)
	- printatomlist() - uses the GromacsAtom printgroatom method to print all the atoms in the object
	- sequence()
	- store()	- takes a gro line and appends a GromacsAtom to atoms
	- test()	- write this!

	"""

	def __init__(self):
		list.__init__(self)	

        def __copy__(self): return GromacsAtoms(self)

        def printatomlist(self):
		# HERE
                for atom in self:
                        atom.printgroatom()
                return

	def sequence(self):
		"""Returns the sequence (all the reside names) of the atom list."""
		sequence=[] 
		oldresnum=0
		for atom in self:
			resname=atom.resname
			resnum=atom.resnum
			if oldresnum != resnum:
				sequence.append(resname)
				oldresnum=resnum
		return sequence

	def store(self,line):
		# 'line' is an atom line in a .gro file. 
		def makefloat(r): # for making a list of float strings into a list of floats
                        try: rf=[float(elem) for elem in r]
                        except:
                                # raise BadCoordinatesInGro
                                print "Bad coordinates in gro!"
                                rf=[0.0,0.0,0.0]
                        return rf

		line.strip()

		try: resnum=int(line[0:5])
		except: raise BadResidueNumber

		resname=line[5:9]
		atomname=line[9:15]

		try: atomnum=int(line[15:20])
		except: raise BadAtomNumber 

		(x,y,z)=makefloat(line[20:44].split())
		try: (vx,vy,vz)=makefloat(line[44:68].split())
		except: (vx,vy,vz)=(0.0,0.0,0.0) # if there aren't any velocities yet

		atom=GromacsAtom(resnum,resname,atomname,atomnum,x,y,z,vx,vy,vz)
		self.append(atom)
		#print "in store method, the number of atoms is",len(self)
		return
	
class GromacsStructure:
	# immutable type? Would this be safer/more efficient?
	"""This class represents .gro files from GROMACS. 

	attributes:
	-----------
	- name		- the file name that the gro file came from
	- header	- first line of gro file; the title
	- natoms	- second line of the gro file (the number of atoms) represented as int
	- boxstring	- the last line of the gro file (the box information) represented as a string
	- atoms		- a GromacsAtoms object which contains all the atoms in the gro file

	methods:
	--------
        - appendatoms(atom)	- takes a GromacsAtoms object and appends it to the GromacsAtoms object in the structure
        - load(filename)	- load up data from a .gro file into self
        - getheader()		- returns the header (which should be hidden)
        - getboxstring()	- returns the box string (which should be hidden)
        - getboxvector()	- returns the box information as a tuple NOT IMPLEMENTED
	- remove(atom)		- removes atom from the list 
	- trim(a1,a2)		- get rid of atoms from a1 to a2
	- write(filename)	- writes the object to 'filename'

	"""

	def __init__(self, name="conf.gro", header="title", natoms=0, boxstring=" 10.0 10.0 10.0" ):
		self.name      = name 
		self.header    = header
		self.natoms    = natoms
		self.boxstring = boxstring
		#self.boxvector=[]
		self.atoms     = GromacsAtoms()	
	
	def appendatoms(self,newatoms):
		"""newatoms should be a GromacsAtoms object. This method appends this atomsobject to self.atoms. The atoms
		are renumbered to be consistent with the GromacsStructure object. The residue numbers are also modified,
		but several residues may be added to the GromacsStructure in this way."""

		try:
			max_residue=self.atoms[-1].resnum
		except:
			max_residue=0

		# have to modify the newatoms object in the following ways:
		# 	1. start the atom numbering from the original number of atoms in self
		#	2. start the residue numbering from the original number of residues in self
		oldresnum=0
		#print "appending %d atoms to %s" %(len(newatoms),self.name)
		for atom in newatoms:
			self.natoms=self.natoms+1
			atom.atomnum=self.natoms
			newresnum=atom.resnum	
			if newresnum != oldresnum:
				max_residue=max_residue+1
				oldresnum=newresnum
			atom.resnum=max_residue
			self.atoms.append(atom)	
		return

	def formatfilename(self,filename):
		if filename[-4:] != ".gro": filename=filename+".gro"
		return filename	

	def getheader(self): return self.header

        def getboxstring(self): return self.boxstring

        def getboxvector(self): pass # box in vector form not yet implemented

	def load(self,filename):
		"""Takes 'filename' and, if found, loads the data in the file to self. If not found, 
		raise a "GroFileNotFound" exception."""
		# if 'filename' doesn't end in '.gro', append it.
		filename=self.formatfilename(filename)
		if exists(filename):
			### build the gro file ###
			grofile=open(filename)
			structurefile=grofile.readlines()
			grofile.close()
			self.list2gro(structurefile)
		else: raise GroFileNotFound
		self.name      = filename
		return 

	def test(self,filename):
	        gro=GromacsStructure()
 	        gro.load(filename)
	        print gro.getheader()
	        print gro.natoms
	        print gro.atoms[0].resname
	        print gro.getboxstring()
	        print "Test complete."
		return

	def list2gro(self,list):
		"""The main GromacsStructure creation method. Takes a gro file which has been loaded
		as a list and parses it into a GromacsStructure."""

		self.header=list[0].strip("\n")
		try: self.natoms=int(list[1])
		except: raise BadAtomNumber
		self.boxstring=list[-1].strip("\n")
		atomlist=list[2:-1]
		for atom in atomlist:
			self.atoms.store(atom)
		self.name="structure-assigned-from-list.gro" 
		return

	def trim(self,atom1,atom2):
		"""A trimmer: removes atoms from atom1 to atom2."""
		self.atoms=self.atoms[0:atom1]+self.atoms[atom2+1:]
		self.natoms=len(self.atoms)
		return
	def write(self,filename):
		"""Writes the gro object to a file 'filename.gro'. Raise an EmptyGroObject
		exception if there are no atoms (self.natoms == 0)? or should it be okay to
		write an empty gro file?"""
		if self.natoms==0:
			raise EmptyGroObject
		else:
			filename=self.formatfilename(filename)
			f=open(filename,'w')
			f.write(self.header+"\n")
			f.write(str(self.natoms)+"\n")
			for atom in self.atoms:
				f.write(atom.printgroatom()+"\n")
			f.write(self.boxstring+"\n")
			f.close()
		return

### global module functions ###
"""
Global Functions
----------------
- author 	- print module authorship information
- gro2gss_atom_sorter - takes a GromacsStructure, a list of atom numbers, an optional list of charges,
 		 and an optional list of "special" atoms, and returns the atoms in a format suitable for
		 a calculation with Gaussian03 -- the numbered atoms will be included in the QM;
		 the "special" atoms will retain their names, whereas others will be truncated to the 
		 chemical symbol (CAREFUL! gamma hydrogen - HG - will be treated as mercury atom!). 
		 Atoms not in the list will be treated as point charges.
- test		- returns a gro object if run0.gro is in the current directory ('gro=GROMACS.test()')
- version	- print module version 
- xtc2gro	- convert an xtc file (using a tpr and trjconv) to gro; return a list of gro objects
"""

def author():
	print "GROMACS module brought to you by", AUTHOR
	print "(You're welcome.)"
	return

def sizeofxtc(xtc):
	"""Returns the number of frames in the xtc file."""
	gmxdump_cmd="gmxdump -f %s 1>dump 2>/dev/null" % xtc
	system(gmxdump_cmd)
	dump=open("dump")
	gmxdump=dump.readlines()
	dump.close()
	# start reading backwards for the word "frame"
	index=-1
	nframes=0 # default
	while (abs(index)<len(gmxdump)):
		if "frame" in gmxdump[index]:
			line=gmxdump[index].split()
			nframes=line[2].strip(":")
			nframes=int(nframes)
			break
		index-=1
	return nframes

def gro2gss_atom_sorter(gro, includelist, charges=[], specialatoms=[]):
        natoms=gro.natoms
        atomsection=""
        chargesection=""
        i=0
        while i<natoms:
                atom=gro.atoms[i]
                if i in includelist:
                        line="%s % 3.3f % 3.3f % 3.3f\n"
                        atomname=atom.atomname.strip()
                        if atomname in specialatoms:
                                line=line % (atomname, atom.x, atom.y, atom.z)
                                atomsection = line + atomsection
                        else:
                                line=line % (atomname[0], atom.x, atom.y, atom.z)
                                atomsection += line
                else:
			# assume that if the reference to charges[i] doesn't work,
			# then we don't want to include the atom in the calculation
			try:
                        	charge=charges[i]
	                        line="% 3.3f % 3.3f % 3.3f % 3.3f\n"
       	                	line=line % (atom.x, atom.y, atom.z, charge)
                        	chargesection += line
			except: pass
                i += 1

        # put it all together
        sorted_atom_list=atomsection+"\n"+chargesection+"\n"
        return sorted_atom_list

def test():
	version()
	author()
	print
	gro=GromacsStructure()
	try:
		gro.load("test.gro")
		print "loaded test.gro for testing"
		print "use 'gro=test()' or 'gro=GROMACS.test()' to play with a fun GromacsStructure object"
		print
	except GroFileNotFound:
		print "please place test.gro in the local directory"
	return gro

def version():
	print "GROMACS module: version %s, %s" % (VERSION, DATE)	
	return

def xtc2gro( xtcfile, tprfile, index, prefix = "echo 0", verbose = False, group = 0, prec = 3 ):
	# take a tpr and an xtc, split into a list of gro files, return a GromacsStructure object from that 
	# list with index 'index' 
	# the "group" variable if true says which group in the index file to get; default is zero
	# We keep track of the 't=' in the outputted gro file and only keep/consider the latest time in order
	# to deal with repeated frames.

	# This little function splits a list 'lst' into chunks of 'chunksize' elements. A list of the chunks
	# is returned. For instance, 
	#
	# 	>> a = splitlst( ( 1, 2, 3, 4) , 2 )
	#       >> print a
	# 	[(1, 2), (3, 4)]
	#
	def splitlist( lst, chunksize ):
	        min = 0
	        max = chunksize
	       	splitted = []
	        while max <= len( lst ):
	                splitted.append( lst[ min:max ] )
	                min = max
	                max = max + chunksize
        	return splitted

	grolist = []
	tmpgrofilename = xtcfile + ".TMP.gro"

	if verbose: print "looking for tpr file", tprfile
	if not exists( tprfile ): raise TprFileNotFound

	if verbose: print "looking for xtc file", xtcfile
	if not exists( xtcfile ): raise XtcFileNotFound

	if not exists( tmpgrofilename ): # no need to run trjconv multiple times
		try:
			args = ( prefix, tprfile, xtcfile, tmpgrofilename, prec )
			trjconv = "%s | trjconv -s %s -f %s -o %s -ndec %d" % args
			trjconv += ">& /dev/null"

			if verbose: print "trying command '%s'" % trjconv

			system( trjconv )

		except:
			raise TrjconvError
	
	biggro = open( tmpgrofilename )
	biggrolisti = biggro.readlines()
	biggro.close()

	unlink( tmpgrofilename )
	natoms = int( biggrolist[1] )	
	splittedgrolist = splitlist( biggrolist, natoms + 3 )
	tmpdir = "tmp/"

	# Grab frames which are contiguous in time, not repeats. Use only most recent
	# frames.
	golist = []
	for gro in splittedgrolist :
		title = gro[0]
		print title

	#print "my list is:"
	#for item in splittedgrolist: print item
	#print "looking for index",index
	gro=splittedgrolist[index]
	groobj=GromacsStructure()
	groobj.list2gro(gro)
	groobj.name=str(index)+".gro"	

	return groobj

### test functions ###
"""
test functions
--------------
test_xtc2gro(v)	- v=1 means verbose, v=0 means quiet. Runs xtc2gro on test.xtc, test.tpr
test_atomlist()	- are the atom lists for multiple objects independent?
"""
def test_xtc2gro(verbose):
	print "this test creates gro files from test.xtc using test.tpr"
	grolist=xtc2gro("test.xtc","test.tpr",verbose)
	return grolist

def test_atomlist():
	print "this test is for atom lists"
	print "runs the following:\n"
	print "gro1=GromacsStructure()"
	gro1=GromacsStructure()
	print "gro1.load('test.gro')"
	gro1.load("test.gro")
	print "gro2=GromacsStructure()"
	gro2=GromacsStructure()
	print "print gro2.atoms[0].printgroatom()"
	print gro2.atoms[0].printgroatom()
	print "print gro2.header"
	print gro2.header
	print "print gro1.header"
	print gro1.header
	return

### global module Exceptions ###
class BadAtomNumber(Exception): pass
class BadCoordinatesInGro(Exception): pass
class BadResidueNumber(Exception): pass
class EmptyGroObject(Exception): pass
class FileNotFound(Exception): pass
class GroFileNotFound(FileNotFound): pass
class GromacsProgramError(Exception): pass
class TprFileNotFound(FileNotFound): pass
class TrjconvError(GromacsProgramError): pass # perhaps the program isn't in the path?
class XtcFileNotFound(FileNotFound): pass
