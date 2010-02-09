from os.path import exists
# from numpy import *	# do I need this?
from os import system,unlink

from GromacsAtom import *


class GromacsAtomList( list ):
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

        def __copy__(self): return GromacsAtomList(self)

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

		try:
                    resnum=int(line[0:5])
		except:
                    raise BadResidueNumber

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
	
