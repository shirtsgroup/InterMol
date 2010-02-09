from os.path import exists
# from numpy import *	# do I need this?
from os import system,unlink


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

