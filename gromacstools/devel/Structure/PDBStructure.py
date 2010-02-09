from os.path import exists
# from numpy import *	# do I need this?
from os import system,unlink

from PDBAtom import *


# HISTORY

# VAV: 2/05/2010  Created; adapated from pdbtools PDB routines 


class PDBStructure:
    """A PDB object class."""

    def __init__(self, filename=None):
	
	self.filename  = filename
	self.atoms     = []
	
	if self.filename != None:
	    self.readAtomsFromPDB(self.filename)

	    
    def readAtomsFromPDB(self, pdbfilename):
	"""Read atom records from the PDB and return them in a list.
    
	REQUIRED ARGUMENTS
	  pdbfilename - the filename of the PDB file to import from
        
	The ATOM records are read, and the sequence for which there are atomic coordinates is stored.
    
	"""
    
	# Read the PDB file into memory.
	pdbfile = open(pdbfilename, 'r')
	lines = pdbfile.readlines()
	pdbfile.close()
    
    
	# Read atoms.
	self.atoms = []
	for line in lines:
	    if line[0:5] == "ATOM ":
		# Parse line into fields.
		atom = { }
		atom["serial"] = int(line[5:11])
		atom["name"] = line[12:16]
		atom["altLoc"] = line[16:17]
		atom["resName"] = line[17:21]
		atom["chainID"] = line[21:22]
		try:
		    atom["resSeq"] = int(line[22:26])
		except:
		    atom["resSeq"] = int(line[22:27])
		atom["iCode"] = line[26:27]
		atom["x"] = float(line[30:38])
		atom["y"] = float(line[38:46])
		atom["z"] = float(line[46:54])
    
		atom["occupancy"] = 1.0
		if (line[54:60].strip() != ''):
		  atom["occupancy"] = float(line[54:60])
		  
		atom["tempFactor"] = 0.0
		try:
		  if (line[60:66].strip() != ''):
		      atom["tempFactor"] = float(line[60:66])
		except:
		    pass
		atom["segID"] = line[72:76]
		atom["element"] = line[76:78]
		atom["charge"] = line[78:80]
		
		self.atoms.append(atom)
    
    
    def writeAtomsToPDB(self, pdbfilename, renumber = False, ignh = False):
	"""Write atom records to PDB file.
      
	REQUIRED ARGUMENTS
	  pdbfilename - the name of the PDB file to write

	OPTIONAL ARGUMENTS
	  if renumber is True, then the atom and residue numbers will be renumbered starting with 1
      
	RETURN VALUES
	  none
      
	EXAMPLE
	writeAtomsToPdb(pdbfilename)
	
	"""
      
	# Renumber if desired.
	if (renumber):
	  first_residue = self.atoms[0]["resSeq"]
	  serial = 1
	  resSeq = 0
	  last_resSeq = None
	  for atom in self.atoms:
	    atom["serial"] = serial
	    serial += 1
      
	    if(atom["resSeq"] != last_resSeq):
	      resSeq += 1
	      last_resSeq = atom["resSeq"]
	    atom["resSeq"] = resSeq      
	    
	# Read the PDB file into memory.
	pdbfile = open(pdbfilename, 'w')
	
	# Write atoms.
	for atom in self.atoms:
	  if (not ignh) or (not self.isAtomNameHydrogen(atom["name"] )):   # Ignore hydrogens if desired
	    pdbfile.write('ATOM %(serial)6d %(name)4s%(altLoc)c%(resName)4s%(chainID)c%(resSeq)4d%(iCode)c   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f%(element)2s%(charge)2s\n' % atom)
	pdbfile.close()
    
    def new_atom(self): 
	"""Return a new 'blank' atom dictionary, filled with dummy values."""
    
	atom = {}
	atom["x"] = 0.0
	atom["y"] = 0.0
	atom["z"] = 0.0
	atom["tempFactor"] = 1.0
	atom["resName"] = 'DUM '
	atom["serial"] = 1        
	atom["name"] = ' C  '
	atom["altLoc"] = ' '
	atom["chainID"] = 'A'
	atom["resSeq"] = 1
	atom["iCode"] = ' '
	atom["occupancy"] = 1.0
	atom["segID"] = ''
	atom["element"] = ''
	atom["charge"] = ''
    
	return atom
    
    
    def isAtomNameHydrogen(self, atomName):
	"""Returns True if the atomName denotes a non-water hydrogen."""
	
	if atomName[1] == 'H':
	    if atomName.count('W') == 0:    # water hydrogens don't count
		return True
	return False
    
    
    def getChargedNames(self):
	"""Find out which residues are charged based on their AMBER-style residue names.
	
	REQUIRED ARGUMENTS
	atoms - a list of atom dictionaries (if filetype='pdb') - see readAtomsFromPdb
	   OR
	      - a GromacsAtoms object  (if filetype='gro') - see mmtools.gromacstools.GromacsData
	   
	
	RETURN VALUES
	chargedRes, a dictionary { resSeq: ('+' or '-' or '0', "resName") } depending on whether the titratible side chain residue charged or netural
	
	NOTES
	The criteria here come from David Mobley renaming fixes in mmtools/mccetools."""
    
	chargedRes = {}    # { resSeq: '+' or '-'}
	
	for atom in self.atoms:
	    
	    # tritratible residues             
	    if atom["resName"].strip() == 'HIP':
		chargedRes[ atom["resSeq"] ] = ('+', atom["resName"])
	    elif atom["resName"].strip() == 'HIS':
		chargedRes[ atom["resSeq"] ] = ('0', atom["resName"])
	    elif atom["resName"].strip() == 'HID':
		chargedRes[ atom["resSeq"] ] = ('0', atom["resName"])
	    elif atom["resName"].strip() == 'HIE':
		chargedRes[ atom["resSeq"] ] = ('0', atom["resName"])
	    elif atom["resName"].strip() == 'LYN':
		chargedRes[ atom["resSeq"] ] = ('0', atom["resName"])
	    elif atom["resName"].strip() == 'LYP':
		chargedRes[ atom["resSeq"] ] = ('+', atom["resName"])
	    elif atom["resName"].strip() == 'CYN':
		chargedRes[ atom["resSeq"] ] = ('0', atom["resName"])
	    elif atom["resName"].strip() == 'CYM':
		chargedRes[ atom["resSeq"] ] = ('-', atom["resName"])
	    elif atom["resName"].strip() == 'ASH':
		chargedRes[ atom["resSeq"] ] = ('0', atom["resName"])
	    elif atom["resName"].strip() == 'GLH':
		chargedRes[ atom["resSeq"] ] = ('0', atom["resName"])
		 
    
	    # remember all 4-letter residue names, which are the termini
	    if len(atom["resName"].strip()) == 4:
		chargedRes[ atom["resSeq"] ] = ('0', atom["resName"])
	    
	    # Additonally, remember the charged names 
	    #
	    # N-termini names     
	    if atom["resName"].strip() == 'NHIP':
		chargedRes[ atom["resSeq"] ] = ('+', atom["resName"])
	    # gmx ffamber does not have CLYN, only CLYP 
	    elif atom["resName"].strip() == 'NLYP':
		chargedRes[ atom["resSeq"] ] = ('+', atom["resName"])
	    elif atom["resName"].strip() == 'NCYM':
		chargedRes[ atom["resSeq"] ] = ('-', atom["resName"])
	    # C-termini names     
	    elif atom["resName"].strip() == 'CHIP':
		chargedRes[ atom["resSeq"] ] = ('+', atom["resName"])
	    # gmx ffamber does not have CLYN, only CLYP 
	    elif atom["resName"].strip() == 'CLYP':
		chargedRes[ atom["resSeq"] ] = ('+', atom["resName"])
	    elif atom["resName"].strip() == 'CCYM':
		chargedRes[ atom["resSeq"] ] = ('-', atom["resName"]) 
		    
	return chargedRes        
    
	
    def setChargedNames(chargedRes):
	"""Change the names of the residues based on the charged state indicated in the chargedRes dict  
    
	REQUIRED ARGUMENTS
	chargedRes - a dictionary { resSeq: '+' or '-' or '0' }
	
	RETURN VALUES
	atoms      - a list of atom dictionaries, renamed to residues in chargedRes dictionary
	
	NOTE: This only works if we ignore the hydrogens when writing the PDB."""
	
	
	for i in range(len(self.atoms)):
	    if chargedRes.has_key(self.atoms[i]["resSeq"] ):
		self.atoms[i]["resName"] = chargedRes[ self.atoms[i]["resSeq"] ][1]
	

