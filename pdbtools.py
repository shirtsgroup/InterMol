# pdbtools.py
#
# PDB reading and writing functions from packmoltools/resolvate, which provide
# rudimentry input and output for pdbfiles.
#
# written by John Chodera, modified by Vincennt Voelz, Stanford University (c) 2007
# -----------------------------------------------------
#
# -----------------------------------------------------
# TO DO
# - This should be revamped with the mmLib module
# -----------------------------------------------------
# MOD HISTORY
# 07/10/2007 - Created pdbtools.py
#            - added ignh=True flag to writeAtomsToPDB() 
#            - enabled 4-letter gmx-style ffamber residue names (for termini e.g.) 


# GLOBAL imports

import sys, os, os.path
import commands, shutil, tempfile, math

from mmtools.gromacstools.GromacsData import *


def readAtomsFromPDB(pdbfilename):
    """Read atom records from the PDB and return them in a list.

    present_sequence = getPresentSequence(pdbfilename, chain=' ')
    
    REQUIRED ARGUMENTS
      pdbfilename - the filename of the PDB file to import from

    OPTIONAL ARGUMENTS
      chain - the one-character chain ID of the chain to import (default ' ')

    RETURN VALUES
      sequence (dictionary) - sequence[residue_id] is the one-letter code corresponding to residue index residue_id

    The ATOM records are read, and the sequence for which there are atomic coordinates is stored.

    """

    # Read the PDB file into memory.
    pdbfile = open(pdbfilename, 'r')
    lines = pdbfile.readlines()
    pdbfile.close()

    # Read atoms.
    atoms = [ ]
    for line in lines:
        if line[0:6] == "ATOM  ":
            # Parse line into fields.
            atom = { }
            atom["serial"] = int(line[6:11])
            atom["name"] = line[12:16]
            atom["altLoc"] = line[16:17]
            atom["resName"] = line[17:21]
            atom["chainID"] = line[21:22]
            atom["resSeq"] = int(line[22:26])
            atom["iCode"] = line[26:27]
            atom["x"] = float(line[30:38])
            atom["y"] = float(line[38:46])
            atom["z"] = float(line[46:54])

            atom["occupancy"] = 1.0
            if (line[54:60].strip() != ''):
              atom["occupancy"] = float(line[54:60])
              
            atom["tempFactor"] = 0.0
            if (line[60:66].strip() != ''):
              atom["tempFactor"] = float(line[60:66])
            
            atom["segID"] = line[72:76]
            atom["element"] = line[76:78]
            atom["charge"] = line[78:80]
            
            atoms.append(atom)
            
    # Return list of atoms.
    return atoms

def writeAtomsToPDB(pdbfilename, atoms, renumber = False, ignh = False):
    """Write atom records to PDB file.
  
    REQUIRED ARGUMENTS
      pdbfilename - the name of the PDB file to write
      atoms - a list of atom dictionaries -- see readAtomsFromPdb
  
    OPTIONAL ARGUMENTS
      if renumber is True, then the atom and residue numbers will be renumbered starting with 1
  
    RETURN VALUES
      none
  
    EXAMPLE
    writeAtomsToPdb(pdbfilename, atoms)
    
    """
  
    # Renumber if desired.
    if (renumber):
      first_residue = atoms[0]["resSeq"]
      serial = 1
      resSeq = 0
      last_resSeq = None
      for atom in atoms:
        atom["serial"] = serial
        serial += 1
  
        if(atom["resSeq"] != last_resSeq):
          resSeq += 1
          last_resSeq = atom["resSeq"]
        atom["resSeq"] = resSeq      
        
    # Read the PDB file into memory.
    pdbfile = open(pdbfilename, 'w')
    
    # Write atoms.
    for atom in atoms:
      if (not ignh) or (not isAtomNameHydrogen(atom["name"] )):   # Ignore hydrogens if desired
        pdbfile.write('ATOM  %(serial)5d %(name)4s%(altLoc)c%(resName)4s%(chainID)c%(resSeq)4d%(iCode)c   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f%(element)2s%(charge)2s\n' % atom)
    pdbfile.close()

def isAtomNameHydrogen( atomName ):
    """Returns True if the atomName denotes a non-water hydrogen."""
    
    if atomName[1] == 'H':
        if atomName.count('W') == 0:    # water hydrogens don't count
            return True
    return False


def getChargedNames(atoms):
    """Find out which residues are charged based on their AMBER-style residue names.
    
    REQUIRED ARGUMENTS
    atoms - a list of atom dictionaries - see readAtomsFromPdb
    
    RETURN VALUES
    chargedRes, a dictionary { resSeq: ('+' or '-' or '0', "resName") } depending on whether the titratible side chain residue charged or netural
    
    NOTES
    The criteria here come from David Mobley renaming fixes in mmtools/mccetools."""

    chargedRes = {}    # { resSeq: '+' or '-'}
    for atom in atoms:
        
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

    
def setChargedNames(atoms, chargedRes):
    """Change the names of the residues based on the charged state indicated in the chargedRes dict  

    REQUIRED ARGUMENTS
    atoms      - a list of atom dictionaries -- see readAtomsFromPdb
    chargedRes - a dictionary { resSeq: '+' or '-' or '0' }
    
    RETURN VALUES
    atoms      - a list of atom dictionaries, renamed to residues in chargedRes dictionary
    
    NOTE: This only works if we ignore the hydrogens when writing the PDB."""
    
    
    for atom in atoms:
        if chargedRes.has_key(atom["resSeq"] ):
            atom["resName"] = chargedRes[ atom["resSeq"] ][1]
    
    return atoms
    


def writeGrofileFromPDB(pdbfile, grofile):
    """Take in a PDB file and write it as a *.gro file."""
    
    atoms = readAtomsFromPDB(pdbfile)
    gstruct = GromacsStructure(name=grofile, header="title" )
    for atom in atoms:
        gatomlist = []
        gatomlist.append( GromacsAtom(self, resnum=atom["resSeq"], resname=atom["resName"], atomname=atom["name"], atomnum=atom["serial"],
                                      x=(atom["x"]/10.), y=(atom["y"]/10.), z=(atom["z"]/10.), vx=0.0, vy=0.0, vz=0.0) )
        gatoms = GromacsAtoms( gatomlist )
        gstruct.appendatoms( gatoms )
    gstruct.write(grofile)
    
    
def writePDBFromGrofile(grofile, pdbfile):
    """Take in a *.gro Grofile file and write it as a PDB file."""
    
    gstruct = GromacsStructure(name=grofile, header="title" )
    gstruct.load(grofile)
    atoms = []
    for gatom in gstruct.atoms:

        atom = { }
        atom["serial"] = gatom.atomnum
        atom["name"] = '%-4s'%(gatom.atomname.strip())
        atom["altLoc"] = ' '
        atom["resName"] = '%-4s'%(gatom.resname.strip())
        atom["chainID"] = ' '        
        atom["resSeq"] = gatom.resnum
        atom["iCode"] = ' '
        atom["x"] = gatom.x*10.    # *.gro files are in nm
        atom["y"] = gatom.y*10.
        atom["z"] = gatom.z*10.
        atom["occupancy"] = 1.0
        atom["tempFactor"] = 0.0
        atom["segID"] = '    '
        atom["element"] = '  '
        atom["charge"] = '   '

        atoms.append(atom)

    writeAtomsToPDB(pdbfile, atoms, renumber = True)
        
        
    
    

    