# pdbtools.py
#
# PDB reading and writing functions from packmoltools/resolvate, which provide
# rudimentry input and output for pdbfiles.
#
# written by John Chodera
# modified by Vincent Voelz
# Stanford University (c) 2007
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
import commands, shutil, tempfile, math, string

from mmtools.gromacstools.GromacsData import *


def readAtomsFromPDB(pdbfilename):
    """Read atom records from the PDB and return them in a list.

    present_sequence = getPresentSequence(pdbfilename, chain=' ')
    contents of protein.seqfile
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


def getChargedNames(atoms, type='pdb'):
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
    
    if type == 'pdb':
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
            
            
    elif type == 'gro':
      for atom in atoms:
        
        # tritratible residues             
        if atom.resname.strip() == 'HIP':
            chargedRes[ atom.resnum ] = ('+', atom.resname)
        elif atom.resname.strip() == 'HIS':
            chargedRes[ atom.resnum ] = ('0', atom.resname)
        elif atom.resname.strip() == 'HID':
            chargedRes[ atom.resnum ] = ('0', atom.resname)
        elif atom.resname.strip() == 'HIE':
            chargedRes[ atom.resnum ] = ('0', atom.resname)
        elif atom.resname.strip() == 'LYN':
            chargedRes[ atom.resnum ] = ('0', atom.resname)
        elif atom.resname.strip() == 'LYP':
            chargedRes[ atom.resnum ] = ('+', atom.resname)
        elif atom.resname.strip() == 'CYN':
            chargedRes[ atom.resnum ] = ('0', atom.resname)
        elif atom.resname.strip() == 'CYM':
            chargedRes[ atom.resnum ] = ('-', atom.resname)
        elif atom.resname.strip() == 'ASH':
            chargedRes[ atom.resnum ] = ('0', atom.resname)
        elif atom.resname.strip() == 'GLH':
            chargedRes[ atom.resnum ] = ('0', atom.resname)
             

        # remember all 4-letter residue names, which are the termini
        if len(atom.resname.strip()) == 4:
            chargedRes[ atom.resnum ] = ('0', atom.resname)
        
        # Additonally, remember the charged names 
        #
        # N-termini names     
        if atom.resname.strip() == 'NHIP':
            chargedRes[ atom.resnum ] = ('+', atom.resname)
        # gmx ffamber does not have CLYN, only CLYP 
        elif atom.resname.strip() == 'NLYP':
            chargedRes[ atom.resnum ] = ('+', atom.resname)
        elif atom.resname.strip() == 'NCYM':
            chargedRes[ atom.resnum ] = ('-', atom.resname)
        # C-termini names     
        elif atom.resname.strip() == 'CHIP':
            chargedRes[ atom.resnum ] = ('+', atom.resname)
        # gmx ffamber does not have CLYN, only CLYP 
        elif atom.resname.strip() == 'CLYP':
            chargedRes[ atom.resnum ] = ('+', atom.resname)
        elif atom.resname.strip() == 'CCYM':
            chargedRes[ atom.resnum ] = ('-', atom.resname) 
            
            
    else:
        print 'Error in getChargedNames():  must have filetype=\'pdb\' or \'gro\''

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
    
    gstruct = GromacsStructureFromPDB(pdbfile)
    gstruct.write(grofile)
   
def GromacsStructureFromPDB(pdbfile):
    """Take in a PDB file a GromacsStructure class."""
    atoms = readAtomsFromPDB(pdbfile)
    gstruct = GromacsStructure(name=pdbfile, header="title" )
    gatoms = GromacsAtoms()
    for atom in atoms:
        gatoms.append( GromacsAtom(resnum=atom["resSeq"], resname=atom["resName"], atomname=atom["name"], atomnum=atom["serial"],
                                      x=(atom["x"]/10.), y=(atom["y"]/10.), z=(atom["z"]/10.), vx=0.0, vy=0.0, vz=0.0) )
    gstruct.appendatoms( gatoms )
    return gstruct 

def GromacsStructureFromGrofile(grofile):

    gstruct = GromacsStructure(name=grofile, header="title" )
    gstruct.load(grofile)
    return gstruct

def oneLetterSeqFromGrofile(grofile):
    """Given a *.gro file, read in the residue sequence, and return a one-letter sequence code
    
    ***NOTE***:  This will *NOT* extract the correct one-letter code for threading D-amino acids with Modeller!!
    This is because the residue names and topologies are the same for D- and L-aminoi acids. 
    
    The naming conversion is defined in the dictionary  GromacsData.tlc2olc
    """
    
    # get the sequence from the grofile
    g = GromacsStructureFromGrofile(grofile)
    g.test(grofile)  # testing   
    return oneLetterSeqFromGromacsStructure(g)

def oneLetterSeqFromGromacsStructure(g):
    """Given a GromacsStructure, read in the residue sequence, and return a one-letter sequence code
        ***NOTE***:  This will *NOT* extract the correct one-letter code for threading D-amino acids with Modeller!!
    This is because the residue names and topologies are the same for D- and L-aminoi acids. 
    
    The naming conversion is defined in the dictionary  GromacsData.tlc2olc    """

    seq = g.atoms.sequence()
    print 'BEFORE', seq[1:100],'...'
    # the sequence may contain termini caps, ions and SOL residues, so filter out only the protein residues      
    i =0
    while i < len(seq):
        if tlc2olc.keys().count( seq[i].strip() ) > 0:
            i=i+1
        else:
            seq.pop(i)
    print 'AFTER', seq
    # convert the three-letter seq to one-letter:
    oneletterseq = []
    for s in seq:
        oneletterseq.append( tlc2olc[ s.strip() ] )
    return string.join(oneletterseq,'')
    
    
def mapSequenceOntoGromacsStructure(newseq, gstruct):
    """Map new GromacsStructure.atoms.sequence resnames onto a similar grofile, overriding any of the previous protonation states.
    ***NOTE***:  This will leave/remove spurious hydrogens that don't match  the prescribed protonation state!
    A 'pdb2gmx -ignh' should take care of this.
    ALSO:  It is assumed that the residue numbering the *.gro file starts at 1.
    
    RETURN VALUE
    gstruct                   The update GromacsStructure object"""
    
    oldseq = gstruct.atoms.sequence()
    
    oldindices = []
    newresnames = []

   
    oldresnum_start = gstruct.atoms[0].resnum
    newresnum_start = 1
    
    print 'oldresnum_start',oldresnum_start
    
    # find residue numbering of *just* the protein residues
    for i in range(0,len(oldseq)):
        if tlc2olc.keys().count( oldseq[i].strip() ) > 0:
            oldindices.append( i+oldresnum_start )
    for i in range(0,len(newseq)):
        if tlc2olc.keys().count( newseq[i].strip() ) > 0:
            newresnames.append( newseq[i] )
            
    # make sure we get the same number of residues
    if len(oldindices) != len(newresnames):
        print 'In mapSequenceOntoGromacsStructure(): len(oldindices) != len(newresnames)'
        raise Exception
        
    # go through the indices and change all the residue names of all
    for atom in gstruct.atoms:
        if oldindices.count(atom.resnum) > 0:
            atom.resname = newresnames[ oldindices.index(atom.resnum) ]

    return gstruct


def countSolventResiduesFromGrofile(grofile):
    """Given a *.gro file, returns the number of Cl ions, Na ions, and the number solvent residues"""
      
    g = GromacsStructureFromGrofile(grofile)
    return countSolventResiduesFromGromacsStructure(g)
      
def countSolventResiduesFromGromacsStructure(gstruct):
    """Given a *.gro file, returns the number of Cl ions, Na ions, and the number solvent residues"""
      
    seq = gstruct.atoms.sequence()
    nClres = 0
    nNares = 0
    nsolres = 0
    for res in seq:
	if res.count('Cl') > 0:
	    nClres += 1
        if res.count('Na') > 0:
	    nNares += 1	 	
	if res.count('SOL') > 0:
	    nsolres += 1
    return [ nClres, nNares, nsolres ]
      

def writePDBFromGrofile(grofile, pdbfile, stripWaters=False, stripIons=False):
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

        if (atom["resName"].count('HOH')>0) or (atom["resName"].count('SOL')>0):
          if not stripWaters:
            atoms.append(atom)
        elif (atom["resName"].count('Cl')>0) or (atom["resName"].count('Na')>0):
          if not stripIons:
            atoms.append(atom)
        else:
            atoms.append(atom)

    writeAtomsToPDB(pdbfile, atoms, renumber = True)
        
        
def rebuildPatchedTermini(outPdbFile):
    """Takes in a PDB file with 'patched' termini from Modeller, and rewrites it with the ffAMBER-style
    ACE and NH2 and residues."""
    
    # The Modeller PDBs look like this:
    """
     EXPDTA    THEORETICAL MODEL, MODELLER 9v1 2007/08/20 13:48:30 
     REMARK   6 MODELLER OBJECTIVE FUNCTION:      4986.1367
     REMARK   6 MODELLER BEST TEMPLATE % SEQ ID: 100.000
     ATOM      1  N   LEU     1      -7.461   4.812   8.221  1.00387.74       1SG   2
     ATOM      2  CA  LEU     1      -7.033   6.174   8.059  1.00387.74       1SG   3
     ATOM      3  CB  LEU     1      -7.447   7.022   9.265  1.00387.74       1SG   4
     ATOM      4  CG  LEU     1      -8.557   8.090   8.962  1.00387.74       1SG   5
     ATOM      5  CD1 LEU     1      -8.754   8.471   7.469  1.00387.74       1SG   6
     ATOM      6  CD2 LEU     1      -9.863   7.852   9.710  1.00387.74       1SG   7
     ATOM      7  C   LEU     1      -5.625   6.415   7.559  1.00387.74       1SG   8
     ATOM      8  O   LEU     1      -4.618   5.810   7.920  1.00387.74       1SG   9
cap* ATOM      9  CAY LEU     1      -9.681   4.051   7.379  1.00387.74       1SG  10
cap* ATOM     10  CY  LEU     1      -8.330   4.670   7.236  1.00387.74       1SG  11
cap* ATOM     11  OY  LEU     1      -8.069   5.220   6.173  1.00387.74       1SG  12
     ATOM     12  N   SER     2      -5.488   7.387   6.667  1.00239.19       1SG  13
     ATOM     13  CA  SER     2      -5.907   8.699   7.035  1.00239.19       1SG  14
     ....
     ATOM    279  C   LEU    34     -14.712  13.829   3.507  1.00350.49       1SG 280
     ATOM    280  O   LEU    34     -13.790  13.944   2.701  1.00350.49       1SG 281
     ATOM    281  N   PHE    35     -14.937  14.531   4.654  1.00213.39       1SG 282
     ATOM    282  CA  PHE    35     -14.196  14.588   5.890  1.00213.39       1SG 283
     ATOM    283  CB  PHE    35     -14.683  14.059   7.205  0.50213.39       1SG 284
     ATOM    284  CG  PHE    35     -15.624  14.612   8.170  0.50213.39       1SG 285
     ATOM    285  CD1 PHE    35     -15.385  15.592   9.112  0.50213.39       1SG 286
     ATOM    286  CD2 PHE    35     -16.780  13.907   8.118  0.50213.39       1SG 287
     ATOM    287  CE1 PHE    35     -16.409  15.856   9.992  0.50213.39       1SG 288
     ATOM    288  CE2 PHE    35     -17.793  14.160   8.976  0.50213.39       1SG 289
     ATOM    289  CZ  PHE    35     -17.589  15.139   9.915  0.50213.39       1SG 290
     ATOM    290  C   PHE    35     -13.342  13.381   5.900  1.00213.39       1SG 291
     ATOM    291  O   PHE    35     -12.488  13.238   5.039  1.00213.39       1SG 292
*cap ATOM    292  NT  PHE    35     -13.553  12.516   6.942  1.00213.39       1SG 293
     END

"""

            
    # Copy the original PDB to a temp file, to preserve a record off it
    tmpPdbFile = outPdbFile + '.before'
    fin = open(outPdbFile,'r')
    fout = open(tmpPdbFile, 'w')
    fout.write( fin.read() )
    fin.close()
    fout.close()
    
    # Read in the atoms from the PDB
    atoms = readAtomsFromPDB( outPdbFile )
    
    # Find the atoms in the termini caps
    nter_atoms = []
    cter_atoms = []
    i=0
    while i < len(atoms):
        if atoms[i]['name'].strip() == 'CY':
            atoms[i]['name'] = ' C  '
            nter_atoms.insert(0, atoms.pop(i))
        elif atoms[i]['name'].strip() == 'CAY':
            atoms[i]['name'] = 'CH3 '
            nter_atoms.insert(0, atoms.pop(i) )
        elif atoms[i]['name'].strip() == 'OY':
            atoms[i]['name'] = ' O  '
            nter_atoms.insert(0, atoms.pop(i) )
        elif atoms[i]['name'].strip() == 'NT':
            atoms[i]['name'] = ' N  '
            cter_atoms.insert(0, atoms.pop(i) )                
        else:
            i += 1 
            
    # print 'nter_atoms', nter_atoms, 'cter_atoms', cter_atoms
    
    # Rename, shift numbering, and tack back onto the atoms list
    for atom in nter_atoms:
        atom['resName'] = 'ACE '
        atom['resSeq'] = atom['resSeq'] - 1
        atoms.insert(0, atom)
    for atom in cter_atoms:
        atom['resName'] = 'NH2 '
        atom['resSeq'] = atom['resSeq'] + 1
        atoms.append( atom )
    
    # write he atoms to PDB -- this will take care of renumbering too
    writeAtomsToPDB(outPdbFile, atoms, renumber = True)


def isGrofileProtonationSame( grofile1, grofile2 ):
    """Checks all the charged residues in a grofile to see if the protonation state is the same.
    NOTE: The protontation state check does NOT include the charge state of the termini.
    Example: {1: ('0',NLEU)} means that the leucine sidechain of residue 1 is not charged ('0'), but says nothing about the termini"""
    
    g1 = GromacsStructureFromGrofile(grofile1)
    g2 = GromacsStructureFromGrofile(grofile2)
    
    state1 = getChargedNames(g1.atoms, type='gro')
    state2 = getChargedNames(g2.atoms, type='gro')
            
    print grofile1,':',len(g1.atoms), 'atoms'
    print state1
    print 
    print grofile2,':',len(g2.atoms), 'atoms'
    print state2
    
    return (state1 == state2)
    
    
def isPDBProtonationSame( pdbfile1, pdbfile2 ):
    """Checks all the charged residues in a pdbfile to see if the protonation state is the same.
    NOTE: The protontation state check does NOT include the charge state of the termini.
    Example: {1: ('0',NLEU)} means that the leucine sidechain of residue 1 is not charged ('0'), but says nothing about the termini"""
    
    atoms1 = readAtomsFromPDB(pdbfile1)
    atoms2 = readAtomsFromPDB(pdbfile2)
    
    state1 = getChargedNames(atoms1, type='pdb')
    state2 = getChargedNames(atoms2, type='pdb')
            
    print pdbfile1,':',len(atoms1), 'atoms'
    print state1
    print 
    print pdbfile2,':',len(atoms2), 'atoms'
    print state2
    
    return (state1 == state2)

    
    
    
    
    
