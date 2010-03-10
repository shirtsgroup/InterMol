#############################
# structureTools.py
# 
# a set of tools for reading, writing, and manipulating *.gro and PDB structures
#
# -----------------------------------------------------
# MOD HISTORY
# 02/05/2010 - VAV Created structureTools.py
#              based on pdbtools.py
	# PDB reading and writing functions from packmoltools/resolvate, which provide
	# rudimentry input and output for pdbfiles.
#
# pdbtools written by John Chodera
# modified by Vincent Voelz
# Stanford University (c) 2007
#
# 07/10/2007 - Created pdbtools.py
#            - added ignh=True flag to writeAtomsToPDB() 
#            - enabled 4-letter gmx-style ffamber residue names (for termini e.g.) 
##############################


# GLOBAL imports

import sys, os, os.path
import commands, shutil, tempfile, math, string

from GromacsStructure import *
from PDBStructure import *


# global variables

# three letter amino acid code to one letter amino acid code
tlc2olc = {
"ALA" : "A" ,      # Non-terminal amino acids, and their ionized names too
"CYS" : "C" ,
"CYN" : "C" ,
"CYX" : "C" ,
"ASP" : "D" ,
"GLU" : "E" ,
"PHE" : "F" ,
"GLY" : "G" ,
"HIS" : "H" ,
"HIE" : "H" ,
"HID" : "H" ,
"HIP" : "H" ,
"ILE" : "I" ,
"NLE" : "J" ,
"LYS" : "K" ,
"LYN" : "K" ,
"LYP" : "K" ,
"LEU" : "L" ,
"MET" : "M" ,
"ASN" : "N" ,
"PRO" : "P" ,
"GLN" : "Q" ,
"ARG" : "R" ,
"SER" : "S" ,
"THR" : "T" ,
"VAL" : "V" ,
"TRP" : "W" ,
"TYR" : "Y" ,
"CALA" : "A" ,      # C-terminal amino acids, and their ionized names too
"CCYS" : "C" ,
"CCYN" : "C" ,
"CCYX" : "C" ,
"CASP" : "D" ,
"CGLU" : "E" ,
"CPHE" : "F" ,
"CGLY" : "G" ,
"CHIS" : "H" ,
"CHIE" : "H" ,
"CHID" : "H" ,
"CHIP" : "H" ,
"CILE" : "I" ,
"CNLE" : "J" ,
"CLYS" : "K" ,
"CLYN" : "K" ,
"CLYP" : "K" ,
"CLEU" : "L" ,
"CMET" : "M" ,
"CASN" : "N" ,
"CPRO" : "P" ,
"CGLN" : "Q" ,
"CARG" : "R" ,
"CSER" : "S" ,
"CTHR" : "T" ,
"CVAL" : "V" ,
"CTRP" : "W" ,
"CTYR" : "Y" ,
"NALA" : "A" ,      # N-terminal amino acids, and their ionized names too
"NCYS" : "C" ,
"NCYN" : "C" ,
"NCYX" : "C" ,
"NASP" : "D" ,
"NGLU" : "E" ,
"NPHE" : "F" ,
"NGLY" : "G" ,
"NHIS" : "H" ,
"NHIE" : "H" ,
"NHID" : "H" ,
"NHIP" : "H" ,
"NILE" : "I" ,
"NNLE" : "J" ,
"NLYS" : "K" ,
"NLYN" : "K" ,
"NLYP" : "K" ,
"NLEU" : "L" ,
"NMET" : "M" ,
"NASN" : "N" ,
"NPRO" : "P" ,
"NGLN" : "Q" ,
"NARG" : "R" ,
"NSER" : "S" ,
"NTHR" : "T" ,
"NVAL" : "V" ,
"NTRP" : "W" ,
"NTYR" : "Y" 
} 



def GromacsStructureFromGrofile(grofile):
    """Create a GromacsStructure objects from a grofile."""
    
    gstruct = GromacsStructure(name=grofile, header="title" )
    gstruct.load(grofile)
    return gstruct


def GromacsStructureFromPDB(pdbfile):
    """Take in a PDB file a GromacsStructure class."""
    p = PDBStructure(filename=pdbfile)
    gstruct = GromacsStructure(name=pdbfile, header="title" )
    gatoms = GromacsAtoms()
    for atom in p.atoms:
        gatoms.append( GromacsAtom(resnum=atom["resSeq"], resname=atom["resName"], atomname=atom["name"], atomnum=atom["serial"],
                                      x=(atom["x"]/10.), y=(atom["y"]/10.), z=(atom["z"]/10.), vx=0.0, vy=0.0, vz=0.0) )
    gstruct.appendatoms( gatoms )
    return gstruct 


def PDBFromGrofile(grofile, stripWaters=False, stripIons=False, stripDummy=False):
    """Take in a *.gro Grofile file and write it as a PDB file."""
    
    gstruct = GromacsStructure(name=grofile, header="title" )
    gstruct.load(grofile)
    
    p = PDBStructure()
    
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
        elif (atom["resName"].count('DMM')>0):
          if not stripDummy:
            p.atoms.append(atom)
        else:
            p.atoms.append(atom)

    return p


def writeGrofileFromPDB(pdbfile, grofile):
    """Take in a PDB file and write it as a *.gro file."""
    
    gstruct = GromacsStructureFromPDB(pdbfile)
    gstruct.write(grofile)
    
    
def writePDBFromGrofile(grofile, pdbfile, stripWaters=False, stripIons=False, stripDummy=False):
    """Take in a *.gro Grofile file and write it as a PDB file."""
    
    p = PDBFromGrofile(grofile, stripWaters=False, stripIons=False, stripDummy=False)
    p.writeAtomsToPDB(pdbfile, renumber = True)
    

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
    
    state1 = getChargedNames(g1.atoms)
    state2 = getChargedNames(g2.atoms)
            
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
    
    p1 = PDBStructure(filename=pdbfile1)
    p2 = PDBStructure(filename=pdbfile2)
        
    state1 = p1.getChargedNames()
    state2 = p2.getChargedNames()
            
    print pdbfile1,':',len(atoms1), 'atoms'
    print state1
    print 
    print pdbfile2,':',len(atoms2), 'atoms'
    print state2
    
    return (state1 == state2)


    
    
    
