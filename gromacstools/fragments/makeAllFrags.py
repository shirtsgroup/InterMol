#!/usr/bin/env python

import sys, os, glob
from zam import protein

# needed for tlc2olc
import mmtools.pdbtools as pdbtools

usage = """makeAllFrags.py OUTDIR [TYPE]

    Will build PDB (*.pdb) and gmx topology files (*.itp) for every amino acid fragment

    OUTDIR    The output directory

    TYPE      "sidechain" - Replaces all backbone atoms as extra C_beta hydrogen.  Default. 
           or "backbone"  - caps N- and C- termini with neutral hydrogen.
"""

# constants

olc2tlc = {'A':'ALA',
           'C':'CYS',
           'D':'ASP',
           'E':'GLU',
           'F':'PHE',
           'G':'GLY',
           'H':'HIS',
           'I':'ILE',
           'K':'LYS',
           'L':'LEU',
           'M':'MET',
           'N':'ASN',
           'P':'PRO',
           'R':'ARG',
           'S':'SER',
           'T':'THR',
           'V':'VAL',
           'W':'TRP',
           'Y':'TYR'}



# functions


def MakeFragment(seq, Outdir, Name, FragType='sidechain'):
    """Makes a fragment PDB and topoolgy files for a given sequence,
    with asdfghj
    
    REQUIRED
    
    Seq         one-letter amino acid seqiuence
    Outdir      Directory to write output files
    Name        basename for filenames
    
    PARAMETERS
    
    FragType    'sidechain' or 'backbone'
    """
    
    # make outdir if it doesn't exist yet
    if not os.path.exists(Outdir):
        os.mkdir(Outdir)
        
    p = protein.ProteinClass(Seq=seq)
    pdbfile = os.path.join(Outdir, Name+'.pdb')
    p.WritePdb(pdbfile)
    
    # read in the PDB file and collect the atom groups
    pdbtools.readAtomsFromPDB(pdbfile)
    print 
    
    
    
    
    
if __name__ == '__main__':
    
    
    
    if len(sys.argv) < 2:
        print usage
        sys.exit()
        
    outdir = sys.argv[1]
    fragtype = 'sidechain'
    if len(sys.argv) >= 3:
         fragtype = sys.argv[2]
         

    Verbose = True
    
    if Verbose:
        print 'outdir', outdir
        print 'fragtype', fragtype
        
    for seq in olc2tlc.keys():
        MakeFragment(seq, outdir, seq, FragType=fragtype)
        
        