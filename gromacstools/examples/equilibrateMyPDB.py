# EquilibrateMyPDB.py 
#
# This example script takes the pdbfile "test.pdb" as input, and outputs
# everything you need to perform a GROMACS molecular dynamics equilibration:
#
#     GROMACS *.tpr
#     GROMACS *.gro
#     ./equilibrate        script that calls mdrun
#     prepare.log          logfile
#
# In this case, the PDB used for testing is a the first 13 residues of 3gb1, which has been
# pre-"puurified" by the mccetools scriptss, with the correct atom and residue names for the AMBER forcefields ports
# The PDB file is the NMR structure
# of the B1 domain from the staphloccocus IgG-binding protein G.  The protonation state for
# this test PDB has a net charge of +1, but the correct number of counterions consistent with
# the salt concentration are automatically added to form a neutral system
#
# Vincent Voelz
# May 21, 2007
#

# import mmtools.gromacstools.system as system     # the old gromacstools way
from mmtools.gromacstools.System import *

import os
import os.path
import sys
#=============================================================================================

# DEBUG flag: if true then save outputs
DEBUG = True


#----------------------------------
# CLASSES
#----------------------------------

class PipelineProtein:

    def __init__(self, pdb=None, seq=None, salt='NaCl', saltconc=0.100, pH=7.0):
        """Initialize a default data in protein object."""
        
        self.pdbfile  = pdb
        self.seqfile  = seq
        self.salt     = salt         # Molar concentration of salt (ionic strength)
        self.saltconc = saltconc
        self.pH       = pH

    
#----------------------------------
# FUNCTIONS
#----------------------------------
        
 
def shoveItThrough(protein, forcefield):
    
    thisdir = os.curdir
    
    # determine base name
    print "PDB File: " + protein.pdbfile
    print "Sequence File: " + protein.seqfile
    extInd = protein.seqfile.rindex(".")
    baseName = protein.seqfile[0:extInd]
        
    # run gromacs setup
    if (1):
        outname = baseName + "_final"
        # g = system.GromacsSystem(mcceOut, useff=forcefield)    # the old gromacstools way
        g = System(protein.pdbfile, useff=forcefield)
        g.setup.setSaltConditions(protein.salt, protein.saltconc)
        g.setup.set_boxSoluteDistance(0.9)   # <--- ***Greg!!!*** <--- periodic box margin distance, in nanometers 
        outdir = os.path.join(thisdir, baseName)
        print 'Writing equilibration directory to',outdir,'...'
        if os.path.exists(outdir) == False:
            os.mkdir(outdir)
        g.prepare(outname=outname, outdir=outdir, verbose=True, cleanup=False, debug=DEBUG, protocol='racecar2', checkForFatalErrors=True)
    

#----------------------------------
# MAIN
#----------------------------------

if __name__ == '__main__':

    pdbfile = 'test.pdb'            # PDB file containing protein to solvate
    sequence = 'test.seq'           # file or string containing the one-letter sequence
    forcefield = 'ffamber99p'       # gromacs forcefield name to use for grompp
    salt = 'NaCl'	            # salt pair for counterions - supported types are in the <forcefield>.rtp file
    saltconc = 0.150	
    pH = 7.0
        
    shoveItThrough( PipelineProtein( pdb=pdbfile, seq=sequence, salt=salt, saltconc=saltconc, pH=pH), forcefield=forcefield  )
        

