#=============================================================================================
# shoveit.py
#
# Shoves a protein through the entire MODELLER -> MCCE -> gromacs pipeline
#=============================================================================================
# AUTHORS
#
# Written by Vincent Voelz, Stanford University, 2007-07-06
#=============================================================================================
# REQUIREMENTS
#=============================================================================================
# VERSION CONTROL INFORMATION
__version__ = "$Revision: $"                                                                     

#=============================================================================================
# IMPORTS
#=============================================================================================

import os, os.path, sys

# GLOBAL IMPORTS
import mmtools.mccetools.mcce as mcce

# import mmtools.gromacstools.system as system     # the old gromacstools way
from mmtools.gromacstools.System import *

from optparse import OptionParser    # For parsing of command line arguments

# pipeline tools
from pipeline import * 
from thread_model import *


#=============================================================================================

# DEBUG flag: if true then save outputs
DEBUG = True


#---------------------------------
# FUNCTIONS
#---------------------------------


def shoveit(protein, outdir, forcefield='ffamber99p'):
    """Shoves a pipelineProtein object through the entire
       MODELLER -> MCCE --> gromacs pipeline."""


    print 'before:'
    protein.print_info()

    protein.setup(outdir)

    print 'after:'
    protein.print_info()

    # Build a PDB model from the pdbTemplate using MODELLER   
    thread_model(protein.pdbfile, protein.seqfile, protein.modelPDBout) 

    # Find the best protonation state using MCCE  
    # protonate(protein.pdbfile, protein.seqfile, protein.modelPDBout) 

    # Prepare a Gromacs simulation according to the pipelineProtein specs 
    # prepare_gmx(protein.pdbfile, protein.seqfile, protein.modelPDBout) 


    # run mcce
    if (0):
        mcceOut = os.path.abspath(mcceOut)
        prmFile = '../../mccetools/prmfiles/run.prm.quick'
        prmFile = os.path.abspath(prmFile)
        mcce.protonatePDB(modelPDBOut, mcceOut, protein.pH, os.environ['MCCE_LOCATION'], cleanup=True, prmfile=prmFile)

    # run gromacs setup
    if (0):
        gromacsOut = baseName + "_final.pdb"
        # g = system.GromacsSystem(mcceOut, useff=forcefield)    # the old gromacstools way
        g = System(mcceOut, useff=forcefield)
        g.setup.setSaltConditions(protein.salt, protein.saltconc)
        if protein.boxProtocol == 'small':
            g.setup.set_boxType = 'octahedron'
            g.setup.set_boxSoluteDistance(1.5)   # periodic box margin distance, in nanometers
        elif protein.boxProtocol == 'big':
            g.setup.set_boxType = 'octahedron'
            g.setup.setUseAbsBoxSize(True)
            g.setup.setAbsBoxSize('7.0')   # periodic box absolute size, in nanometers (string)

        thisOutDir = os.path.join(thisdir, baseName)
        print 'Writing equilibration directory to',thisOutDir,'...'
        if os.path.exists(thisOutDir) == False:
            os.mkdir(thisOutDir)
        g.prepare(outname=gromacsOut, outdir=thisOutDir, verbose=True, cleanup=False, debug=DEBUG, protocol='racecar2', checkForFatalErrors=True)

    if(0):
        # cleanup    if not DEBUG:
        os.remove(modelPDBOut)
        os.remove(mcceOut)



#----------------------------------
# MAIN
#----------------------------------


if __name__ == '__main__':

    # Create command-line argument options.
    usage_string = """%prog --pdb PDBFILE --seq SEQUENCE --outdir OUTDIR--salt SALT_TYPE --saltconc SALT_CONC --pH PHVALUE --boxprotocol PROTOCOL 
Shoves a protein all the way through the MODELLER -> MCCE -> gromacs pipeline.
For detailed help:  shoveit.py -h
"""

    version_string = "%prog %__version__"

    parser = OptionParser(usage=usage_string, version=version_string)

    parser.add_option("-p", "--pdb", metavar='PDBFILE',
            action="store", type="string", dest='pdbfile', default=None,
            help="Input PDB file to serve as the template protein structure.")
    parser.add_option("-s", "--seq", metavar='SEQUENCE',
            action="store", type="string", dest='sequence', default=None,
            help="The amino acid sequence of the chain to be modeled.  \
            This can be either a one-letter code string, or a file containing that string.")

    parser.add_option("-o", "--outdir", metavar='OUTDIR',
            action="store", type="string", dest='outdir', default=None,
            help="Output directory to putr output and project working files.")
    parser.add_option("-i", "--salt", metavar='SALT_TYPE',
            action="store", type="string", dest='salt', default='NaCl',
            help="The type salt ions to add the the solvated box.  Optional.  Default is 'NaCl'")  
    parser.add_option("-c", "--saltconc", metavar='SALT_CONC',
            action="store", dest='saltconc', default=0.0,
            help="The molar (M) concentration (i.e. ionic strength) of salt ions desired.  Optional.  Default is none.") 
    parser.add_option("-a", "--pH", metavar='PHVALUE',
            action="store", dest='pH', default=7.0,
            help="Desired pH of the simulation.  Optional.  Default is 7.0 (Careful with histidines!)")
    parser.add_option("-b", "--boxprotocol", metavar='PROTOCOL',
            action="store", type="string", dest='boxProtocol', default='small',
            help="Keyword for the box size protocol.  Optional.  Default is 'small'. \
            Current options are 'small'=15A from protein edges; 'big'=70A; 'big2'=80A")
    parser.add_option("-f", "--forcefield", metavar='FORCEFIELD',
            action="store", type="string", dest='forcefield', default='ffamber99p', \
            help="Forcefield to build gromacs projects with.  Optional.  Default='ffamber99p'.") 

    # Parse command-line arguments.
    (options,args) = parser.parse_args()

    # Perform minimal error checking.
    if ((not options.pdbfile) or (not options.sequence) or (not options.outdir)):
        parser.error("Must specify the pdbfile, the sequence and the outdir.\n")
   
    pProtein  = PipelineProtein( pdb=options.pdbfile, seq=options.sequence, salt=options.salt, saltconc=float(options.saltconc), pH=float(options.pH), boxProtocol=options.boxProtocol)

    shoveit( pProtein, options.outdir, forcefield=options.forcefield )

