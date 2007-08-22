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
from optparse import OptionParser    # For parsing of command line arguments
import pipeline 

#=============================================================================================

# DEBUG flag: if true then save outputs
DEBUG = True


#---------------------------------
# FUNCTIONS
#---------------------------------


#----------------------------------
# MAIN
#----------------------------------


if __name__ == '__main__':

    # Parse command-line argument options.
    version_string = "%prog %__version__"

    parser = OptionParser(version=version_string)


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
    parser.add_option("-t", "--termini_caps", action="store_true", dest="captermini", default=False)

    # Parse command-line arguments.
    (options,args) = parser.parse_args()

    # Perform minimal error checking.
    if ((not options.pdbfile) or (not options.sequence) or (not options.outdir)):
        parser.print_help()
        parser.error("Must specify the pdbfile, the sequence and the outdir.\n")
        
   
    pProtein  = pipeline.PipelineProtein( pdb=options.pdbfile, seq=options.sequence, salt=options.salt, saltconc=float(options.saltconc), pH=float(options.pH), boxProtocol=options.boxProtocol)

    pipeline.shoveit( pProtein, options.outdir, forcefield=options.forcefield, captermini=options.captermini)

