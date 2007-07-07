#=============================================================================================
# protonate.py 
#
# Generate a protonated PDB file using MCCE
#=============================================================================================
# AUTHORS
#
# Written by Vincent Voelz, Stanford University, 2007-07-06
#=============================================================================================
# REQUIREMENTS
#=============================================================================================
# TODO
# - Need a replacement PDB + dihedral manipulation scheme to replace the ZAM modules
#=============================================================================================
# VERSION CONTROL INFORMATION
__version__ = "$Revision: $"                                                                         
                     
#=============================================================================================
# IMPORTS

import sys, os
import random, math, string, tempfile

from optparse import OptionParser    # For parsing of command line arguments

from mmtools.utilities import Units
import mmtools.mccetools.mcce as mcce 


#=============================================================================================
# FUNCTIONS 
#=============================================================================================


def protonate(inpdbfile, outpdbfile, pH):
    """Generate a protonated PDB file from a template PDB using MCCE."""

    mcceOut = os.path.abspath(outpdbfile)
    prmFile = '../mccetools/prmfiles/run.prm.quick'
    prmFile = os.path.abspath(prmFile)
    mcce.protonatePDB(inpdbfile, mcceOut, pH, os.environ['MCCE_LOCATION'], cleanup=True, prmfile=prmFile)



#===========================================================================
# MAIN
#===========================================================================

if __name__ == '__main__':

    # Create command-line argument options.
    usage_string = """%prog --inpdb PDBFILE --outpdb PDBFILE 

Generate a PDB structure with the correct protonation state using MCCE.
"""

    version_string = "%prog %__version__"

    parser = OptionParser(usage=usage_string, version=version_string)

    parser.add_option("-i", "--inpdb", metavar='PDBFILE',
            action="store", type="string", dest='inpdb', default=None,
            help="A template PDB file to calculate the protonation state from.")
    parser.add_option("-o", "--outpdb", metavar='PDBFILE',
            action="store", type="string", dest='outpdb', default=None,
            help="The output PDBFILE.") 
    parser.add_option("-a", "--pH", metavar='PHVALUE',
            action="store", dest='pH', default=7.0, help="The pH value.  Optional. Default=7.0 ")

    # Parse command-line arguments.
    (options,args) = parser.parse_args()

    # Perform minimal error checking.
    if ((not options.inpdb) or (not options.outpdb)):
        parser.error("Must specify all options.\n")

    # Send arguments to gencoil()
    protonate( options.inpdb, options.outpdb, options.pH )
 
