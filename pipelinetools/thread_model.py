#=============================================================================================
# thread.py 
#
# Threads a protein sequence onto a template PDB structure
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
import mmtools.modellertools.modelPDB as modelPDB 


#=============================================================================================
# FUNCTIONS 
#=============================================================================================


def thread_model(pdbTemplate, sequence, outPdbFile):


    # use either a sequence file or a sequence string
    if os.path.exists(sequence):
        sequenceFilename = sequence
    else:
        tmpdir = tempfile.mkdtemp()
        sequenceFilename = os.path.join(tmpdir,'sequence')
        fseq = open(sequenceFilename,'w')
        fseq.write(sequence)
        fseq.close()

    myModel = modelPDB.ModelPDB()
    myModel.makeModel(pdbTemplate, sequenceFilename, outPdbFile)
 
    # cleanup
    if not os.path.exists(sequence):
        os.remove(sequenceFilename)
        os.rmdir(tmpdir)

#===========================================================================
# MAIN
#===========================================================================

if __name__ == '__main__':

    # Create command-line argument options.
    usage_string = """%prog --pdb PDB_TEMPLATE --seq SEQUENCE --outpdb PDBFILENAME 

Thread a template PDB sturcture with the sequence given.

REQUIREMENTS

Requires that the Python modules for MODELLER (/Library/modeller-9v1/modlib)
are in your PYTHONPATH.
 
"""

    version_string = "%prog %__version__"

    parser = OptionParser(usage=usage_string, version=version_string)

    parser.add_option("-p", "--pdb", metavar='PDB_TEMPLATE',
            action="store", type="string", dest='pdb_template', default=None,
            help="A template PDB file to serve as a template for the model.")
    parser.add_option("-s", "--seq", metavar='SEQUENCE',
            action="store", type="string", dest='sequence', default=None,
            help="The amino acid sequence to be threaded onto the template structure.  \
            This can be either a one-letter code string, or a file containing that string.")
    parser.add_option("-o", "--outpdb", metavar='PDBFILENAME',
            action="store", type="string", dest='outpdb', default=None,
            help="The output filename of the threaded PDB structure.")

    # Parse command-line arguments.
    (options,args) = parser.parse_args()

    # Perform minimal error checking.
    if ((not options.pdb_template) or (not options.sequence) or (not options.outpdb)):
        parser.error("Must specify all options.\n")

    # Send arguments to gencoil()
    thread( options.pdb_template, options.sequence, options.outpdb )
 
