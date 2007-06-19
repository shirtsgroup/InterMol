#!/usr/bin/python
#=============================================================================================
# modelmissing.py
#
# Model missing atoms and residues from a PDB file.
#
# Written 2007-06-19 by
#
# John D. Chodera <jchodera@stanford.edu>
# Department of Chemistry, Stanford University
# Pande group
#=============================================================================================
# TODO
# - Have tool build all chains if no --chain argument is specified.
#=============================================================================================
"""Command-line tool for modeling missing atoms and residues in a PDB file.
"""

# Set up parser for command-line options.
import optparse
version_string = "Version string goes here."
usage_string = """\
USAGE

%prog --sourcepdb source.pdb [--chainid chainid] --output[db output.pdb

Model missing atoms and residues of the specified chain of <source.pdb>.
If no chain is specified, chainid = ' '.

ARGUMENTS
  --sourcepdb source.pdb - Source PDB file to use as coordinate template and target sequence.
  --chainid chainid - one letter chain ID of chain from source.pdb to model (default ' ')
  --outputpdb output.pdb - name of final output PDB file 

EXAMPLE

  python2.4 ./modelmissing.py --sourcepdb examples/misc/1NHV.pdb --chainid A --outputpdb 1nhv_model.pdb

NOTES

  MODELLER must be in PYTHONPATH.
"""

parser = optparse.OptionParser(usage=usage_string, version=version_string)

parser.add_option("-s", "--sourcepdb", metavar='SOURCEPDB',
                  action="store", type="string", dest='sourcepdb_filename', default = None,
                  help="Source PDB filename.")
parser.add_option("-c", "--chainid", metavar='CHAINID',
                  action="store", type="string", dest='chainid', default = ' ',
                  help="One-letter chain ID.")
parser.add_option("-o", "--outputpdb", metavar='OUTPUTPDB',
                  action="store", type="string", dest='outputpdb_filename', default = None,
                  help="Output PDB filename.")

# Parse command-line arguments.
(options,args) = parser.parse_args()

# Perform minimal error checking.
if (not options.sourcepdb_filename) or (not options.outputpdb_filename):
    parser.print_help()
    parser.error("Both sourcepdb_filename and outputpdb_filename must be specified.")

import modelPDB
import sys

# Initialize modelPDB class.
myModel = modelPDB.ModelPDB()

# Model missing atoms/residues.
myModel.modelMissingAtoms(options.sourcepdb_filename, options.outputpdb_filename, chain=options.chainid, debug = True)


