#!/usr/bin/python

#=============================================================================================
# extract-snapshot-to-pdb
#
# Extract a single snapshot INDEX from netcdf trajectory.
#=============================================================================================
# AUTHORS
#
# Written by John Chodera, Stanford University, 2007-03-05
#=============================================================================================
# REQUIREMENTS
# - SQLite
#=============================================================================================
# TODO
#=============================================================================================
# VERSION CONTROL INFORMATION
__version__ = "$Revision: $"                                                                                              
#=============================================================================================
# IMPORTS
import sys
import os.path
from optparse import OptionParser # For parsing of command line arguments
#=============================================================================================

# Create command-line argument options.
usage_string = """
usage: %prog --prmtop PRMTOP_FILE --trajectory TRAJECTORY_FILE --frame FRAME_NUMBER --output OUTPUT_FILE.pdb

example: %prog --prmtop solute.prmtop --trajectory md.netcdf --frame 100 --output extracted.pdb
"""

version_string = "%prog %__version__"

parser = OptionParser(usage=usage_string, version=version_string)

parser.add_option("-p", "--prmtop", metavar='PRMTOP_FILE',
                  action="store", type="string", dest='prmtop_filename', default=None,
                  help="AMBER prmtop filename.")
parser.add_option("-t", "--trajectory", metavar='TRAJECTORY_FILE',
                  action="store", type="string", dest='trajectory_filename', default=None,
                  help="Trajectory file to extract snapshot from.")
parser.add_option("-f", "--frame", metavar='FRAME_NUMBER',
                  action="store", type="int", dest='frame_number', default=None,
                  help="Frame number (starting from 1) to extract PDB file from.")
parser.add_option("-o", "--output", metavar='OUTPUT_FILE',
                  action="store", type="string", dest='output_filename', default=None,
                  help="Name of PDB file to extract file to.")

# Parse command-line arguments.
(options,args) = parser.parse_args()

# Perform minimal error checking.
if ((not options.trajectory_filename) or (not options.frame_number) or (not options.output_filename)):
  parser.print_help()
  parser.error("All options must be specified.\n")

# Create temporary directory.
import tempfile
import os.path
tmpdir = tempfile.mkdtemp()

# Construct ptraj input file.
ptraj_input_filename = os.path.join(tmpdir,'ptraj.in')
ptraj_input = """\
trajin %(trajectory_file)s %(index)d %(index)d\n\
trajout %(output)s pdb\n\
""" % { 'index' : options.frame_number,
        'trajectory_file' : options.trajectory_filename,
        'output' : options.output_filename }

outfile = open(ptraj_input_filename, 'w')
outfile.write(ptraj_input)
outfile.close()

# Form command.
import commands
command = 'ptraj %(prmtop)s < %(ptraj_input_filename)s' % { 'ptraj_input_filename' : ptraj_input_filename,
                                                           'prmtop' : options.prmtop_filename }
print command
output = commands.getoutput(command)
print output

# Move file to proper filenmame
import shutil
source_filename = '%s.%d' % (options.output_filename, options.frame_number)
dest_filename = options.output_filename
shutil.move(source_filename, dest_filename)

# Clean up temporary directory.
for filename in os.listdir(tmpdir):
  os.remove(os.path.join(tmpdir,filename))
os.rmdir(tmpdir)





