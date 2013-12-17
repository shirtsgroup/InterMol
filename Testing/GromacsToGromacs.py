import sys
import os.path
from optparse import OptionParser

from intermol.System import System
from intermol.GromacsExt.GromacsStructureParser import readStructure, writeStructure
from intermol.GromacsExt.GromacsTopologyParser import GromacsTopologyParser

parser = OptionParser()
parser.add_option('-p', type='str', dest='top', default='',
        help="Topology .top file")
parser.add_option('-c', type='str', dest='gro', default='',
        help="Structure .gro file")
parser.add_option('-n', type='str', dest='name', default='system',
        help="Name of system")

(options, args) = parser.parse_args()
gro = options.gro
gro_file = os.path.join('Inputs/GromacsInputs/', gro)
if not os.path.isfile(gro_file):
    raise Exception("File not found: {0}!".format(gro_file))

top = options.top
top_file = os.path.join('Inputs/GromacsInputs/', top)
if not os.path.isfile(top_file):
    raise Exception("File not found: {0}!".format(top_file))

name = options.name
#--- end of options ---

System._sys = System(name)
print "System initialized"

print "Reading in Gromacs topology '{0}'...".format(top_file)
GromacsTopologyParser._GroTopParser = GromacsTopologyParser()
GromacsTopologyParser._GroTopParser.parseTopology(top_file)

print "Reading in Gromacs structure '{0}'...".format(gro_file)
readStructure(gro_file)

top_out = os.path.join('Outputs/GromacsToGromacs/', name + '.top')
print "Writing out Gromacs topology '{0}'".format(top_out)
GromacsTopologyParser = GromacsTopologyParser()
GromacsTopologyParser.writeTopology(top_out)

gro_out = os.path.join('Outputs/GromacsToGromacs/', name + '.gro')
print "Writing in Gromacs structure '{0}'".format(gro_out)
writeStructure(gro_out)
