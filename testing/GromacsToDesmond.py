import os.path
from optparse import OptionParser

import intermol.Driver as Driver


#--- cmd line options ---
parser = OptionParser()
parser.add_option('-p', type='str', dest='top', default='system2_GMX.top',
        help="Topology .top file")
parser.add_option('-c', type='str', dest='gro', default='system2_GMX.gro',
        help="Structure .gro file")
parser.add_option('-n', type='str', dest='name', default='system',
        help="Name of system")

(options, args) = parser.parse_args()
gro = options.gro
gro_in = os.path.join('Inputs/GromacsInputs/', gro)
if not os.path.isfile(gro_in):
    raise Exception("File not found: {0}!".format(gro_in))

top = options.top
top_in = os.path.join('Inputs/GromacsInputs/', top)
if not os.path.isfile(top_in):
    raise Exception("File not found: {0}!".format(top_in))

name = options.name

cms_out = os.path.join('Outputs/GromacsToDesmond/', name + '.cms')
#--- end of options ---

Driver.initSystem(name)

Driver.load(top_in, gro_in)
Driver.write(cms_out)
