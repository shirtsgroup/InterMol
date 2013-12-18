import pdb
import os.path
from optparse import OptionParser
import numpy as np

import intermol.Driver as Driver
from gromacs_energies import gromacs_energies

#--- cmd line options ---
parser = OptionParser()
parser.add_option('-p', type='str', dest='top', default='system2_GMX.top',
        help="Topology .top file")
parser.add_option('-c', type='str', dest='gro', default='system2_GMX.gro',
        help="Structure .gro file")
parser.add_option('-n', type='str', dest='name', default='system',
        help="Name of system")
parser.add_option('-g', type='str', dest='gropath', default='',
        help="path for GROMACS binary")
parser.add_option('-s', type='str', dest='grosuff', default='',
        help="suffix for GROMACS binary")

(options, args) = parser.parse_args()
gro = options.gro
gro_in = os.path.join('Inputs/Gromacs/', gro)
if not os.path.isfile(gro_in):
    raise Exception("File not found: {0}!".format(gro_in))

top = options.top
top_in = os.path.join('Inputs/Gromacs/', top)
if not os.path.isfile(top_in):
    raise Exception("File not found: {0}!".format(top_in))

name = options.name

top_out = os.path.join('Outputs/GromacsToGromacs/', name + '.top')
gro_out = os.path.join('Outputs/GromacsToGromacs/', name + '.gro')
#--- end of options ---

# calc input energies
e_in = gromacs_energies(top_in, gro_in, 'in',
        gropath=options.gropath, grosuff=options.grosuff)

# where the magic happens
Driver.initSystem(name)
Driver.load(top_in, gro_in)
Driver.write(top_out, gro_out)

# calc output energies
e_out = gromacs_energies(top_out, gro_out, 'GtoG',
        gropath=options.gropath, grosuff=options.grosuff)

print "======================================================================="
print "Summary statistics"
types = ['Bond', 'Angle', 'Proper Dih.', 'Ryckaert-Bell.', 'LJ-14', 'Coulomb-14',
        'LJ (SR)', 'Disper. corr.', 'Coulomb (SR)', 'Coul. recip.', 'Potential',
        'Kinetic En.', 'Total Energy', 'Temperature']

e_in = np.asarray(e_in)
e_out = np.asarray(e_out)
print "%20s %12s %12s %12s" % ("Type", "Input", "Output", "Diff")
for name, i, o in zip(types, e_in, e_out):
    print "%20s %12.5f %12.5f %12.5f" % (name, i, o, i-o)

print " "
diff = e_out - e_in
rms_signed = np.sqrt(np.mean(diff ** 2))
print "RMS signed error: %10.5f" % (rms_signed)
