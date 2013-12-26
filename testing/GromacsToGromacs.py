import pdb
import os
import numpy as np

import intermol.Driver as Driver
from gromacs_energies import gromacs_energies

def gromacs_to_gromacs(top, gro, name='2PPN', gropath='', grosuff='',
        energy=True, clean=True):
    """Test gromacs to gromacs conversion
    """

    gro_in = os.path.join('Inputs/Gromacs/', gro)
    if not os.path.isfile(gro_in):
        raise Exception("File not found: {0}!".format(gro_in))

    top_in = os.path.join('Inputs/Gromacs/', top)
    if not os.path.isfile(top_in):
        raise Exception("File not found: {0}!".format(top_in))

    top_out = os.path.join('Outputs/GromacsToGromacs/', name + '.top')
    gro_out = os.path.join('Outputs/GromacsToGromacs/', name + '.gro')

    # calc input energies
    if energy:
        e_in = gromacs_energies(top_in, gro_in, 'in', gropath, grosuff)

    # where the magic happens
    Driver.initSystem(name)
    Driver.load(top_in, gro_in)
    Driver.write(top_out, gro_out)

    # calc output energies
    if energy:
        e_out = gromacs_energies(top_out, gro_out, 'GtoG', gropath, grosuff)

    # delete gromacs backup files
    if clean:
        import glob
        filelist = glob.glob("Inputs/Gromacs/#*#")
        filelist += glob.glob("Outputs/GromacsToGromacs/#*#")
        for f in filelist:
            os.remove(f)
    if energy:
        e_in = np.asarray(e_in)
        e_out = np.asarray(e_out)
        diff = e_out - e_in
        rms = np.sqrt(np.mean(diff ** 2))
        return rms, e_in, e_out

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-p', type='str', dest='top', default='2PPN/2PPN.top',
            help="Topology .top file")
    parser.add_option('-c', type='str', dest='gro', default='2PPN/2PPN.gro',
            help="Structure .gro file")
    parser.add_option('-n', type='str', dest='name', default='2PPN',
            help="Name of system")
    parser.add_option('-g', type='str', dest='gropath', default='',
            help="path for GROMACS binary")
    parser.add_option('-s', type='str', dest='grosuff', default='',
            help="suffix for GROMACS binary")
    parser.add_option('-x', type='int', dest='clean', default=1,
            help="Clean backup files produced by GROMACS")
    parser.add_option('-e', type='int', dest='energy', default=1,
            help="Evaluate energies")



    (options, args) = parser.parse_args()
    top = options.top
    gro = options.gro
    name = options.name
    gropath = options.gropath
    grosuff = options.grosuff
    energy = options.energy
    clean = options.clean

    rms, e_in, e_out = gromacs_to_gromacs(top, gro, name, gropath, grosuff, energy, clean)

    print "======================================================================="
    print "Summary statistics"
    types = ['Bond', 'Angle', 'Proper Dih.', 'Ryckaert-Bell.', 'LJ-14', 'Coulomb-14',
            'LJ (SR)', 'Disper. corr.', 'Coulomb (SR)', 'Coul. recip.', 'Potential',
            'Kinetic En.', 'Total Energy', 'Temperature']

    print "%20s %12s %12s %12s" % ("Type", "Input", "Output", "Diff")
    for name, i, o in zip(types, e_in, e_out):
        print "%20s %15.8f %15.8f %15.8f" % (name, i, o, i-o)

    print " "
    print "RMS signed error: %10.5f" % (rms)
