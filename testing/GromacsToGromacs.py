import pdb
import os
import numpy as np

import intermol.Driver as Driver
from gromacs_energies import gromacs_energies

def gromacs_to_gromacs(name='2PPN', top=None, gro=None, gropath='', grosuff='',
        energy=True, clean=True):
    """Test gromacs to gromacs conversion
    """

    if top == None:
        top = os.path.join(name,'topol.top')

    if gro == None:
        gro = os.path.join(name,'conf.gro')

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
        e_in = gromacs_energies(name, top_in, gro_in, 'in', gropath, grosuff)

    # where the magic happens
    Driver.initSystem(name)
    Driver.load(top_in, gro_in)
    Driver.write(top_out, gro_out)

    # calc output energies
    if energy:
        e_out = gromacs_energies(name, top_out, gro_out, 'GtoG', gropath, grosuff)

    # delete gromacs backup files
    if clean:
        import glob
        filelist = glob.glob("Inputs/Gromacs/{name}/#*#".format(name=name))
        filelist += glob.glob("Outputs/GromacsToGromacs/{name}/#*#".format(name=name))
        for f in filelist:
            os.remove(f)
    if energy:
        diff = dict()
        for e_type, value in e_in.iteritems():
            diff[e_type] = value - e_out[e_type]

        diff_num = list()
        for value in diff.values():
            diff_num.append(value._value)
        diff_num = np.asarray(diff_num)
        rms = np.sqrt(np.mean(diff_num ** 2))
        return rms, e_in, e_out

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-n', type='str', dest='name', default='2PPN',
            help="Name of system")
    parser.add_option('-p', type='str', dest='top', default=None,
            help="Topology .top file")
    parser.add_option('-c', type='str', dest='gro', default=None,
            help="Structure .gro file")
    parser.add_option('-g', type='str', dest='gropath', default='',
            help="path for GROMACS binary")
    parser.add_option('-s', type='str', dest='grosuff', default='',
            help="suffix for GROMACS binary")
    parser.add_option('-x', type='int', dest='clean', default=1,
            help="Clean backup files produced by GROMACS")
    parser.add_option('-e', type='int', dest='energy', default=1,
            help="Evaluate energies")



    (options, args) = parser.parse_args()
    name = options.name
    top = options.top
    gro = options.gro
    gropath = options.gropath
    grosuff = options.grosuff
    energy = options.energy
    clean = options.clean

    rms, e_in, e_out = gromacs_to_gromacs(name,top, gro, gropath, grosuff, energy, clean)

    print "======================================================================="
    print "Summary statistics"
    types = ['Bond', 'Angle', 'Proper Dih.', 'Ryckaert-Bell.', 'LJ-14', 'Coulomb-14',
            'LJ (SR)', 'Disper. corr.', 'Coulomb (SR)', 'Coul. recip.', 'Potential',
            'Kinetic En.', 'Total Energy', 'Temperature']

    print "%20s %12s %12s %12s" % ("Type", "Input", "Output", "Diff")
    for name, i, o in zip(types, e_in.values(), e_out.values()):
        print "%20s %15.8f %15.8f %15.8f" % (name, i._value, o._value, (i-o)._value)

    print " "
    print "RMS signed error: %10.5f" % (rms)
