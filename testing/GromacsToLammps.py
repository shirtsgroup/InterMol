import pdb
import os
import numpy as np

import intermol.Driver as Driver
from gromacs_energies import gromacs_energies
from lammps_energies import lammps_energies

def gromacs_to_lammps(name, gropath='', grosuff='', lmppath='', lmpbin='lmp_openmpi',
        energy=True, clean=True):
    """Test gromacs to lammps conversion
    """

    gro_in = os.path.join('Inputs/Gromacs/', name, 'conf.gro')
    if not os.path.isfile(gro_in):
        raise Exception("File not found: {0}!".format(gro_in))

    top_in = os.path.join('Inputs/Gromacs/', name, 'topol.top')
    if not os.path.isfile(top_in):
        raise Exception("File not found: {0}!".format(top_in))

    lmp = os.path.join('Outputs/GromacsToLammps/', name, 'data.lmp')
    inp = os.path.join('Outputs/GromacsToLammps/', name, 'data.input')

    # calc input energies
    if energy:
        e_in = gromacs_energies(name, top_in, gro_in, 'in', gropath, grosuff)

    # where the magic happens
    Driver.initSystem(name)
    Driver.load(top_in, gro_in)
    Driver.write(lmp)

    # calc output energies
    if energy:
        e_out = lammps_energies(name, 'GtoL', lmppath, lmpbin)

    # delete gromacs backup files
    if clean:
        import glob
        filelist = glob.glob("Inputs/Gromacs/{name}/*#".format(name=name))
        for f in filelist:
            os.remove(f)

    if energy:
        diff = dict()
        for e_type, value in e_in.iteritems():
            if e_type in e_out:
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
    parser.add_option('-n', type='str', dest='name', default='system2_GMX',
            help="Name of system")
    parser.add_option('-g', type='str', dest='gropath', default='',
            help="path for GROMACS binary")
    parser.add_option('-s', type='str', dest='grosuff', default='',
            help="suffix for GROMACS binary")
    parser.add_option('-l', type='str', dest='lmppath', default='',
            help="path for LAMMPS binary")
    parser.add_option('-b', type='str', dest='lmpbin', default='lmp_openmpi',
            help="name of for LAMMPS binary")
    parser.add_option('-e', type='int', dest='energy', default=1,
            help="Evaluate energies")
    parser.add_option('-x', type='int', dest='clean', default=1,
            help="Clean backup files produced by GROMACS")
 

    (options, args) = parser.parse_args()
    name = options.name
    gropath = options.gropath
    grosuff = options.grosuff
    lmppath = options.lmppath
    lmpbin = options.lmpbin
    energy = options.energy
    clean = options.clean

    results = gromacs_to_lammps(name, gropath, grosuff, lmppath, lmpbin, energy, clean)

    if energy:
        rms, e_in, e_out = results
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
