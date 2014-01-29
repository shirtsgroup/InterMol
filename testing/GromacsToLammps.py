import pdb
import os
import numpy as np

import intermol.Driver as Driver
from gromacs_energies import gromacs_energies
#from lammps_energies import lammpss_energies

def gromacs_to_lammps(top, gro, name='system2_GMX', gropath='', grosuff='',
        energy=True, clean=True):
    """Test gromacs to lammps conversion
    """

    gro_in = os.path.join('Inputs/Gromacs/', gro)
    if not os.path.isfile(gro_in):
        raise Exception("File not found: {0}!".format(gro_in))

    top_in = os.path.join('Inputs/Gromacs/', top)
    if not os.path.isfile(top_in):
        raise Exception("File not found: {0}!".format(top_in))

    lmp_out = os.path.join('Outputs/GromacsToLammps/', '{0}.lmp'.format(name))

    # calc input energies
    if energy:
        e_in = gromacs_energies(top_in, gro_in, 'in', gropath, grosuff)

    # where the magic happens
    Driver.initSystem(name)
    Driver.load(top_in, gro_in)
    Driver.write(lmp_out)

    # calc output energies
    if energy:
        e_out = lammps_energies(lmp_out, inp_out, 'GtoL')

    # delete gromacs backup files
    if clean:
        import glob
        filelist = glob.glob("Inputs/Gromacs/#*#")
        #filelist += glob.glob("Outputs/LammpsToGromacs/*.log")
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
    parser.add_option('-p', type='str', dest='top', default='system2_GMX/system2_GMX.top',
            help="Topology .top file")
    parser.add_option('-c', type='str', dest='gro', default='system2_GMX/system2_GMX.gro',
            help="Structure .gro file")
    parser.add_option('-n', type='str', dest='name', default='system2_GMX',
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

    results = gromacs_to_lammps(top, gro, name, gropath, grosuff, energy, clean)

    if energy:
        rms, e_in, e_out = results
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