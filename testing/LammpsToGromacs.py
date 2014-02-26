import pdb
import os
import numpy as np

import intermol.Driver as Driver
from gromacs_energies import gromacs_energies
from lammps_energies import lammps_energies
from helper_functions import *

def lammps_to_gromacs(name, gropath='', grosuff='', lmppath='', lmpbin='lmp_openmpi',
        energy=True, clean=True):
    """Test lammps to gromacs conversion
    """
    lmp_in = os.path.join('Inputs/Lammps/', name, 'data.lmp')
    if not os.path.isfile(lmp_in):
        raise Exception("File not found: {0}!".format(lmp_in))

    top_out = os.path.join('Outputs/LammpsToGromacs/', name, 'topol.top')
    gro_out = os.path.join('Outputs/LammpsToGromacs/', name, 'conf.gro')

    # calc input energies
    if energy:
        e_in = lammps_energies(name, 'in', lmppath, lmpbin)

    # where the magic happens
    Driver.initSystem(name)
    Driver.load(lmp_in)
    Driver.write(top_out, gro_out)

    # calc output energies
    if energy:
        e_out = gromacs_energies(name, top_out, gro_out, 'LtoG', gropath, grosuff)

    # delete gromacs backup files
    if clean:
        import glob
        filelist = glob.glob("Inputs/Gromacs/{name}/#*#".format(name=name))
        filelist += glob.glob("Outputs/GromacsToGromacs/{name}/#*#".format(name=name))
        for f in filelist:
            os.remove(f)
    if energy:
       return combine_energy_results(e_in, e_out)

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-n', type='str', dest='name', default='bmim',
            help="Name of system")
    parser.add_option('-g', type='str', dest='gropath', default='',
            help="path for GROMACS binary")
    parser.add_option('-s', type='str', dest='grosuff', default='',
            help="suffix for GROMACS binary")
    parser.add_option('-l', type='str', dest='lmppath', default='',
            help="path for LAMMPS binary")
    parser.add_option('-b', type='str', dest='lmpbin', default='lmp_openmpi',
            help="name of for LAMMPS binary")
    parser.add_option('-e', dest='energy', action="store_true",
            help="Evaluate energies",default=False)
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

    results = lammps_to_gromacs(name, gropath, grosuff, lmppath, lmpbin, energy, clean)

    if energy:
        print_energy_summary(results)
