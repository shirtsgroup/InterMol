import pdb
import os
import numpy as np

import intermol.Driver as Driver
#import intermol.convert_md as Driver
from lammps_energies import lammps_energies
from helper_functions import *

def lammps_to_lammps(name, lmppath='', lmpbin='lmp_openmpi', energy=True,
        verbose=False):
    """Test lammps to lammps conversion
    """

    lmp_in = os.path.join('Inputs/Lammps/', name, 'data.lmp')
    if not os.path.isfile(lmp_in):
        raise Exception("File not found: {0}!".format(lmp_in))

    lmp_out = os.path.join('Outputs/LammpsToLammps/', name, 'data.lmp')

    # calc input energies
    if energy:
        e_in = lammps_energies(name, 'in', lmppath, lmpbin, verbose)

    # where the magic happens
    Driver.initSystem(name)
    Driver.load(lmp_in)
    Driver.write(lmp_out)

    # calc output energies
    if energy:
        e_out = lammps_energies(name, 'LtoL', lmppath, lmpbin)

    if energy:
       return combine_energy_results(e_in, e_out)

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-n', type='str', dest='name', default='bmim',
            help="Name of system")
    parser.add_option('-l', type='str', dest='lmppath', default='',
            help="path for LAMMPS binary")
    parser.add_option('-b', type='str', dest='lmpbin', default='lmp_openmpi',
            help="name of for LAMMPS binary")
    parser.add_option('-e', dest='energy', action="store_true",
            help="Evaluate energies",default=False)
    parser.add_option('-v', dest='verbose', action="store_true",
            help="Write LAMMPS output to screen",default=False)



    (options, args) = parser.parse_args()
    name = options.name
    lmppath = options.lmppath
    lmpbin = options.lmpbin
    energy = options.energy
    verbose = options.verbose

    results = lammps_to_lammps(name, lmppath, lmpbin, energy, verbose)

    if energy:
        print_energy_summary(results)
