import os
import subprocess

import intermol.unit as units

import pdb

def lammps_energies(name, in_out='in', lmppath='', lmpbin='lmp_openmpi'):
    """Evaluate energies of LAMMPS files

    Args:
        lmppath = path to LAMMPS binaries
        lmpbin = name of LAMMPS binary
    """

    if in_out == 'in':
        base = 'Inputs/Lammps'
    elif in_out == 'GtoL':
        base = 'Outputs/GromacsToLammps'
    elif in_out == 'LtoL':
        base = 'Outputs/LammpsToLammps'
    else:
        raise Exception("Unknown flag: {0}".format(in_out))

    lmpbin = os.path.join(lmppath, lmpbin)
    sim_dir = os.path.join(base, name)
    log = os.path.join(sim_dir, 'log.lammps')

    # mdrunin'
    saved_path = os.getcwd()
    os.chdir(sim_dir)
    run_lammps = "{lmpbin} < data.input".format(lmpbin=lmpbin)
    #run_lammps = "{lmpbin} < input_file.out".format(lmpbin=lmpbin)
    os.system(run_lammps)
    os.chdir(saved_path)

    # energizin'
    proc = subprocess.Popen(["awk '/E_bond/{getline; print}' %s" % (log)],
            stdout=subprocess.PIPE, shell=True)
    (energies, err) = proc.communicate()

    data = map(float, energies.split())

    # give everything units
    #temp = data[-1] * units.kelvin
    data = [value * units.kilocalories_per_mole for value in data]
    #data.append(temp)

    # pack it all up in a dictionary
    types = ['Bond', 'Angle', 'Proper Dih.', 'Improper', 'Pairs', 'vdW', 'Coulomb', 'Potential']

    e_out = dict(zip(types, data))
    return e_out
