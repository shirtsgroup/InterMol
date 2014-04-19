import os
import subprocess

import intermol.unit as units

import pdb

def lammps_energies(input_file, lmppath='lmp_openmpi',
        verbose=False):
    """Evaluate energies of LAMMPS files

    Args:
        input_file = path to input file (expects data file in same folder)
        lmppath = path to LAMMPS binaries
    """

    directory, input_file = os.path.split(input_file)
    log = os.path.join(directory, 'log.lammps')

    # mdrunin'
    saved_path = os.getcwd()
    os.chdir(directory)

    run_lammps = "{lmppath} < {input_file}".format(
            lmppath=lmppath, input_file=input_file)

    os.system(run_lammps)
    os.chdir(saved_path)

    # energizin'
    proc = subprocess.Popen(["awk '/E_bond/{getline; print}' %s" % (log)],
            stdout=subprocess.PIPE, shell=True)
    (energies, err) = proc.communicate()
    if not energies:
        raise Exception("Unable to read LAMMPS energy output")


    # give everything units
    data = map(float, energies.split())
    data = [value * units.kilocalories_per_mole for value in data]

    # pack it all up in a dictionary
    types = ['Bond', 'Angle', 'Proper Dih.', 'Improper', 'Non-bonded',
            'Dispersive', 'Electrostatic', 'Coul. recip.', 'Disper. corr.',
            'Potential']
    e_out = dict(zip(types, data))

    #groupings
    #e_out['Electrostatic'] += e_out['Coul. recip.']
    e_out['All dihedrals'] = e_out['Proper Dih.'] + e_out['Improper']
    return e_out, log
