import logging
import os
import subprocess

import simtk.unit as units
from intermol.lammps.lammps_parser import load_lammps, write_lammps

logger = logging.getLogger('InterMolLog')


def read_file(in_file):
    logger.info('Reading LAMMPS files {0}'.format(in_file))
    system = load_lammps(in_file)
    logger.info('...loaded.')
    return system


def write_file(in_file, system, unit_set='real'):
    logger.info("Writing LAMMPS file '{0}'".format(in_file))
    write_lammps(in_file, system, unit_set)
    logger.info('...done.')


def lammps_energies(input_file, lmppath='lmp_openmpi'):
    """Evaluate energies of LAMMPS files

    Args:
        input_file = path to input file (expects data file in same folder)
        lmppath = path to LAMMPS binaries
    """
    logger.info('Evaluating energy of {0}'.format(input_file))

    directory, input_file = os.path.split(input_file)

    # mdrunin'
    saved_path = os.getcwd()
    os.chdir(directory)

    cmd = "{lmppath} < {input_file}".format(
            lmppath=lmppath, input_file=input_file)
    logger.debug('Running LAMMPS with command:\n    %s' % cmd)
    with open('lammps_stdout.txt', 'w') as out, open('lammps_stderr.txt', 'w') as err:
        exit = subprocess.call(cmd, stdout=out, stderr=err, shell=True)
    os.chdir(saved_path)
    if exit:
        logger.error('Energy evaluation failed. See %s/lammps_stderr.txt' % directory)
        raise Exception('Energy evaluation failed for {0}'.format(input_file))

    # energizin'
    proc = subprocess.Popen(["awk '/E_bond/{getline; print}' %s/lammps_stdout.txt" % (directory)],
            stdout=subprocess.PIPE, shell=True)
    (energies, err) = proc.communicate()
    if not energies:
        raise Exception('Unable to read LAMMPS energy output')


    # give everything units
    data = map(float, energies.split())
    data = [value * units.kilocalories_per_mole for value in data]

    # pack it all up in a dictionary
    types = ['Bond', 'Angle', 'Proper Dih.', 'Improper', 'Non-bonded',
            'Dispersive', 'Electrostatic', 'Coul. recip.', 'Disper. corr.',
            'Potential']
    e_out = dict(zip(types, data))

    # groupings
    e_out['Electrostatic'] += e_out['Coul. recip.']
    e_out['All dihedrals'] = e_out['Proper Dih.'] + e_out['Improper']
    return e_out, '%s/lammps_stdout.txt' % directory
