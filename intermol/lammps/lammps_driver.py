from collections import OrderedDict
import logging
import os
from subprocess import Popen, PIPE

import simtk.unit as units
from intermol.lammps.lammps_parser import load_lammps, write_lammps
from intermol.tests.testing_tools import run_subprocess, which

# Python 2/3 compatibility.
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = OSError

logger = logging.getLogger('InterMolLog')


def read_file(in_file):
    logger.info('Reading LAMMPS files {0}'.format(in_file))
    system = load_lammps(in_file)
    logger.info('...loaded.')
    return system


def write_file(in_file, system, unit_set='real', nonbonded_style = None):
    logger.info("Writing LAMMPS file '{0}'".format(in_file))
    write_lammps(in_file, system, unit_set, nonbonded_style = nonbonded_style)
    logger.info('...done.')


def lammps_energies(input_file, lmp_path='lmp_openmpi'):
    """Evaluate energies of LAMMPS files

    Args:
        input_file = path to input file (expects data file in same folder)
        lmp_path = path to LAMMPS binaries
    """

    # look at potential paths
    if not lmp_path:
        for exe in ['lmp_openmpi', 'lmp_mpi', 'lmp_serial']:
            if which(exe):
                lmp_path = exe
                break
        else:
            logger.exception('Found no LAMMPS executable.')

    logger.info('Evaluating energy of {0}'.format(input_file))
    directory, input_file = os.path.split(os.path.abspath(input_file))
    stdout_path = os.path.join(directory, 'lammps_stdout.txt')
    stderr_path = os.path.join(directory, 'lammps_stderr.txt')
    # TODO: Read energy info from stdout in memory instead of from log files.
    try:
        os.remove(stdout_path)
    except FileNotFoundError:
        pass
    try:
        os.remove(stderr_path)
    except FileNotFoundError:
        pass

    # Step into the directory.
    saved_path = os.getcwd()
    os.chdir(directory)

    cmd = [lmp_path, '-in', input_file]
    proc = run_subprocess(cmd, 'lammps', stdout_path, stderr_path)
    if proc.returncode != 0:
        logger.error('LAMMPS failed. See %s/lammps_stderr.txt' % directory)

    # Step back out.
    os.chdir(saved_path)

    return _group_energy_terms(stdout_path)


def _group_energy_terms(stdout_path):
    """Parse LAMMPS stdout to extract and group the energy terms in a dict. """
    proc = Popen(["awk '/E_bond/{getline; print}' %s" % stdout_path], stdout=PIPE, shell=True)
    energies, err = proc.communicate()
    if not energies:
        raise Exception('Unable to read LAMMPS energy output')

    energy_values = [float(x) * units.kilocalories_per_mole for x in energies.split()]
    energy_types = ['Bond', 'Angle', 'Proper Dih.', 'Improper', 'Non-bonded',
                    'Dispersive', 'Electrostatic', 'Coul. recip.',
                    'Disper. corr.', 'Potential']
    e_out = OrderedDict(zip(energy_types, energy_values))

    e_out['Electrostatic'] += e_out['Coul. recip.']
    e_out['All dihedrals'] = e_out['Proper Dih.'] + e_out['Improper']
    return e_out, stdout_path
