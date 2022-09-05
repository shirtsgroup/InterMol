from collections import OrderedDict
import logging
import os
from subprocess import Popen, PIPE
import warnings

import parmed.unit as units

from intermol.utils import run_subprocess, which
from intermol.lammps.lammps_parser import load, save

# Python 2/3 compatibility.
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = OSError

logger = logging.getLogger('InterMolLog')

to_canonical = {
    'Bond': ['bond'],

    'Angle': ['angle'],

    'Proper Dih.': ['dihedral', 'proper'],
    'Improper': ['dihedral', 'improper'],

    'Dispersive': ['vdw total'],
    'Disper. corr.': ['vdw total', 'vdw (LR)'],
    'Electrostatic': ['coulomb total'],
    'Coul.recip.': ['coulomb total','coulomb (LR)'],

    'Non-bonded': ['nonbonded'],
    'Potential': ['potential']
}


for exe in ['lammps', 'lmp_mpi', 'lmp_serial', 'lmp_openmpi',
            'lmp_mac_mpi', '/home/mish4610/software/lammps/src/lmp_serial']:
    if which(exe):
        LMP_PATH = exe
        break
else:
    LMP_PATH = None


def energies(input_file, lmp_path=None):
    """Evaluate energies of LAMMPS files

    Args:
        input_file = path to input file (expects data file in same folder)
        lmp_path = path to LAMMPS binaries
    """
    if lmp_path is None and LMP_PATH is not None:
        lmp_path = LMP_PATH
    elif LMP_PATH is None:
        raise IOError('Unable to find LAMMPS executables.')

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
    return e_out, stdout_path
