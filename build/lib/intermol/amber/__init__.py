from collections import OrderedDict
import logging
import os

import parmed.unit as units

from intermol.exceptions import AmberError
from intermol.utils import which, run_subprocess


AMB_PATH = ''

logger = logging.getLogger('InterMolLog')


to_canonical = {
    'BOND': ['bond'],

    'ANGLE': ['angle'],
    'UB': ['angle','urey-bradley'],

    'DIHED': ['dihedral','proper'],
    'IMP': ['dihedral','improper'],
    'CMAP': ['dihedral','cmap'],

    'HBOND': ['h-bond'],

    'VDWAALS': ['vdw total'],
    '1-4 VDW': ['vdw total', 'vdw-14'],

    'EEL': ['coulomb total'],
    '1-4 EEL': ['coulomb total', 'coulomb-14'],

    'ENERGY': ['potential']
}


def energies(prmtop, crd, input, amb_path):
    """Compute single-point energies using AMBER.

    Args:
        prmtop (str):
        crd (str):
        input (str)
        amb_path (str):

    Returns:
        e_out:
        ener_xvg:
    """

    logger.info('Evaluating energy of {0}'.format(crd))

    directory, _ = os.path.split(os.path.abspath(prmtop))

    if not prmtop:
        prmtop = os.path.join(directory, 'parm.prmtop')
    if not crd:
        crd = os.path.join(directory, 'ener.edr')
    mdout = os.path.join(directory, 'amber.out')
    stdout_path = os.path.join(directory, 'amber_stdout.txt')
    stderr_path = os.path.join(directory, 'amber_stderr.txt')

    # Did they give a path, or the name of the file?
    islastbin = os.path.basename(os.path.normpath(amb_path))
    if islastbin == 'sander':
        amber_bin = amb_path
    else:
        amber_bin = os.path.join(amb_path, 'sander')
    if not which(amber_bin):
        raise IOError('Unable to find AMBER executable (sander).')

    # Run sander.
    cmd = [amber_bin, '-i', input, '-c', crd, '-p', prmtop, '-o', mdout, '-O']
    proc = run_subprocess(cmd, 'amber', stdout_path, stderr_path)
    if proc.returncode != 0:
        logger.error('sander failed. See %s' % stderr_path)

    return _group_energy_terms(mdout)


def _group_energy_terms(mdout):
    """Parse AMBER output file and group the energy terms in a dict. """

    with open(mdout) as f:
        all_lines = f.readlines()

    # Find where the energy information starts.
    for i, line in enumerate(all_lines):
        if line[0:8] == '   NSTEP':
            startline = i
            break
    else:
        raise AmberError('Unable to detect where energy info starts in AMBER '
                         'output file: {}'.format(mdout))

    # Strange ranges for amber file data.
    ranges = [[1, 24], [26, 49], [51, 77]]

    e_out = OrderedDict()
    potential = 0 * units.kilocalories_per_mole
    for line in all_lines[startline+3:]:
        if '=' in line:
            for i in range(3):
                r = ranges[i]
                term = line[r[0]:r[1]]
                if '=' in term:
                    energy_type, energy_value = term.split('=')
                    energy_value = float(energy_value) * units.kilocalories_per_mole
                    potential += energy_value
                    energy_type = energy_type.rstrip()
                    e_out[energy_type] = energy_value
        else:
            break
    e_out['ENERGY'] = potential
    return e_out, mdout
