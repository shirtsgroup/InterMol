from collections import OrderedDict
import logging
import os

import simtk.unit as units

from intermol.utils import which, run_subprocess

AMB_PATH = ''

logger = logging.getLogger('InterMolLog')

# --------- energy evaluation methods ---------- #
key_dict = {'ENERGY': 'Potential',
            'BOND': 'Bond',
            'ANGLE': 'Angle',
            'DIHED': 'All dihedrals',
            '1-4 VDW': 'LJ-14',
            '1-4 EEL': 'Coulomb-14',
            'EEL': 'Coulomb',
            'VDWAALS': 'LJ'
            }


def standardize_key(in_key):
    if in_key in key_dict:
        out_key = key_dict[in_key]
    else:
        out_key = in_key
    return out_key


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

    amber_bin = os.path.join(amb_path, 'sander')
    if not which(amber_bin):
        raise IOError('Unable to find AMBER executable (sander).')

    # run sander
    cmd = [amber_bin, '-i', input, '-c', crd, '-p', prmtop, '-o', mdout, '-O']
    proc = run_subprocess(cmd, 'amber', stdout_path, stderr_path)
    if proc.returncode != 0:
        logger.error('sander failed. See %s' % stderr_path)

    # Extract energies from amber output

    return _group_energy_terms(mdout)


def _group_energy_terms(mdout):
    """Parse AMBER output file to extract and group the energy terms in a dict. """

    with open(mdout) as f:
        all_lines = f.readlines()

    # find where the energy information starts
    i = 0
    for line in all_lines:
        if line[0:8] == '   NSTEP':
            startline = i
            break
        i+=1

    energy_types = []
    energy_values = []

    ranges = [[1,24],[26,49],[51,77]]  # strange ranges for amber file data.

    potential = 0*units.kilocalories_per_mole
    for line in all_lines[startline+3:]:
        if '=' in line:
            for i in range(3):
                r = ranges[i]
                term = line[r[0]:r[1]]
                if '=' in term:
                    energy_type, energy_value = term.split('=')
                    energy_value = float(energy_value)*units.kilocalories_per_mole
                    potential += energy_value
                    energy_type = energy_type.rstrip()
                    if energy_type in key_dict.keys():# remove the whitespace before storing as key.
                        energy_values.append(energy_value)
                        energy_types.append(key_dict[energy_type])
        else:
            break
    energy_types.append('Potential')
    energy_values.append(potential)

    e_out = OrderedDict(zip(energy_types, energy_values))
    # now total up other terms.
    # Dispersive energies.
    dispersive = ['LJ-14', 'LJ']
    e_out['Dispersive'] = 0 * units.kilocalories_per_mole
    for group in dispersive:
        if group in e_out:
            e_out['Dispersive'] += e_out[group]

    # Electrostatic energies.
    electrostatic = ['Coulomb-14', 'Coulomb']
    e_out['Electrostatic'] = 0 * units.kilocalories_per_mole

    for group in electrostatic:
        if group in e_out:
            e_out['Electrostatic'] += e_out[group]

    e_out['Non-bonded'] = e_out['Electrostatic'] + e_out['Dispersive']

    return e_out, mdout
