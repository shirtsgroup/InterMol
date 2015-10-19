from collections import OrderedDict
import logging
import os
import subprocess

import simtk.unit as units

import intermol.tests
from intermol.tests.testing_tools import which, run_subprocess

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


def amber_energies(prmtop, crd, input, amb_path):
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

    i = 0
    import pdb
    pdb.set_trace()
    for line in all_lines:
        if line[0:8] == '   NSTEP':
            startline = i
            break
    sp = 24
    energy_types = []
    energy_values = []
    potential = 0*units.kilocalories_per_mole
    for line in all_lines[startline+3,:]:
        if '=' in line:
            for i in range(3):
                term = line[sp*i,sp*(i+1)]
                if '=' in term:
                    energy_type, energy_value = line.split('=')
                    energy_value = float(energy_value)*units.kilocalories_per_mole
                    potential += energy_value
                    energy_values.append(energy_value)
                    if energy_type in key_dict.keys():
                        energy_types.append(key_dict[energy_type])
    energy_types.append('Potential')
    energy_values.append(potential)
    
    # Dispersive energies.
    dispersive = ['1-4 VDW', 'VDWAALS']
    e_out['Dispersive'] = 0 * units.kilocalories_per_mole
    for group in dispersive:
        if group in e_out:
            e_out['Dispersive'] += e_out[group]

    # Electrostatic energies.
    electrostatic = ['1-4 EEL', 'EEL']
    e_out['Electrostatic'] = 0 * units.kilocalories_per_mole
    for group in electrostatic:
        if group in e_out:
            e_out['Electrostatic'] += e_out[group]

    e_out['Non-bonded'] = e_out['Electrostatic'] + e_out['Dispersive']

    # All the various dihedral energies.
    # TODO: What else goes in here?
    all_dihedrals = ['DIHED']
    e_out['All dihedrals'] = 0 * units.kilocalories_per_mole
    for group in all_dihedrals:
        if group in e_out:
            e_out['All dihedrals'] += e_out[group]

    e_out = OrderedDict(zip(energy_types, energy_values))

    return e_out, mdout
