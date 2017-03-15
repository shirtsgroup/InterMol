from collections import OrderedDict
import logging
import os

import parmed.unit as units

from intermol.utils import which, run_subprocess
from intermol.gromacs.gromacs_parser import load, save

GMX_PATH = ''

logger = logging.getLogger('InterMolLog')


def binaries(gmxpath, gmxsuff):
    """Locate the paths to the best available gromacs binaries. """
    def gmx_path(binary_path):
        return os.path.join(gmxpath, binary_path + gmxsuff)

    if which('gmx_d'):
        logger.debug("Using double precision binaries for gromacs")
        main_binary = gmx_path('gmx_d')
        grompp_bin = [main_binary, 'grompp']
        mdrun_bin = [main_binary, 'mdrun']
        genergy_bin = [main_binary, 'energy']
    elif which('grompp_d') and which('mdrun_d') and which('g_energy_d'):
        logger.debug("Using double precision binaries")
        grompp_bin = [gmx_path('grompp_d')]
        mdrun_bin = [gmx_path('mdrun_d')]
        genergy_bin = [gmx_path('g_energy_d')]
    elif which('gmx'):
        logger.debug("Using double precision binaries")
        main_binary = gmx_path('gmx')
        grompp_bin = [main_binary, 'grompp']
        mdrun_bin = [main_binary, 'mdrun']
        genergy_bin = [main_binary, 'energy']
    elif which('grompp') and which('mdrun') and which('g_energy'):
        logger.debug("Using single precision binaries")
        grompp_bin = [gmx_path('grompp')]
        mdrun_bin = [gmx_path('mdrun')]
        genergy_bin = [gmx_path('g_energy')]
    else:
        raise IOError('Unable to find gromacs executables.')
    return grompp_bin, mdrun_bin, genergy_bin

# energy terms we are ignoring
unwanted = ['Kinetic En.', 'Total Energy', 'Temperature', 'Pressure',
            'Volume', 'Box-X', 'Box-Y', 'Box-Z', 'Box-atomic_number',
            'Pres. DC', 'Vir-XY', 'Vir-XX', 'Vir-XZ', 'Vir-YY', 'Vir-YX',
            'Vir-YZ', 'Vir-ZX', 'Vir-ZY', 'Vir-ZZ', 'pV', 'Density', 'Enthalpy']


def energies(top, gro, mdp, gmx_path=GMX_PATH, grosuff='', grompp_check=False):
    """Compute single-point energies using GROMACS.

    Args:
        top (str):
        gro (str):
        mdp (str):
        gmx_path (str):
        grosuff (str):
        grompp_check (bool):

    Returns:
        e_out:
        ener_xvg:
    """

    if not os.path.isfile(mdp):
        logger.error("Can't find mdp file %s to compute energies" % (mdp))
    mdp = os.path.abspath(mdp)

    directory, _ = os.path.split(os.path.abspath(top))

    tpr = os.path.join(directory, 'topol.tpr')
    ener = os.path.join(directory, 'ener.edr')
    ener_xvg = os.path.join(directory, 'energy.xvg')
    conf = os.path.join(directory, 'confout.gro')
    mdout = os.path.join(directory, 'mdout.mdp')
    state = os.path.join(directory, 'state.cpt')
    traj = os.path.join(directory, 'traj.trr')
    log = os.path.join(directory, 'md.log')
    stdout_path = os.path.join(directory, 'gromacs_stdout.txt')
    stderr_path = os.path.join(directory, 'gromacs_stderr.txt')

    grompp_bin, mdrun_bin, genergy_bin = binaries(gmx_path, grosuff)

    # Run grompp.
    grompp_bin.extend(['-f', mdp, '-c', gro, '-p', top, '-o', tpr, '-po', mdout, '-maxwarn', '5'])
    proc = run_subprocess(grompp_bin, 'gromacs', stdout_path, stderr_path)
    if proc.returncode != 0:
        logger.error('grompp failed. See %s' % stderr_path)

    # Run single-point calculation with mdrun.
    mdrun_bin.extend(['-nt', '1', '-s', tpr, '-o', traj, '-cpo', state, '-c', conf, '-e', ener, '-g', log])
    proc = run_subprocess(mdrun_bin, 'gromacs', stdout_path, stderr_path)
    if proc.returncode != 0:
        logger.error('mdrun failed. See %s' % stderr_path)

    # Extract energies using g_energy
    select = " ".join(map(str, range(1, 20))) + " 0 "
    genergy_bin.extend(['-f', ener, '-o', ener_xvg, '-dp'])
    proc = run_subprocess(genergy_bin, 'gromacs', stdout_path, stderr_path, stdin=select)
    if proc.returncode != 0:
        logger.error('g_energy failed. See %s' % stderr_path)

    return _group_energy_terms(ener_xvg)


def _group_energy_terms(ener_xvg):
    """Parse energy.xvg file to extract and group the energy terms in a dict. """
    with open(ener_xvg) as f:
        all_lines = f.readlines()
    energy_types = [line.split('"')[1] for line in all_lines if line[:3] == '@ s']
    energy_values = [float(x) * units.kilojoule_per_mole for x in all_lines[-1].split()[1:]]
    e_out = OrderedDict(zip(energy_types, energy_values))

    # Discard non-energy terms.
    for group in unwanted:
        if group in e_out:
            del e_out[group]

    # Dispersive energies.
    # TODO: Do buckingham energies also get dumped here?
    dispersive = ['LJ (SR)', 'LJ-14', 'Disper. corr.']
    e_out['Dispersive'] = 0 * units.kilojoules_per_mole
    for group in dispersive:
        if group in e_out:
            e_out['Dispersive'] += e_out[group]

    # Electrostatic energies.
    electrostatic = ['Coulomb (SR)', 'Coulomb-14', 'Coul. recip.']
    e_out['Electrostatic'] = 0 * units.kilojoules_per_mole
    for group in electrostatic:
        if group in e_out:
            e_out['Electrostatic'] += e_out[group]

    e_out['Non-bonded'] = e_out['Electrostatic'] + e_out['Dispersive']

    all_angles = ['Angle', 'U-B', 'G96Angle', 'Restricted Angles', 'Bond-Cross',
                  'BA-Cross', 'Quartic Angles']
    e_out['All angles'] = 0 * units.kilojoules_per_mole
    for group in all_angles:
        if group in e_out:
            e_out['All angles'] += e_out[group]

    all_dihedrals = ['Ryckaert-Bell.', 'Proper Dih.', 'Improper Dih.']
    e_out['All dihedrals'] = 0 * units.kilojoules_per_mole
    for group in all_dihedrals:
        if group in e_out:
            e_out['All dihedrals'] += e_out[group]

    return e_out, ener_xvg

