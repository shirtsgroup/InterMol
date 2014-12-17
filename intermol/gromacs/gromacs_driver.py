from collections import OrderedDict
import logging
import os
import subprocess

import simtk.unit as units

import intermol.tests
from intermol.tests.testing_tools import which
from intermol.gromacs.gromacs_parser import load_gromacs, write_gromacs


logger = logging.getLogger('InterMolLog')


def read_file(top_in, gro_in, gropath):
    # Ensure .gro and .top are a valid match.
    tests_path = os.path.dirname(intermol.tests.__file__)
    mdp_path = os.path.join(tests_path, 'gromacs', 'grompp.mdp')
    gromacs_energies(top_in, gro_in, mdp_path, gropath, '',
            grompp_check=True)

    logger.info("Reading Gromacs files '{0}', '{1}'.".format(top_in, gro_in))
    system = load_gromacs(top_in, gro_in)
    logger.info('...loaded.')
    return system


def write_file(system, top_out, gro_out):
    logger.info("Writing Gromacs files '{0}', '{1}'.".format(top_out, gro_out))
    write_gromacs(top_out, gro_out, system)
    logger.info('...done.')


def gromacs_energies(top=None, gro=None, mdp=None,
                     gropath=None, grosuff=None, grompp_check=False):
    """
    gropath = path to gromacs binaries
    grosuff = suffix of gromacs binaries, usually '' or '_d'
    """
    if not grompp_check:
        logger.info('Evaluating energy of {0}'.format(gro))
    if not gropath:
        gropath = ''
    if not grosuff:
        grosuff = ''

    directory, _ = os.path.split(top)

    tpr = os.path.join(directory, 'topol.tpr')
    ener = os.path.join(directory, 'ener.edr')
    ener_xvg = os.path.join(directory, 'energy.xvg')
    conf = os.path.join(directory, 'confout.gro')
    mdout = os.path.join(directory, 'mdout.mdp')
    state = os.path.join(directory, 'state.cpt')
    traj = os.path.join(directory, 'traj.trr')
    log = os.path.join(directory, 'md.log')
    stdout = os.path.join(directory, 'gromacs_stdout.txt')
    stderr = os.path.join(directory, 'gromacs_stderr.txt')

    if which('gmx'):
        grompp_bin = os.path.join(gropath, 'gmx grompp' + grosuff)
        mdrun_bin = os.path.join(gropath, 'gmx mdrun' + grosuff)
        genergy_bin = os.path.join(gropath, 'gmx energy' + grosuff)
    else:
        grompp_bin = os.path.join(gropath, 'grompp' + grosuff)
        mdrun_bin = os.path.join(gropath, 'mdrun' + grosuff)
        genergy_bin = os.path.join(gropath, 'g_energy' + grosuff)

    # grompp'n it up
    cmd = [grompp_bin, '-f', mdp, '-c', gro, '-p', top, '-o', tpr, '-po', mdout, '-maxwarn', '5']
    logger.debug('Running Gromacs with command:\n    %s' % ' '.join(cmd))
    with open(stdout, 'w') as out, open(stderr, 'w') as err:
        exit = subprocess.call(cmd, stdout=out, stderr=err)
    if exit:
        logger.error('grompp failed. See %s' % stderr)
        raise Exception('grompp failed for {0}'.format(gro.split('/')[-1]))
    if grompp_check:
        return

    # mdrunin'
    cmd = [mdrun_bin, '-nt', '1', '-s', tpr, '-o', traj, '-cpo', state, '-c', 
        conf, '-e', ener, '-g', log]
    logger.debug('Running Gromacs with command:\n    %s' % ' '.join(cmd))
    with open(stdout, 'wa') as out, open(stderr, 'a') as err:
        exit = subprocess.call(cmd, stdout=out, stderr=err)
    if exit:
        logger.error('mdrun failed. See %s' % stderr)
        raise Exception('mdrun failed for {0}'.format(gro.split('/')[-1]))

    # energizin'
    select = " ".join(map(str, range(1, 20))) + " 0 "
    cmd = 'echo {select} | {genergy_bin} -f {ener} -o {ener_xvg} -dp'.format(
            select=select, genergy_bin=genergy_bin, ener=ener, ener_xvg=ener_xvg)
    logger.debug('Running Gromacs with command:\n    %s' % cmd)
    with open(stdout, 'wa') as out, open(stderr, 'a') as err:
        exit = subprocess.call(cmd, stdout=out, stderr=err, shell=True)
    if exit:
        logger.error('g_energy failed. See %s' % stderr)
        raise Exception('g_energy failed for {0}'.format(gro.split('/')[-1]))

    # extract g_energy output and parse initial energies
    with open(ener_xvg) as f:
        all_lines = f.readlines()

    types = []
    for line in all_lines:
        if line[:3] == '@ s':
            types.append(line.split('"')[1])

    # take last line
    data = map(float, all_lines[-1].split()[1:])  # [0] is the time

    # give everything units
    data = [value * units.kilojoules_per_mole for value in data]

    # pack it up in a dictionary
    e_out = OrderedDict(zip(types, data))

    # discard non-energy terms
    unwanted = ['Kinetic En.', 'Total Energy', 'Temperature', 'Pressure',
            'Volume', 'Box-X', 'Box-Y', 'Box-atomic_number', 'Pres. DC']
    for group in unwanted:
        if group in e_out:
            del e_out[group]

    # dispersive energies - do buckingham energies also get dumped here?
    dispersive = ['LJ (SR)', 'LJ-14', 'Disper.corr.']
    e_out['Dispersive'] = 0 * units.kilojoules_per_mole
    for group in dispersive:
        if group in e_out:
            e_out['Dispersive'] += e_out[group]

    # electrostatic energies
    electrostatic = ['Coulomb (SR)', 'Coulomb-14', 'Coul. recip.']
    e_out['Electrostatic'] = 0 * units.kilojoules_per_mole
    for group in electrostatic:
        if group in e_out:
            e_out['Electrostatic'] += e_out[group]

    e_out['Non-bonded'] = e_out['Electrostatic'] + e_out['Dispersive']

    # all the various dihedral energies - what else goes in here?
    all_dihedrals = ['Ryckaert-Bell.', 'Proper Dih.', 'Improper Dih.']
    e_out['All dihedrals'] = 0 * units.kilojoules_per_mole
    for group in all_dihedrals:
        if group in e_out:
            e_out['All dihedrals'] += e_out[group]

    return e_out, ener_xvg
