import sys
import os
import pdb

import intermol.unit as units

def gromacs_energies(top=None, gro=None, mdp=None, gropath='',grosuff=''):
    """

    gropath = path to gromacs binaries
    grosuff = suffix of gromacs binaries, usually '' or '_d'

    """
    directory, _ = os.path.split(top)

    tpr  = os.path.join(directory , 'topol.tpr')
    ener  = os.path.join(directory , 'ener.edr')
    ener_xvg  = os.path.join(directory , 'energy.xvg')
    conf  = os.path.join(directory , 'confout.gro')
    mdout = os.path.join(directory , 'mdout.mdp')
    state  = os.path.join(directory , 'state.cpt')
    traj  = os.path.join(directory , 'traj.trr')
    log  = os.path.join(directory , 'md.log')

    grompp_bin = os.path.join(gropath, 'grompp' + grosuff)
    mdrun_bin = os.path.join(gropath, 'mdrun' + grosuff)
    genergy_bin = os.path.join(gropath, 'g_energy' + grosuff)

    # grompp'n it up
    cmd = ("{grompp_bin} -f {mdp} -c {gro} -p {top} -o {tpr} -po {mdout} -maxwarn 1".format(
        grompp_bin=grompp_bin, mdp=mdp, top=top, gro=gro, tpr=tpr, mdout=mdout))
    print 'Running GROMACS with command:'
    print cmd
    exit = os.system(cmd)
    if exit:
        print 'grompp failed for {0}'.format(top)
        sys.exit(1)

    # mdrunin'
    cmd = ("{mdrun_bin} -s {tpr} -o {traj} -cpo {state} -c {conf} -e {ener} -g {log}".format(
        mdrun_bin=mdrun_bin, tpr=tpr, traj=traj, state=state,
        conf=conf, ener=ener, log=log))
    print cmd
    exit = os.system(cmd)
    if exit:
        print 'mdrun failed for {0}'.format(top)
        sys.exit(1)

    # energizin'
    select = " ".join(map(str, range(1, 20))) + " 0 "
    cmd = ("echo {select} | ".format(select=select) + "{genergy_bin} -f {ener} -o {ener_xvg} -dp".format(
            genergy_bin=genergy_bin, ener=ener, ener_xvg=ener_xvg))
    if exit:
        print 'g_energy failed for {0}'.format(top)
        sys.exit(1)

    # extract g_energy output and parse initial energies
    with open(ener_xvg) as f:
        all_lines = f.readlines()

    types = []
    for line in all_lines:
        if line[:3] == '@ s':
            types.append(line.split('"')[1])

    # take last line
    data = map(float, all_lines[-1].split()[1:])  # what is [0] of that line?

    # give everything units
    data = [value * units.kilojoules_per_mole for value in data]

    # pack it up in a dictionary
    e_out = dict(zip(types, data))
    # fix temperature unit (or we could just toss it...)
    e_out['Temperature'] = e_out['Temperature']._value * units.kelvin

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

    # non-bonded energies
    e_out['Non-bonded'] = e_out['Electrostatic'] + e_out['Dispersive']

    return e_out, ener_xvg
