import os
import pdb

import intermol.unit as units

def gromacs_energies(name, top=None, gro=None, in_out='in', gropath='',grosuff=''):
    """

    gropath = path to gromacs binaries
    grosuff = suffix of gromacs binaries, usually '' or '_d'

    """
    mdp = 'Inputs/Gromacs/grompp.mdp'

    if in_out == 'in':
        base = 'Inputs/Gromacs'
        if top == None:
            top = os.path.join(base, name, 'topol.top')
        if gro == None:
            gro = os.path.join(base, name, 'conf.gro')
    elif in_out == 'GtoG':
        base = 'Outputs/GromacsToGromacs'
        if top == None:
            base = os.path.join(base, name, 'topol.top')
        if gro == None:
            base = os.path.join(base, name, 'conf.gro')
    else:
        raise Exception("Unknown flag: {0}".format(in_out))

    tpr  = os.path.join(base, name, 'topol.tpr')
    ener  = os.path.join(base, name, 'ener.edr')
    ener_xvg  = os.path.join(base, name, 'energy.xvg')
    conf  = os.path.join(base, name, 'confout.gro')
    mdout = os.path.join(base, name, 'mdout.mdp')
    state  = os.path.join(base, name, 'state.cpt')
    traj  = os.path.join(base, name, 'traj.trr')
    log  = os.path.join(base, name, 'md.log')

    grompp_bin = os.path.join(gropath, 'grompp' + grosuff)
    mdrun_bin = os.path.join(gropath, 'mdrun' + grosuff)
    genergy_bin = os.path.join(gropath, 'g_energy' + grosuff)

    # grompp'n it up
    os.system(grompp_bin + " -f {mdp} -c {gro} -p {top} -o {tpr} -po {mdout} -maxwarn 1".format(mdp=mdp,
            top=top, gro=gro, tpr=tpr, mdout=mdout))

    # mdrunin'
    os.system(mdrun_bin + " -s {tpr} -o {traj} -cpo {state} -c {conf} -e {ener} -g {log}".format(tpr=tpr,
            traj=traj, state=state, conf=conf, ener=ener, log=log))

    # energizin'
    select = " ".join(map(str, range(1, 15))) + " 0 "
    os.system("echo {select} | ".format(select=select) + genergy_bin + " -f {ener} -o {ener_xvg} -dp".format(ener=ener,
            ener_xvg=ener_xvg))

    # extract g_energy output and parse initial energies
    with open(ener_xvg) as f:
        all_lines = f.readlines()

    # take last line
    sec_last = all_lines[-1].split()[1:]
    data = map(float, sec_last)

    # give everything units
    temp = data[-1] * units.kelvin
    data = [value * units.kilojoules_per_mole for value in data[:-1]]
    data.append(temp)

    # pack it all up in a dictionary
    types = ['Bond', 'Angle', 'Proper Dih.', 'Ryckaert-Bell.', 'LJ-14', 'Coulomb-14',
            'LJ (SR)', 'Disper. corr.', 'Coulomb (SR)', 'Coul. recip.', 'Potential',
            'Kinetic En.', 'Total Energy', 'Temperature']
    e_out = dict(zip(types, data))
    return e_out
