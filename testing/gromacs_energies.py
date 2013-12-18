import os
import pdb

def gromacs_energies(top, gro, in_out='in', gropath='',grosuff=''):
    """

    gropath = path to gromacs binaries
    grosuff = suffix of gromacs binaries, usually '' or '_d'

    """
    mdp = 'Inputs/Gromacs/grompp.mdp'

    if in_out == 'in':
        base = 'Inputs/Gromacs'
    elif in_out == 'GtoG':
        base = 'Outputs/GromacsToGromacs'
    else:
        raise Exception("Unknown flag: {0}".format(in_out))

    tpr  = os.path.join(base, 'topol.tpr')
    ener  = os.path.join(base, 'ener.edr')
    ener_out  = os.path.join(base, 'energy.xvg')
    conf  = os.path.join(base, 'confout.gro')
    mdout = os.path.join(base, 'mdout.mdp')
    state  = os.path.join(base, 'state.cpt')
    traj  = os.path.join(base, 'traj.trr')
    log  = os.path.join(base, 'md.log')

    grompp_bin = os.path.join(gropath, 'grompp' + grosuff)
    mdrun_bin = os.path.join(gropath, 'mdrun' + grosuff)
    genergy_bin = os.path.join(gropath, 'g_energy' + grosuff)

    # grompp'n it up
    os.system(grompp_bin + " -f {mdp} -c {gro} -p {top} -o {tpr} -po {mdout}".format(mdp=mdp,
            top=top, gro=gro, tpr=tpr, mdout=mdout))

    # mdrunin'
    os.system(mdrun_bin + " -s {tpr} -o {traj} -cpo {state} -c {conf} -e {ener} -g {log}".format(tpr=tpr,
            traj=traj, state=state, conf=conf, ener=ener, log=log))

    # energizin'
    select = " ".join(map(str, range(1, 15))) + " 0 "
    os.system("echo {select} | ".format(select=select) + genergy_bin + " -f {ener} -o {ener_out} -dp".format(ener=ener,
            ener_out=ener_out))

    # extract g_energy output and parse initial energies
    with open(ener_out) as f:
        all_lines = f.readlines()

    # take last line
    sec_last = all_lines[-1].split()
    data = map(float, sec_last)
    return data
