import os
import pdb

def gromacs_energies(top, gro, in_out='in'):
    """
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

    os.system("grompp -f {mdp} -c {gro} -p {top} -o {tpr} -po {mdout}".format(mdp=mdp, 
            top=top, gro=gro, tpr=tpr, mdout=mdout))

    os.system("mdrun -s {tpr} -o {traj} -cpo {state} -c {conf} -e {ener} -g {log}".format(tpr=tpr,
            traj=traj, state=state, conf=conf, ener=ener, log=log))

    select = " ".join(map(str, range(1, 15))) + " 0 "
    os.system("echo {select} | g_energy -f {ener} -o {ener_out}".format(select=select,
        ener=ener, ener_out=ener_out))
    with open(ener_out) as f:
        all_lines = f.readlines()

    # take second to last line - last line is after taking a timestep
    sec_last = all_lines[-2].split()
    data = map(float, sec_last)
    return data 
