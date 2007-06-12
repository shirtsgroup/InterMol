from mmtools.gromacstools.System import *
import os, sys

# SET UP A GROMACS SIMULATION
pdbfile = 'test.pdb'            # PDB file containing protein to solvate
forcefield = 'ffamber99p'       # gromacs forcefield name to use for grompp
salt = 'NaCl'			# salt pair for counterions - supported types are in the <forcefield>.rtp file
saltconc = 0.150                # Molar units
min_mdpfile = '/User/vincentvoelz/mdp/racecar2_min1.mdp'
equil_mdpfile = '../mdp/racecar2_equil1.mdp'

g = System(pdbfile, useff=forcefield, verbose=True)
g.setup.setSaltConditions(salt, saltconc)


# prepare a system, writing TPR and GRO files, and *mdrun* script in the current directory
thisdir = os.path.abspath( os.curdir )
g.prepare(outname='prod',outdir=thisdir, debug=False, protocol='racecar2', checkForFatalErrors=True)

