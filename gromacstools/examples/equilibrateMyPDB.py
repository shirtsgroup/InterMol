# EquilibrateMyPDB.py 
#
# This example script takes the pdbfile "test.pdb" as input, and outputs
# everything you need to perform a GROMACS molecular dynamics equilibration:
#
#     GROMACS *.tpr
#     GROMACS *.gro
#     ./equilibrate        script that calls mdrun
#     prepare.log          logfile
#
# In this case, the PDB used for testing is a the first 13 residues of 3gb1, which has been
# pre-"puurified" by the mccetools scriptss, with the correct atom and residue names for the AMBER forcefields ports
# The PDB file is the NMR structure
# of the B1 domain from the staphloccocus IgG-binding protein G.  The protonation state for
# this test PDB has a net charge of +1, but the correct number of counterions consistent with
# the salt concentration are automatically added to form a neutral system
#
#
# Vincent Voelz
# May 21, 2007
#
import mmtools.gromacstools.system as system
import os, sys

# SET UP A GROMACS SIMULATION
pdbfile = 'test.pdb'            # PDB file containing protein to solvate
forcefield = 'ffamber99p'       # gromacs forcefield name to use for grompp
salt = 'NaCl'			# salt pair for counterions - supported types are in the <forcefield>.rtp file
saltconc = 0.150		# salt concentration in solution (Molarity) # JDC: Should use Units class for this after modifying ioncalc.py to use Units too.

g = system.GromacsSystem(pdbfile, useff=forcefield)
# g.setup.setSaltConditions('NaCl', 0.150 )
# JDC: Presumably you meant to use the variables you defined up there.
g.setup.setSaltConditions(salt, saltconc)

# prepare a system, writing TPR and GRO files, and *mdrun* script in the current directory
thisdir = os.path.abspath( os.curdir )
g.prepare(outname='prod',outdir=thisdir)

