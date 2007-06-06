import mmtools.mccetools.mcce as mcce
import os, sys

### INPUT LINES ####
#
# Specify the input and output PDB filenames
# NOTE: pdbfile and outpdbfile should be local (not absolute) paths for this project  
pdbfile = '1e8l_model1.pdb'
outfile = os.path.join(os.curdir,'titrateHEWL_e8_pK.out')

# Specify the course of the pH titration
pHstart = 1.0
pHstep = 1.0
pHiters = 14

# Specify a MCCE parameter file with the desired set of parameters for calculating the pKa 
prmfile = '../prmfiles/run.prm.quick'
prmfile = os.path.abspath(prmfile)

# Specify additional or different parameters than prmfile,  if desired.
# xtraprms = {}
xtraprms = {'EPSILON_PROT':'8.0'}    # NOTE: only eps=4.0 and eps=8.0 are supported!!!

# Write output PDB file with the correct protonation state 
### work is done in a temporary dir; Setting cleanup=True will erase these temporary files 
mcce.titratePDB(pdbfile, outfile, pHstart, pHstep, pHiters, os.environ['MCCE_LOCATION'], cleanup=False, prmfile=prmfile, xtraprms=xtraprms)



