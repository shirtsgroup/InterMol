# protonateMyPDB.py 
#
# This example script takes the 3gb1_testing.pdb pdbfile as input, and outputs
# a PDB file with the protonation state at a desired pH calculated by MCCE.
#
# In this case, the PDB used for testing is the first 13 residues of 3gb1, the NMR structure
# of the B1 domain from the staphloccocus IgG-binding protein G.
#
# The parameters used in the MCCE run can be found in the prmfile specified below
#
# Vincent Voelz
# May 19 2007
#

import mmtools.mccetools.mcce as mcce
import os, sys

### INPUT LINES ####
#
# Specify the input and output PDB filenames
# NOTE: pdbfile and outpdbfile should be local (not absolute) paths for this project  
pdbfile = '3gb1_testing.pdb'
outpdbfile = os.path.join(os.curdir,'3gb1_out.pdb')

# Specify the pH
pH = 7.0

# Specify a MCCE parameter file with the desired set of parameters for calculating the pKa 
prmfile = '../prmfiles/run.prm.quick'
prmfile = os.path.abspath(prmfile)

# Specify additional or different parameters than prmfile,  if desired.
# xtraprms = {'TITR_PH0':'3.0'}

# Read in the parameters
params = mcce.read_paramfile(prmfile)
mcce.print_prm(params)

# Write output PDB file with the correct protonation state 
### work is done in a temporary dir; Setting cleanup=True will erase these temporary files 
mcce.protonatePDB(pdbfile, outpdbfile, pH, os.environ['MCCE_LOCATION'], cleanup=True, prmfile=prmfile)



