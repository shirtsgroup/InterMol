# titrateMyPDB.py 
#
# This example script performs a titration on the provided 1b0d.pdb (HEWL) pdbfile used as one of the test proteins for MCCE.
#
# The parameters used in the MCCE run can be found in the prmfile specified below
#
# Vincent Voelz
# May 19 2007
#

import mmtools.mccetools.mcce as mcce
import os, sys

### INPUT LINES ####
# Specify the input and output PDB filenames
# NOTE: pdbfile and outpdbfile should be local (not absolute) paths for this project  
#pdbfile = '1b0d_cleaned.pdb'
pdbfile = 'AAAKAAKAAKAA.pdb'

# Specify the pH
pH = 7.0

# Specify a MCCE parameter file with the desired set of parameters for calculating the pKa 
prmfile = '../prmfiles/run.prm.default'
prmfile = os.path.abspath(prmfile)

# Specify additional or different parameters than prmfile,  if desired.
# xtraprms = {'TITR_PH0':'3.0'}

# Read in the parameters
#params = mcce.read_paramfile(prmfile)
#mcce.print_prm(params)

# Write output PDB file with the correct protonation state 
### work is done in a temporary dir; Setting cleanup=True will erase these temporary files 
#mcce.protonatePDB(pdbfile, outpdbfile, pH, os.environ['MCCE_LOCATION'], cleanup=False, prmfile=prmfile)
pHstart = 0
pHstep = 1.0
pHiters = 14
result = mcce.titrate(pdbfile, pHstart, pHstep, pHiters, os.environ['MCCE_LOCATION'], cleanup=False, prmfile=prmfile)
print result





