#=============================================================================================
# Find the maximal common GAFF substructure among a number of different ligands and write
# this structure to a mol2 file.
# 
# Written by John D. Chodera <jchodera@gmail.com> 2008-02-08
#=============================================================================================

#=============================================================================================
# IMPORTS
#=============================================================================================

from openeye.oechem import *
from mmtools.moltools.ligandtools import *
from mmtools.moltools.relativefeptools import *
import os
import os.path
from numpy import *

#=============================================================================================
# PARAMETERS
#=============================================================================================

# base path for all parameterized ligand files
ligand_basepath = 'ligands-parameterized/'

# Target ligand group.
ligand_name_list = [ ('LG%d-0-eq' % index) for index in (2, 5, 6, 8, 9) ]

# Filename of mol2 file for generated intermediate structure to write.
intermediate_filename = 'intermediate.mol2'

#=============================================================================================
# SUBROUTINES
#=============================================================================================

#=============================================================================================
# MAIN
#=============================================================================================

# read all ligands
ligands = list()
for ligand_name in ligand_name_list:
    # attempt to load the ligand with GAFF parameters
    ligand = loadGAFFMolecule(ligand_basepath, ligand_name)
    # skip if no ligand found
    if not ligand: continue
    # append to list
    ligands.append(ligand)
print "%d ligands loaded" % len(ligands)

# find common substructure
common_substructure = determineCommonSubstructure(ligands, debug = True)

# find RMS-fit charges
common_substructure = determineMinimumRMSCharges(common_substructure, ligands, debug = True)

# write out molecule
writeMolecule(common_substructure, intermediate_filename, preserve_atomtypes = True)

