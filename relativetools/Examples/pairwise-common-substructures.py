#=============================================================================================
# Find pairwise common substructures between ligands and report number of shared atoms.
# 
# Written by John D. Chodera <jchodera@gmail.com> 2008-02-08
#=============================================================================================

#=============================================================================================
# IMPORTS
#=============================================================================================

from openeye.oechem import *
from mmtools.moltools.ligandtools import *
from mmtools.moltools.relativefeptools import *
from mmtools.gromacstools.MdpFile import *
import commands
import os
import os.path
import shutil
from numpy import *

#=============================================================================================
# PARAMETERS
#=============================================================================================

# base path for all parameterized ligand files
ligand_basepath = 'ligands-parameterized/'

# list of ligands to consider
ligand_name_list = [ ('jnk.aff-%d' % index) for index in (13, 16, 53) ]

#=============================================================================================
# SUBROUTINES
#=============================================================================================

#=============================================================================================
# MAIN
#=============================================================================================

# Consider all ligand pairs.
nligands = len(ligand_name_list)

#for ligand1_index in range(nligands):
for ligand1_index in range(len(ligand_name_list)):
    # get ligand name
    ligand1_name = ligand_name_list[ligand1_index]
    
    # Load molecule with GAFF atom names and integer bondtypes.
    ligand1 = loadGAFFMolecule(ligand_basepath, ligand1_name)
    if not ligand1: continue     
    
    # Create an OEMCSSearch from this molecule.
    mcss = OEMCSSearch(ligand1, OEExprOpts_StringType, OEExprOpts_IntType)
    # ignore substructures smaller than 4 atoms
    mcss.SetMinAtoms(4)

    # Consider all other ligands.
    for ligand2_index in range(len(ligand_name_list)):
        # skip self
        if ligand1_index == ligand2_index: continue
        
        # get name
        ligand2_name = ligand_name_list[ligand2_index]        

        # Load the molecule
        ligand2 = loadGAFFMolecule(ligand_basepath, ligand2_name)        
        if not ligand2: continue

        # SHOW MATCH
        #print '%(ligand1_name)s : %(ligand2_name)s' % vars()

        # Compare atomsets
        matchcount = 0
        maxmatch = 0
        for match in mcss.Match(ligand2):
#            print "Match:", matchcount, "Size:", match.NumAtoms(), " atoms"
#
#            # determine mutated charges and atomtypes
#            (mutated_charges_1, mutated_charges_2) = determineMutatedCharges(ligand1, ligand2, match)            
#            (mutated_atomtypes_1, mutated_atomtypes_2) = determineMutatedAtomtypes(ligand1, ligand2, match)
#
#            # DEBUG
#            print "ligand1: %s" % ligand1.GetTitle()
#            for atom in ligand1.GetAtoms():
#                atomname = atom.GetName()
#                print "%6s : %16s %6.3f -> %16s %6.3f" % (atom.GetName(), atom.GetType(), atom.GetPartialCharge(), mutated_atomtypes_1[atomname], mutated_charges_1[atomname])
#
#            print "ligand2: %s" % ligand2.GetTitle()
#            for atom in ligand2.GetAtoms():
#                atomname = atom.GetName()
#                print "%6s : %16s %6.3f -> %16s %6.3f" % (atom.GetName(), atom.GetType(), atom.GetPartialCharge(), mutated_atomtypes_2[atomname], mutated_charges_2[atomname])                
                
            
#            print "%4s %6s        %4s %6s   %6s" % ('atom', 'charge', 'atom', 'charge', 'charge')
#            index = 0
#            for mp in match.GetAtoms():
#                atom1 = mp.pattern
#                atom2 = mp.target                        
#                print "%4s %6.3f <----> %4s %6.3f : %6.3f" % (atom1.GetName(), atom1.GetPartialCharge(), atom2.GetName(), atom2.GetPartialCharge(), common_charge[index])
#                index += 1
#            print
            matchcount += 1
            maxmatch = max(maxmatch, match.NumAtoms())

        print '%(ligand1_name)12s : %(ligand2_name)12s    %(maxmatch)8d atoms' % vars()

