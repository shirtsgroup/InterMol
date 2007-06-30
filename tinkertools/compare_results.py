# compare_results.py 
#
# Compare free energy differences for all AMOEBA Tinker free energy calculations.
#
# John D. Chodera
# 21 Jun 2007

import mmtools.tinkertools.tinkertools as tinkertools
import mmtools.utilities.Units as Units # for units
import mmtools.utilities.Constants as Constants # for constants

from numarray import *
from math import *

import AlchemicalTools

import os
import os.path

# PARAMETERS

# Complete list.
#molecules = [ 'methane', 'butane', 'methanol', 'ethanol', 'toluene', 'p-cresol', 'methyl sulfide', 'acetamide',  'N-methylacetamide', 'N,N-dimethylformamide', 'acetic acid', 'formaldehyde', 'acetaldehyde', 'hydrogen sulfide', 'dimethyl sulfide', 'dimethyl disulfide', 'methyl ethyl sulfide', 'benzene', 'ethylbenzene', '1H-imidazole', 'ammonia', 'dimethylamine', 'methylamine', 'pyrrolidine', 'trimethylamine', 'ethylamine']

# partial list
molecules = [ 'methane', 'butane', 'methanol', 'ethanol', 'toluene', 'p-cresol', 'methyl sulfide', 'acetamide',  'N-methylacetamide', 'N,N-dimethylformamide', 'acetic acid', 'formaldehyde', 'acetaldehyde', 'benzene', 'ethylbenzene', '1H-imidazole', 'ammonia', 'dimethylamine', 'methylamine', 'pyrrolidine', 'trimethylamine', 'ethylamine']

dataset_names = ['amoeba', 'exptvals', 'ambervals']
datasets = dict()
import pickle
for dataset_name in dataset_names:
    filename = dataset_name + '.pickle'
    infile = open(filename, 'rb')
    datasets[dataset_name] = pickle.load(infile)
    infile.close()

# get a list of molecules
molecules = datasets['amoeba'].keys()
nmolecules = len(molecules)

# exclude amines because we screwed up and for got to put in the lone pairs 
tmp = []
for molecule in molecules:
    if molecule.find('amine') < 0 and (molecule != 'pyrrolidine'):
        tmp.append(molecule)
molecules = tmp    

# report results
cum_DDG2 = 0.0
cum_d2DDG = 0.0
print "%24s %12s   %12s   %12s" % ('molecule', 'AMOEBA', 'experiment', 'difference')
datafile = open('amoeba-experiment.data','w')
errorfile = open('amoeba-experiment.errors','w')
for molecule in molecules:

    #DeltaG_amoeba = datasets['amoeba'][molecule]['total']['DeltaG'] - datasets['amoeba'][molecule]['dispersion-correction']
    DeltaG_amoeba = datasets['amoeba'][molecule]['total']['DeltaG'] 
    dDeltaG_amoeba = datasets['amoeba'][molecule]['total']['dDeltaG']    

    DeltaG_exp = datasets['exptvals'][molecule]
    dDeltaG_exp = 0.2

    DDG = (DeltaG_amoeba - DeltaG_exp)
    dDDG = sqrt(dDeltaG_amoeba**2 + dDeltaG_exp**2)

    # Report total free energy contribution.
    print "%24s %5.1f +- %3.1f   %5.1f +- %3.1f : %5.1f +- %3.1f kcal/mol" % (molecule, DeltaG_amoeba, dDeltaG_amoeba, DeltaG_exp, dDeltaG_exp, DDG, dDDG)    

    cum_DDG2 += DDG**2 
    cum_d2DDG += dDDG**2 * DDG**2

    datafile.write("%16.8f %16.8f %16.8f %16.8f\n" % (DeltaG_exp, DeltaG_amoeba, dDeltaG_exp, dDeltaG_amoeba))
    errorfile.write("%16.8f %16.8f\n" % (dDeltaG_exp, dDeltaG_amoeba))


datafile.close()
errorfile.close()

print ""
DDG = sqrt(cum_DDG2 / nmolecules)
dDDG = sqrt(cum_d2DDG / (nmolecules * cum_DDG2))
print "%24s %6.2f +- %4.2f" % ('RMS', DDG, dDDG)

# Compute mean size of polarization correction
correction_sum = 0.0
for molecule in molecules:
    correction_sum += datasets['amoeba'][molecule]['dispersion-correction']
correction_sum /= nmolecules
print "average correction: %f kcal/mol" % correction_sum
