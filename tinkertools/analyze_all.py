# analyze_all.py 
#
# Compute free energy differences for all AMOEBA Tinker free energy calculations.
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

temperature = 298.0 * Units.K
beta = 1.0 / (Constants.kB * temperature)
print "(1/beta) = %f kcal/mol" % ((1./beta) / (Units.kcal/Units.mol))

debug = False

# Convert molecle names to directory names (removing commas and spaces).
directories = []
for molecule in molecules:
    molecule = molecule.replace(' ','')
    molecule = molecule.replace(',','')
    molecule += '0'
    directories.append(molecule)

# SUBROUTINES
def estimate_free_energy_difference(basedir):
    """Estimate the cumulative free energy difference and uncertainty for a leg of the thermodynamic cycle.

    ARGUMENTS
      basedir (string) - the directory containing subdirectories for each lambda value

    RETURNS
      DeltaF, dDeltaF (float) - the free energy difference and uncertainty in going from lambda = 0 -> lambda = 1
    """

    # Make an ordered list of all lambda values (which appear as subdirectory names within 'basedir').
    import commands
    lambda_values = commands.getoutput('ls %s' % basedir).split() # strings of sorted lambda values
    nlambda = len(lambda_values) # number of lambda values

    # Determine forward and reverse work values (as applicable) for each lambda value.
    work = dict()
    for n in range(0,nlambda):
        lambda_value = lambda_values[n]
        if debug: print "\nlambda = %s" % lambda_value

        # Create dictionary for storing forward and reverse work values.
        work[lambda_value] = dict()
        
        # Construct directory name.
        subdir = os.path.join(basedir, lambda_value)

        # List all reprocessed values in each directory.
        reprocessed_lambda_values = commands.getoutput('ls %s' % os.path.join(subdir, 'repr*.log')).split()        

        # Read potential energies for this lambda value
        repr_filename = os.path.join(subdir, 'repr%s.log' % lambda_value)
        U_t = tinkertools.read_potential_energy(repr_filename)

        # Read U_{n+1}
        if n < nlambda-1:
            repr_filename = os.path.join(subdir, 'repr%s.log' % lambda_values[n+1])
            if debug: print "Reading %s..." % repr_filename
            U_t_next = tinkertools.read_potential_energy(repr_filename)

            # Compute correlated work values.
            common_indices = range(0, min(len(U_t), len(U_t_next)))
            DeltaU = U_t_next[common_indices] - U_t[common_indices]
            
            # Determine uncorrelated work measurements.
            g = AlchemicalTools.statisticalInefficiency(DeltaU, DeltaU)
            uncorrelated_indices = range(0, len(common_indices), int(ceil(g)))
            if debug: print "There are %d uncorrelated forward work measurements." % len(uncorrelated_indices)                        

            # Store work values
            work[lambda_value]['forward'] = DeltaU[uncorrelated_indices]
        
        # Read U_{n-1}
        if n > 0:
            repr_filename = os.path.join(subdir, 'repr%s.log' % lambda_values[n-1])
            if debug: print "Reading %s..." % repr_filename
            U_t_prev = tinkertools.read_potential_energy(repr_filename)

            # Compute correlated work values.
            common_indices = range(0, min(len(U_t), len(U_t_prev)))
            DeltaU = U_t_prev[common_indices] - U_t[common_indices]
            
            # Determine uncorrelated work measurements.
            g = AlchemicalTools.statisticalInefficiency(DeltaU, DeltaU)
            uncorrelated_indices = range(0, len(common_indices), int(ceil(g)))
            if debug: print "There are %d uncorrelated reverse work measurements." % len(uncorrelated_indices)                        

            # Store work values
            work[lambda_value]['reverse'] = DeltaU[uncorrelated_indices]

    # Compute integrated free energy difference.
    DeltaF = 0.0
    d2DeltaF = 0.0
    if debug: print "\nCumulative free energy difference:"
    if debug: print "%20s   %8s %8s   %8s %8s" % ('change', 'DeltaF', 'dDeltaF', 'DeltaF', 'dDeltaF')
    for n in range(0,nlambda-1):
        # Determine lambda values
        lambda_n = lambda_values[n]
        lambda_n_plus_one = lambda_values[n+1]
        
        # Extract forward and reverse work
        w_F = work[lambda_n]['forward']
        w_R = work[lambda_n_plus_one]['reverse']

        # Compute free energy difference using BAR.
        (DeltaF_n, dDeltaF_n) = AlchemicalTools.BAR(beta / (1.0/(Units.kcal/Units.mol)), w_F, w_R, compute_uncertainty = True)
        
        # Accumulate free energy and statistical uncertainty.
        DeltaF += DeltaF_n
        d2DeltaF += dDeltaF_n**2        

        # Report
        if debug: print "%8s -> %8s   %8.3f %8.3f   %8.2f %8.3f" % (lambda_n, lambda_n_plus_one, DeltaF_n, dDeltaF_n, DeltaF, sqrt(d2DeltaF))
    if debug: print "TOTAL"


    # Compute std dev of error.
    dDeltaF = sqrt(d2DeltaF)
    
    # Return free energy difference and uncertainty.
    return (DeltaF, dDeltaF)

def compute_transfer_free_energy(basedir, return_components = False, debug = False):
    """Compute gas-to-water transfer free energy.

    REQUIRED ARGUMENTS
      basedir (string) - the base directory for the Tinker AMOEBA free energy calculation, containing a subdirectory for each leg of the cycle

    OPTIONAL ARGUMENTS
      return_components (boolean) - if True, will return a dictionary of components instead

    RETURNS
      [DeltaF, dDeltaF] (float) - the free energy of transfer in kcal/mol
    or
      components (dictionary)
    """
    
    # Define the legs of the thermodynamic cycle, where there is a separate subdirectory for each.
    cycle = ['nochg_vac', 'novdw_wat', 'nochg_wat']
    # Set the signs for each leg of the thermodynamic cycle contributing to the overall cycle of interest
    sign = dict()
    sign['nochg_vac'] = 1.0
    sign['novdw_wat'] = -1.0
    sign['nochg_wat'] = -1.0

    # Store components in dict.
    components = dict()
    
    # Accumulate free energy across each leg of the cycle.
    DeltaF = 0.0
    d2DeltaF = 0.0
    for leg in cycle:
        if debug: print leg
        
        # Compute alchemical contribution from this leg of the cycle
        import os.path
        (DeltaF_leg, dDeltaF_leg) = estimate_free_energy_difference(os.path.join(basedir, leg))
        DeltaF_leg *= sign[leg]

        # Compute pV correction to this leg of the cycle if in solvent.
        if leg != 'nochg_vac':
            pV_correction = tinkertools.pV_correction(os.path.join(basedir, leg), pressure = 1.0 * Units.atm) / (Units.kcal/Units.mol) * sign[leg]
            DeltaF_leg += pV_correction

        # Store components in a dict.
        components[leg] = dict()
        components[leg]['DeltaG'] = DeltaF_leg
        components[leg]['dDeltaG'] = dDeltaF_leg

        # Accumulate free energy difference and uncertainty.
        DeltaF += DeltaF_leg
        d2DeltaF += dDeltaF_leg**2
                
    # Compute uncertainty from squared uncertainty.
    dDeltaF = sqrt(d2DeltaF)
        
    # Compute NVT LR correction using the fully-interacting endpoint.
    import os.path
    LR_correction = tinkertools.LR_correction(os.path.join(basedir,'nochg_wat','0.0')) / (Units.kcal/Units.mol)
    if debug: print "LR correction = %f kcal/mol" % LR_correction
    # Accumulate contribution
    DeltaF += LR_correction
    # Store component
    components['dispersion-correction'] = dict()
    components['dispersion-correction'] = LR_correction

    # Store totals.
    components['total'] = dict()
    components['total']['DeltaG'] = DeltaF
    components['total']['dDeltaG'] = dDeltaF

    if return_components:
        return components
    else:
        return (DeltaF, dDeltaF)

# store transfer free energies in a dictionary
transfer_free_energies = dict()

nmolecules = len(molecules)
for n in range(0, nmolecules):
    molecule = molecules[n]
    basedir = directories[n]

    # Compute free energy of transfer from gas to water.
    if debug: print "Processing directory %s..." % basedir
    components = compute_transfer_free_energy(basedir, return_components = True)

    DeltaG = components['total']['DeltaG']
    dDeltaG = components['total']['dDeltaG']    

    # Store values.
    transfer_free_energies[molecule] = components
    
    # Report total free energy contribution.
    print "%32s %7.3f +- %5.3f kcal/mol" % (molecule, DeltaG, dDeltaG)    

# Write pickle of computed values
import pickle
pickle_file = open('amoeba.pickle', 'wb')
pickle.dump(transfer_free_energies, pickle_file)
pickle_file.close()
