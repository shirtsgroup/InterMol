# analyze.py 
#
# Example analysis script for AMOEBA free energy calculations in Tinker.
#
# John D. Chodera
# 21 Jun 2007

import mmtools.tinkertools.tinkertools as tinkertools
import mmtools.utilities.Units as Units # for units
import mmtools.utilities.Constants as Constants # for constants

from numarray import *
from math import *

import AlchemicalTools

# PARAMETERS
basedir = 'methanol0' # base directory for free energy calculation to analyze
temperature = 298.0 * Units.K
beta = 1.0 / (Constants.kB * temperature)
print "(1/beta) = %f kcal/mol" % ((1./beta) / (Units.kcal/Units.mol))

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
        print "\nlambda = %s" % lambda_value

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
            print "Reading %s..." % repr_filename
            U_t_next = tinkertools.read_potential_energy(repr_filename)

            # Compute correlated work values.
            common_indices = range(0, min(len(U_t), len(U_t_next)))
            DeltaU = U_t_next[common_indices] - U_t[common_indices]
            
            # Determine uncorrelated work measurements.
            g = AlchemicalTools.statisticalInefficiency(DeltaU, DeltaU)
            uncorrelated_indices = range(0, len(common_indices), int(ceil(g)))
            print "There are %d uncorrelated forward work measurements." % len(uncorrelated_indices)                        

            # Store work values
            work[lambda_value]['forward'] = DeltaU[uncorrelated_indices]
        
        # Read U_{n-1}
        if n > 0:
            repr_filename = os.path.join(subdir, 'repr%s.log' % lambda_values[n-1])
            print "Reading %s..." % repr_filename
            U_t_prev = tinkertools.read_potential_energy(repr_filename)

            # Compute correlated work values.
            common_indices = range(0, min(len(U_t), len(U_t_prev)))
            DeltaU = U_t_prev[common_indices] - U_t[common_indices]
            
            # Determine uncorrelated work measurements.
            g = AlchemicalTools.statisticalInefficiency(DeltaU, DeltaU)
            uncorrelated_indices = range(0, len(common_indices), int(ceil(g)))
            print "There are %d uncorrelated reverse work measurements." % len(uncorrelated_indices)                        

            # Store work values
            work[lambda_value]['reverse'] = DeltaU[uncorrelated_indices]

    # Compute integrated free energy difference.
    DeltaF = 0.0
    d2DeltaF = 0.0
    print "\nCumulative free energy difference:"
    print "%20s   %8s %8s   %8s %8s" % ('change', 'DeltaF', 'dDeltaF', 'DeltaF', 'dDeltaF')
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
        print "%8s -> %8s   %8.3f %8.3f   %8.2f %8.3f" % (lambda_n, lambda_n_plus_one, DeltaF_n, dDeltaF_n, DeltaF, sqrt(d2DeltaF))
    print "TOTAL"


    # Compute std dev of error.
    dDeltaF = sqrt(d2DeltaF)
    
    # Return free energy difference and uncertainty.
    return (DeltaF, dDeltaF)

# MAIN

# Define the legs of the thermodynamic cycle, where there is a separate subdirectory for each.
cycle = ['nochg_vac', 'novdw_wat', 'nochg_wat']
# Set the signs for each leg of the thermodynamic cycle contributing to the overall cycle of interest
sign = dict()
sign['nochg_vac'] = 1.0
sign['novdw_wat'] = -1.0
sign['nochg_wat'] = -1.0

# Accumulate free energy across each leg of the cycle.
DeltaF = 0.0
d2DeltaF = 0.0
for leg in cycle:
    print leg

    # Compute alchemical contribution from this leg of the cycle
    (DeltaF_leg, dDeltaF_leg) = estimate_free_energy_difference(os.path.join(basedir, leg))
    DeltaF_leg *= sign[leg]
    # Accumulate free energy difference and uncertainty.
    DeltaF += DeltaF_leg
    d2DeltaF += dDeltaF_leg**2

    if leg != 'nochg_vac':
        # Compute pV correction to this leg of the cycle
        pV_correction = tinkertools.pV_correction(os.path.join(basedir, leg), pressure = 1.0 * Units.atm)
        pV_correction *= sign[leg]
        print "pV correction: %e kcal/mol" % (pV_correction / (Units.kcal/Units.mol))
        DeltaF += pV_correction / (Units.kcal/Units.mol)

# Compute uncertainty from squared uncertainty.
dDeltaF = sqrt(d2DeltaF)

# Compute NVT LR correction using the fully-interacting endpoint.
import os.path
LR_correction = tinkertools.LR_correction(os.path.join(basedir,'nochg_wat','0.0'), cutoff = 9.0 * Units.A)
dDeltaF += LR_correction
print "LR correction = %f kcal/mol" % (LR_correction / (Units.kcal/Units.mol))

# Report total free energy contribution.
print "DeltaF = %.3f +- %.3f kcal/mol" % (DeltaF, dDeltaF)

