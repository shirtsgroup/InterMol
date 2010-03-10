#!/usr/bin/env python

from mmtools.gromacstools.devel.SimulationPreparation.PrepWorkerExplicitSolvation import *
from mmtools.gromacstools.devel.SimulationPreparation.SystemSetup import *


# NOTES
#
# In this example, we will prepare an explicit solvent simulation of the lambda repressor protein
# The work will be done in a temporary directory, with status shown in standard output
# The output will be saved to the directory 'output'.

# Initialize the system
ss = SystemSetup()

# Set simulation system parameters
ss.setForcefield('ffamber03')
ss.setSaltConditions('NaCl', 0.100)  # Use 100 mM sodium chloride concentration for
ss.set_boxType('cubic')              # Use a rectangular (cubic) perodic box type
ss.set_boxSoluteDistance(1.2)        # Size the box so that the boundaries are 12 Angstroms away from the lambda

# Initialize the PrepWorker object
pw  = PrepWorkerExplicitSolvation('lambda.gro', finalOutputDir='lambda-ff03-0.1Msalt', finalOutputName='lambda')

# prepare the simulation
pw.prepare(ss)

