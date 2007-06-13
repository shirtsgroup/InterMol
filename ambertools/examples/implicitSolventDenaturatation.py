# imlicitSolventDenaturation.py 
#
# Thermally denature a PDB file in GB/SA implicit solvent with the AMBER tools.
#
# John D. Chodera
# 11 Jun 2007

import mmtools.ambertools.hacks as amberhacks # import unpleasant hacks
import mmtools.utilities.Units as Units # for units

# PARAMETERS
initial_pdb_filename = "1ENH_enhd_mcce_out.pdb" # PDB file of inital AMBER-ready protein
workdir = "1ENH-melt" # name of directory for AMBER files

# MAIN

# Get the space-separated three-letter-code sequence from the PDB file.
print "Extracting sequence from %s" % initial_pdb_filename
sequence = amberhacks.getPresentSequence(initial_pdb_filename, allow_four_letters = True)
# Change LYP to LYS
sequence = sequence.replace("LYP","LYS")
print "sequence = "
print sequence

# Create work directory.
import os
if os.path.exists(workdir)==False:
    os.mkdir(workdir)

import mmtools.ambertools.ImplicitSolventSimulation as amber

# Set up the system.
print "Creating system..."
sim = amber.ImplicitSolventSimulation(sequence = sequence, leaprc_filename = "leaprc.ff03", gbradii = "amber6", workdir = workdir, debug = True)

# Minimize the system.
print "minimizing..."
minimized_pdb_filename = sim.minimize(useropts = { 'maxcyc' : 100 } )
#sim.minimize()

# Run 1 ps of high-temperature dynamics.
print "running dynamics..."
md_pdb_filename = sim.dynamics(temperature = 600.0 * Units.K, duration = 1.0 * Units.ps)

# Resume dynamics for an additional 1 ps
print "resuming dynamics..."
md_pdb_filename = sim.dynamics(temperature = 600.0 * Units.K, duration = 1.0 * Units.ps)




