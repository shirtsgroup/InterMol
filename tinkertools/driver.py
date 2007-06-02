#!/usr/bin/python

from tinkertools import *

#Example driver script

#Get non-water atom types
mol_atomtypes=get_nonwater_types('test/water.xyz')
print mol_atomtypes

print "Generating parameter file..."
#Generate templates with full parameters for regular and vacuum calcs
generate_parameter_file('TEMPLATE.key', '/dmobley1/tinker/newhydration/amoeba.prm', mol_atomtypes, 'test/TEMPLATE_new.key')
generate_parameter_file('TEMPLATE_vac.key', '/dmobley1/tinker/newhydration/amoeba.prm', mol_atomtypes, 'test/TEMPLATE_new_vac.key')

print "Setting up calculation..."
#Default choices for simulation length and lambda values are wired into this but can be altered.
setup_calc('test/TEMPLATE_new.key', 'test/TEMPLATE_new_vac.key', 'test', mol_atomtypes, '/dmobley1/tinker/newhydration/amoeba.prm', 'test/water.xyz', 'test/water.dyn')
