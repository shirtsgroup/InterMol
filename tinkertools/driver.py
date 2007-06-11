#!/usr/bin/python

from tinkertools import *

#Example driver script

#Get non-water atom types
mol_atomtypes=get_nonwater_types('test/water.xyz')
print mol_atomtypes


print "Setting up calculation..."
#Default choices for simulation length and lambda values are wired into this but can be altered.
setup_calc('TEMPLATE.key', 'TEMPLATE_vac.key', 'test', mol_atomtypes, 'amoeba.prm', 'test/water.xyz', 'test/water.dyn')
