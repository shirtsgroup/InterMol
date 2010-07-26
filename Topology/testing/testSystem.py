#!/usr/bin/env python

import os, sys, glob
sys.path.append("/home/ctk3b")
import simtk.unit as units
from mmtools.Topology.System import *


print '##########################################'
print '# Tests of the Topology.System methods #'
print '##########################################'


print "Create a System object from file ...",
s = System()
print "Done."

print "Set the name of the System object to Frank ..." 
s.setName('Frank')
print 'New name:', t.getName()
print "Done."

print "Add a molecule to the System object ..."
s.addMolecule()
print "There are now", s.getNumMolecules(), "molelcules in the Topology."
print "Done."

print "Insert a molecule at position 1 ..."
s.insertMolecule(1)
print "There are now", s.getNumMolecules(), "molecules in the Topology."
print "Done."

print "Delete molecule 0 from the Topology...",
s.delMolecule(0) 
print "Done."



print '#########################################################'
print '# Tests of the Topology.Topology.Topology methods #'
print '#########################################################'

topology = s.molecules[0]

print "Add a particle."
mass = 12.0 * units.amu
topology.addParticle(mass)
print 'Done.'

print "Add a NonbondedForce."
nonbondedForce = NonbondedForce()
topology.addForce(nonbondedForce)
print "Done."

print "Create a deep copy."
import copy
system_copy = copy.deepcopy(topology)
print 'system_copy', system_copy
print 'Done.'

