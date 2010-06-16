#!/usr/bin/env python

import os, sys, glob
sys.path.append("/home/ctk3b")
import simtk.unit as units
from mmtools.Topology.Topology import *


print '##########################################'
print '# Tests of the Topology.Topology methods #'
print '##########################################'


print "Create a Topology object from file ...",
t = Topology()
print "Done."

print "Set the name of the Topology object to Frank ..." 
t.setName('Frank')
print 'New name:', t.getName()
print "Done."

print "Add a molecule to the Topology object ..."
t.addMolecule()
print "There are now", t.getNumMolecules(), "molelcules in the Topology."
print "Done."

print "Insert a molecule at position 1 ..."
t.insertMolecule(1)
print "There are now", t.getNumMolecules(), "molecules in the Topology."
print "Done."

print "Delete molecule 0 from the Topology...",
t.delMolecule(0) 
print "Done."



print '#########################################################'
print '# Tests of the Topology.Topology.TopologySystem methods #'
print '#########################################################'

system = t.molecules[0]

print "Add a particle."
mass = 12.0 * units.amu
system.addParticle(mass)
print 'Done.'

print "Add a NonbondedForce."
nonbondedForce = NonbondedForce()
system.addForce(nonbondedForce)
print "Done."

print "Create a deep copy."
import copy
system_copy = copy.deepcopy(system)
print 'system_copy', system_copy
print 'Done.'

