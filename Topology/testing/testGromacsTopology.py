#!/usr/bin/env python

import os, sys, glob
from mmtools.Topology.GromacsTopology import *


# Tests of the Topology.GromacsTopology tools

print "Create a GromacsTopology object from file ...",
g = GromacsTopology(topfile='lambda.top')
print "Done."

print "Write 'out.top' from the contents of the GromacsTopology object ...",
g.writeTopologyFile('out.top', ExpandIncludes=False, RebuildDirectives=False)
print "Done."

print "Write 'out-expanded.top' with all of the (nested) #include files expanded ...",
g.writeTopologyFile('out-expanded.top', ExpandIncludes=True, RebuildDirectives=False)
print "Done."

print "Show all the directives defined in the GromacsTopologyFileObject:"
all_directives = g.GromacsTopologyFileObject.getAllDirectives([])
for d in all_directives:
   print '\t', d.name.strip()
print "Done."

print "Show only the ParameterDirectives:"
print 'ParameterDirectives:'
for d in g.GromacsTopologyFileObject.ParameterDirectives:
    print '\t',d.name.strip()

print "Show each set of MoleculeDefiniton Directive:"
for mol in range(len(g.GromacsTopologyFileObject.MoleculeDefinitionDirectives)):
  print '\tmol', mol
  for d in g.GromacsTopologyFileObject.MoleculeDefinitionDirectives[mol]:
    print '\t\t', d.name.strip()

print "Show the SystemDirectives:"
for d in g.GromacsTopologyFileObject.SystemDirectives:
    print '\t', d.name.strip()

