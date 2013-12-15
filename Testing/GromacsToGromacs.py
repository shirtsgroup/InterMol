import sys
import os.path
from System import System

System._sys = System("Redone Sample")
print "System initialized\n"
print "ENTER FILENAME?\n"
filename = raw_input("")
while not (os.path.isfile(filename+".gro") or os.path.isfile(filename+".top")):
 print "File doesn't exit, try again\n"
 filename = raw_input("")
filename_out = filename + "_OUT"
from GromacsExt.GromacsTopologyParser import GromacsTopologyParser
print 'Reading in Gromacs topology "%s"...' %(filename)
if not GromacsTopologyParser._GroTopParser:
  GromacsTopologyParser._GroTopParser = GromacsTopologyParser()
GromacsTopologyParser._GroTopParser.parseTopology(filename + '.top')
import GromacsExt.GromacsStructureParser as GromacsStructureParser
print 'Reading in Gromacs structure "%s"...' %(filename)
GromacsStructureParser.readStructure(filename + '.gro')
print "\nWriting out Gromacs topology %s"%(filename_out+".top")
import GromacsExt.GromacsTopologyParser as GromacsTopologyParser
GromacsTopologyParser = GromacsTopologyParser()
GromacsTopologyParser.writeTopology(filename_out+".top")
print "\nWriting in Gromacs structure %s"%(filename_out+".gro")
import GromacsExt.GromacsStructureParser as GromacsStructureParser
GromacsStructureParser.writeStructure(filename_out+".gro")

