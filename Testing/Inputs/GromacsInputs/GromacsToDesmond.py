import sys
import os.path
from ctools.System import System

System._sys = System("Redone Sample")
print "System initialized\n"
print "ENTER FILENAME?\n"
filename = raw_input("")
while not  (os.path.isfile(filename+".gro") or os.path.isFile(filename+".top")):
 print "File doesn't exist, try again\n"
 filename = raw_input("")
filename_out = filename + "_OUT"
from ctools.GromacsExt.GromacsTopologyParser import GromacsTopologyParser
print 'Reading in Gromacs topology "%s"...' %(filename)
if not GromacsTopologyParser._GroTopParser:
  GromacsTopologyParser._GroTopParser = GromacsTopologyParser()
GromacsTopologyParser._GroTopParser.parseTopology(filename + '.top')
import ctools.GromacsExt.GromacsStructureParser as GromacsStructureParser
print 'Reading in Gromacs structure "%s"...' %(filename)
GromacsStructureParser.readStructure(filename + '.gro')
print "Writing out Desmond structure %s"%(filename_out)
import ctools.DesmondExt.DesmondParser as DesmondParser
DesmondParser = DesmondParser()
DesmondParser.writeFile(filename_out)

