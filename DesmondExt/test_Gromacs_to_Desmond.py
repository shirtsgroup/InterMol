import sys
from ctools.System import System

System._sys = System("Redone Sample")
print "System initialized\n"



#READING IN GROMACS--WRITING OUT IN DESMOND
filename = 'system_GMX'
from ctools.GromacsExt.GromacsTopologyParser import GromacsTopologyParser
print 'Reading in Gromacs topology "%s"...' %(filename)
if not GromacsTopologyParser._GroTopParser:
  GromacsTopologyParser._GroTopParser = GromacsTopologyParser()
GromacsTopologyParser._GroTopParser.parseTopology(filename + '.top')

import ctools.GromacsExt.GromacsStructureParser as GromacsStructureParser
print 'Reading in Gromacs structure "%s"...' %(filename)
GromacsStructureParser.readStructure(filename + '.gro')

filename = "system_GMX_out.cms"
print "Writing out Desmond structure %s"%(filename)
import ctools.DesmondExt.DesmondParser as DesmondParser
DesmondParser = DesmondParser()
DesmondParser.writeFile(filename)

