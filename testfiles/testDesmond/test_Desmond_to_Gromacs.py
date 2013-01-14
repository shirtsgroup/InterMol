import sys
from ctools.System import System

System._sys = System("Redone Sample")
print "System initialized\n"


#READING IN DESMOND--WRITING OUT IN GROMACS
filename = "3.cms"
print "Reading in Desmond structure %s"%(filename)
import ctools.DesmondExt.DesmondParser as DesmondParser
DesmondParser = DesmondParser()
DesmondParser.readFile(filename)

filenameout = "3"
print "\nWriting in Gromacs topology %s"%(filenameout+".top")
import ctools.GromacsExt.GromacsTopologyParser as GromacsTopologyParser
GromacsTopologyParser = GromacsTopologyParser()
GromacsTopologyParser.writeTopology(filenameout+".top")

print "\nWriting in Gromacs structure %s"%(filenameout+".gro")
import ctools.GromacsExt.GromacsStructureParser as GromacsStructureParser
#GromacsStructureParser = GromacsStructureParser()
GromacsStructureParser.writeStructure(filenameout+".gro")
