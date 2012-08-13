import sys
from ctools.System import System

System._sys = System("Redone Sample")
print "System initialized\n"

#READING IN DESMOND--WRITING OUT IN DESMOND
filename = "-Redone.cms"
print "Reading in Desmond structure %s"%(filename)
import ctools.DesmondExt.DesmondParser as DesmondParser
DesmondParser = DesmondParser()
DesmondParser.readFile(filename)

filename = "-Redone_out.cms"
print "Writing out Desmond structure %s"%(filename)
DesmondParser.writeFile(filename)


#READING IN GROMACS--WRITING OUT IN DESMOND
#filename = 'system_GMX'
#from ctools.GromacsExt.GromacsTopologyParser import GromacsTopologyParser
#print 'Reading in Gromacs topology "%s"...' %(filename)
#if not GromacsTopologyParser._GroTopParser:
#  GromacsTopologyParser._GroTopParser = GromacsTopologyParser()
#GromacsTopologyParser._GroTopParser.parseTopology(filename + '.top')

#import ctools.GromacsExt.GromacsStructureParser as GromacsStructureParser
#print 'Reading in Gromacs structure "%s"...' %(filename)
#GromacsStructureParser.readStructure(filename + '.gro')

#filename = "system_GMX_out.cms"
#print "Writing out Desmond structure %s"%(filename)
#import ctools.DesmondExt.DesmondParser as DesmondParser
#DesmondParser = DesmondParser()
#DesmondParser.writeFile(filename)


#READING IN DESMOND--WRITING OUT IN GROMACS
#filename = "-Redone.cms"
#print "Reading in Desmond structure %s"%(filename)
#import ctools.DesmondExt.DesmondParser as DesmondParser
#DesmondParser = DesmondParser()
#DesmondParser.readFile(filename)

#filenameout = "-Redone"
#print "\nWriting in Gromacs topology %s"%(filenameout+".top")
#import ctools.GromacsExt.GromacsTopologyParser as GromacsTopologyParser
#GromacsTopologyParser = GromacsTopologyParser()
#GromacsTopologyParser.writeTopology(filenameout+".top")

#print "\nWriting in Gromacs structure %s"%(filenameout+".gro")
#import ctools.GromacsExt.GromacsStructureParser as GromacsStructureParser
##GromacsStructureParser = GromacsStructureParser()
#GromacsStructureParser.writeStructure(filenameout+".gro")
