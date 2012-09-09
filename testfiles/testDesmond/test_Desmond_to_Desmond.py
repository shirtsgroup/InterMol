import sys
from ctools.System import System

System._sys = System("Redone Sample")
print "System initialized\n"

#READING IN DESMOND--WRITING OUT IN DESMOND
#filename = "Redone.cms"
filename = "3.cms"
print "Reading in Desmond structure %s"%(filename)
import ctools.DesmondExt.DesmondParser as DesmondParser
DesmondParser = DesmondParser()
DesmondParser.readFile(filename)

#filename = "Redone_out.cms"
filename = "3_out.cms"
print "Writing out Desmond structure %s"%(filename)
DesmondParser.writeFile(filename)

