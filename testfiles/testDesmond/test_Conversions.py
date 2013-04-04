import sys
import os.path
from ctools.System import System

System._sys = System("Redone Sample")
print "System initialized\n"
num = 0
print "WHICH CONVERSION DO YOU WANT?\n"
print "1. Desmond (Read) --> Desmond (Write)\n"
print "2. Desmond (Read) --> Gromacs (Write)\n"
print "3. Gromacs (Read) --> Desmond (Write)\n"
while True:
  try:
    num = int(raw_input(""))
    if(1 <= num <= 3):
      break
    else:
      print "Not valid response; try again\n"
  except ValueError:
    print "Not a valid response; try again\n"

print "ENTER FILENAME?\n"
filename = raw_input("")
while not (os.path.isfile(filename+".cms") or os.path.isfile(filename+".gro") or os.path.isfile(filename+".top")):
  print "File doesn't exist, try again\n"
  filename = raw_input("")
filename_out = filename + "_OUT"

#READING IN DESMOND--WRITING OUT IN DESMOND
if num == 1:
  print "\nReading in Desmond structure %s"%(filename)
  import ctools.DesmondExt.DesmondParser as DesmondParser
  DesmondParser = DesmondParser()
  DesmondParser.readFile(filename+".cms")
  print "\nWriting out Desmond structure %s"%(filename_out)
  DesmondParser.writeFile(filename_out+".cms")
elif num == 2:
  print "\nReading in Desmond structure %s"%(filename)
  import ctools.DesmondExt.DesmondParser as DesmondParser
  DesmondParser = DesmondParser()
  DesmondParser.readFile(filename+".cms")
  print "\nWriting in Gromacs topology %s"%(filename_out+".top")
  import ctools.GromacsExt.GromacsTopologyParser as GromacsTopologyParser
  GromacsTopologyParser = GromacsTopologyParser()
  GromacsTopologyParser.writeTopology(filename_out+".top")
  print "\nWriting in Gromacs structure %s"%(filename_out+".gro")
  import ctools.GromacsExt.GromacsStructureParser as GromacsStructureParser
  GromacsStructureParser.writeStructure(filename_out+".gro")
else:
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


