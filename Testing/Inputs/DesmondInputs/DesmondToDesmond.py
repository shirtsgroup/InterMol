import sys
sys.path.append('../..')
from System import System
System._sys = System.System("Redone Sample")
print "System initialized\n"
print "ENTER FILENAME?\n"
filename = raw_input("")
#while not (os.path.isfile(filename+".cms")):
# print "File doesn't exist, try again\n"
# filename = raw_intput("")
#filename_out = filename + "_OUT"
print "\nReading in Desmond structure %s"%(filename)
sys.path.append('../../ParserFiles')
from DesmondParser import DesmondParser
DesmondParser = DesmondParser()
DesmondParser.readFile(filename+".cms")
print "\nWriting out Desmond structure %s"%(filename_out)
DesmondParser.writeFile(filename_out+".cms")

