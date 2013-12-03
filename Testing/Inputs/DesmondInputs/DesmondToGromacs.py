import sys
import string
import os
path = '/home/ahy3nz/mmtools/ctools/Testing/Inputs/DesmondInputs/'
sys.path.append('../..')
from System import System
System._sys = System.System("Redone Sample")
print "System initialized\n"
for dir in os.listdir('/home/ahy3nz/mmtools/ctools/Testing/Inputs/DesmondInputs/'):
        if os.path.isdir("%s" % (dir)):
                for file in os.listdir("%s" %(dir)):
                        if ".cms" in file:
                                filename = string.rstrip(file,".cms")
                                print "\nReading in Desmond structure %s"%(filename)
                                sys.path.append('../../ParserFiles')
                                from DesmondParser import DesmondParser
                                DesmondParser = DesmondParser()
                                DesmondParser.readFile(path+filename+"/"+filename+".cms")
				filename_out = filename + "_OUT"
				if not (os.path.exists('/home/ahy3nz/mmtools/ctools/Testing/Inputs/DesmondInputs/DesmondToGromacsInputs/%s' %(filename))):
					os.makedirs('/home/ahy3nz/mmtools/ctools/Testing/Inputs/DesmondInputs/DesmondToGromacsInputs/%s' %(filename))
				filename_out = os.path.abspath('/home/ahy3nz/mmtools/ctools/Testing/Inputs/DesmondInputs/DesmondToGromacsInputs/%s/%s' %(filename, filename_out))
				print "\nWriting out Gromacs topology %s.top "%(filename_out)
				from GromacsParser import GromacsTopologyParser
				GromacsTopologyParser = GromacsTopologyParser()
				GromacsTopologyParser.writeTopology(filename_out+".top")
				print "\nWriting in Gromacs Stucture as %s"%(filename_out+".gro")
				from GromacsParser import GromacsStructureParser
				GromacsStructureParser.writeStructure(filename_out+".gro")

