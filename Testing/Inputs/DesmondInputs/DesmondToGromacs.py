import sys
import string
import os
path = os.getcwd() + "/"
sys.path.append('../..')
from System import System
System._sys = System.System("Redone Sample")
print "System initialized\n"
f = open( "DesmondToGromacsErrors.txt", "w")
for dir in os.listdir(path):
        if os.path.isdir("%s" % (dir)):
                for file in os.listdir("%s" %(dir)):
                        if ".cms" in file:
                                filename = string.rstrip(file,".cms")
                                print "\nReading in Desmond structure %s"%(filename)
                                sys.path.append('../../ParserFiles')
                                from DesmondParser import DesmondParser
                                DesmondParser = DesmondParser()
                                try:
					DesmondParser.readFile(path+filename+"/"+filename+".cms")
				except Exception,e:
					f.write("\nError reading %s -- %s" %(filename, e))
				filename_out = filename + "_OUT"
				try:
					if not (os.path.exists('%sDesmondToGromacsInputs/%s' %(path,filename))):
						os.makedirs('%sDesmondToGromacsInputs/%s' %(path,filename))
					filename_out = os.path.abspath('%sDesmondToGromacsInputs/%s/%s' %(path,filename, filename_out))
					print "\nWriting out Gromacs topology %s.top"%(filename_out)
					from GromacsParser import GromacsTopologyParser
					GromacsTopologyParser = GromacsTopologyParser()
					GromacsTopologyParser.writeTopology(filename_out+".top")
				except Exception, e:
					f.write("\nError writing %s.top -- %s" %(filename_out,e))
				print "\nWriting in Gromacs Stucture as %s"%(filename_out+".gro")
				from GromacsParser import GromacsStructureParser
				try:
					GromacsStructureParser.writeStructure(filename_out+".gro")
				except Exception,e:
					f.write("\nError  writing %s.gro -- %s" (filename_out, e))
f.close()
