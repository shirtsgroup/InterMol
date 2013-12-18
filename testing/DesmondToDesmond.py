import sys
import string
import os
path = os.getcwd()+"/Inputs/DesmondInputs/"
from System import *

System._sys = System("Redone Sample")
print "System initialized\n"
f = open("DesmondToDesmondErrors.txt", "w")
for dir in os.listdir(path):
	if os.path.isdir("%s" % (dir)):
		for file in os.listdir("%s" %(dir)):
			if (( ".cms" in file)):
				filename = string.rstrip(file,".cms")
				print "\nReading in Desmond structure %s"%(filename)
				from DesmondExt.DesmondParser import DesmondParser
				DesmondParser = DesmondParser()
				try:
					DesmondParser.readFile(path+filename+"/"+filename+".cms")
				except Exception,e:
					print "%s" %(e)
					f.write("\nError reading %s.cms --  %s"  %(filename, e))

				filename_out = filename + "_OUT"
				try:
					 if not (os.path.exists('%sOutputs/DesmondToDesmondOutputs/%s' %(path,filename))):
	                                       	 os.makedirs('%sDesmondToDesmondInputs/%s' %(path, filename))
                               		 filename_out = os.path.abspath('%sDesmondToDesmondOutputs/%s/%s' %(path, filename, filename_out))

					 DesmondParser.writeFile(filename_out+".cms")
				except Exception,e:
					print "%s" %(e)
					f.write("\nError reading %s.cms -- %s" %(filename_out,e))

f.close()
