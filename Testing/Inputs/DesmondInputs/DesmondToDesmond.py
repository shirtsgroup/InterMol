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
				if not (os.path.exists('/home/ahy3nz/mmtools/ctools/Testing/Inputs/DesmondInputs/DesmondToDesmondInputs/%s' %(filename))):
					os.makedirs('/home/ahy3nz/mmtools/ctools/Testing/Inputs/DesmondInputs/DesmondToDesmondInputs/%s' %(filename))
				filename_out = os.path.abspath('/home/ahy3nz/mmtools/ctools/Testing/Inputs/DesmondInputs/DesmondToDesmondInputs/%s/%s' %(filename, filename_out))
				print "\nWriting out Desmond structure %s"%(filename_out)
				try:
					DesmondParser.writeFile(filename_out+".cms")
				except Exception as inst:
					print"\nError reading %s" %s(filename)
