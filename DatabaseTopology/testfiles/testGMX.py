import Topology.Driver as Driver
import clearDB
import pdb

#welcome message
print "Welcome to MMTools v0.1\n\nWe will now read in the default topology (system_GMX.top) and structure files (system_GMX.gro) and output to \"system_GMX_out.gro\" and \"system_GMX_out.top\""

username = raw_input("Enter your username (hint: cs4750pmc8pausr or cs4750pmc8palgn): ")
password = raw_input("Enter your password (hint: usr or lgn):")

print "\n\n" + username + "\t" + password

#prompt for user name with hint
#prompt for password with hint



Driver.initSystem("Solvated GMX", "stardock.cs.virginia.edu", username,password)
#read only
#Driver.initSystem("Solvated GMX", "stardock.cs.virginia.edu", "lgn",false)

Driver.loadTopology("system_GMX.top")
#Driver.loadStructure("system_GMX.gro")
Driver.writeTopology("system_GMX.array")
#Driver.writeStructure("system_GMX_out.gro")
#Driver.writeTopology("system_GMX_out.top")


