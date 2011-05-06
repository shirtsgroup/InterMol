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


Driver.initSystem("Solvated GMX", "stardock.cs.virginia.edu", username, password)
Driver.loadTopology("system2_GMX.top")

Driver.writeTopology("system2_GMX.array")
#pdb.set_trace()
#Driver.loadStructure("system2_GMX.gro")
#Driver.writeStructure("system2_GMX_out.gro")
#Driver.writeTopology("system2_GMX_out.top")

