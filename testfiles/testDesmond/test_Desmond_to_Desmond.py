import sys
from ctools.System import System

System._sys = System("Redone Sample")
print "System initialized\n"

#READING IN DESMOND--WRITING OUT IN DESMOND
#filename = "Redone.cms"
#filename = "3.cms"
filename = "a2a_dppc-out.cms"
#filename = "desmond_md_job-in.cms"
#filename = "example.cms"
#filename = "lck_Me_Cl_complex_12_1-out.cms"
#filename = "lck_Me_Cl_solvent_12_1-out.cms"
#filename = "simulated_annealing_example.cms"
print "Reading in Desmond structure %s"%(filename)
import ctools.DesmondExt.DesmondParser as DesmondParser
DesmondParser = DesmondParser()
DesmondParser.readFile(filename)

#filename = "Redone_OUT.cms"
#filename = "3_OUT.cms"
filename = "a2a_dppc-out_OUT.cms"
#filename = "desmond_md_job-in_OUT.cms"
#filename = "example_OUT.cms"
#filename = "lck_Me_Cl_complex_12_1-out_OUT.cms"
#filename = "lck_Me_Cl_solvent_12_1-out_OUT.cms"
#filename = "simulated_annealing_example_OUT.cms"
print "Writing out Desmond structure %s"%(filename)
DesmondParser.writeFile(filename)

