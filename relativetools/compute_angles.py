#!/usr/bin/python

"""Script to compute angles/distances describing the ligand's six degrees of freedom relative to the protein. Prompts for reference atoms."""

#python computeangles.py -f system.g96 -s minimization-constrained.tpr

import os
import re
import random
import sys
from optparse import OptionParser
parser=OptionParser()


def find_atom_nums(atom_locator_strings,g96file):
    """ Takes a list of strings to recognize in the input file (representing atom locators)
        and returns the GROMACS atom number corresponding to them in a list in the same order.
        ONLY WORKS FOR G96 FILES RIGHT NOW"""
    regexes = [re.compile(i) for i in atom_locator_strings]
    atom_nums = [None for i in atom_locator_strings]
    for line in g96file:
        for i in range(0,len(regexes)):
            if regexes[i].search(line):
                atom_nums[i]=line.split()[3]
    return atom_nums

def three_ligand_atoms(g96file):
    """ Takes a g96 file and returns the atom indices of three atoms to use as restraint atoms.
        Currently uses three random atoms because it's easy to code....
        TODO: something not retarded."""
    regex=re.compile("^ *[0-9]* MOL")
    molatoms=[]
    for line in g96file:
        if regex.search(line):
            molatoms.append(line.split()[3])
    return random.sample(molatoms,3)

def three_ligand_atoms_from_topfile(topfile):
    """ Takes a top file and returns the atom indices of three atoms to use as restraint atoms that are marked as "common" in comments.
        Currently uses three random atoms because it's easy to code....
        TODO: something not retarded."""
    regex=re.compile("MOL.*common")	
    molatoms=[]
    for line in topfile:
        if regex.search(line):
            molatoms.append(line.split()[0])
    return random.sample(molatoms,3)

#######
#Read input
#######
#Input trajectory
parser.add_option("-f","--file", dest="infile", help="FILE= Input trajectory file including extension. Can be xtc or trr. tpr file assumed to have same prefix unless supplied by -s option.", metavar="FILE")
#Input tpr file
parser.add_option("-s", "--tprfile", dest="tprfile", help="Name of input tpr file with extension. Default: Same as name of input trr/xtc file.", metavar="FILE")
#Use old anglegps.ndx?
parser.add_option("-o", "--old", action="store_true", dest="useoldangle", default=False, help="Use old anglegps.ndx file? (Useful if don't want to input atom numbers again if we're just re-analyzing or something). Default: False. No argument required.", metavar="NONE")
# random seed
parser.add_option("-d", "--seed", dest="seed", default=None, help="Random seed to initialize RNG. Defaults to system time", metavar="NONE")
#Name of output
parser.add_option("-n", "--name", dest="outname", help="Name for output files. Default: Same as input name.", metavar="NAME")


(options, args)=parser.parse_args()
if not options.infile:
   parser.error("Please supply input filename using -f.")

FileError="ERROR: Cannot find input file." #For use below.
random.seed(options.seed)

infile=options.infile
#xtc or trr
intype=infile[:-4]
#Error check for existence
if not (os.path.isfile(infile)):
   print "ERROR: Cannot find file %(intype)s."
   raise FileError

#Output name?
outname=''
if not options.outname:
  #Same as input name less extension 
  outname=infile[0:-4]
else:
  outname=options.outname


#Tpr name?
tprname=''
if not options.tprfile:
   tprname=infile[0:-4]+'.tpr'
else:
   tprname=options.tprfile
#Error check
if not (os.path.isfile(tprname)):
   print "ERROR: Cannot find file %(tprname)s."
   raise FileError



#Now, unless using old angle file, create an anglegps.ndx index file with user input
# Use particular atoms for SAMPL calculations
#TODO : generalize this method for arbitrary input residue/atoms
protein_atoms=["155 ASN   CA","106 VAL   CA","109 LEU   CA"]

g96file=open(infile,"rt")
atoms=find_atom_nums(protein_atoms,g96file)
g96file.seek(0)
# lowercase=protein, uppercase=ligand
(atoma,atomb,atomc)=atoms
#ligatoms=three_ligand_atoms(g96file)

topfile=open('system.top','rt')
ligatoms=three_ligand_atoms_from_topfile(topfile) # JDC - use top file
topfile.close()
(atomA,atomB,atomC)=ligatoms

#Write anglegps.ndx
file=open('anglegps.ndx','w')
text="""[ atom_a ]
%(atoma)s
[ atom_A ]
%(atomA)s
[ theta_A ]
%(atomb)s %(atoma)s %(atomA)s
[ theta_B ]
%(atoma)s %(atomA)s %(atomB)s
[ phi_A ]
%(atomc)s %(atomb)s %(atoma)s %(atomA)s
[ phi_B ]
%(atomb)s %(atoma)s %(atomA)s %(atomB)s
[ phi_C ]
%(atoma)s %(atomA)s %(atomB)s %(atomC)s
""" % vars()

file.write(text)
file.close()
#Done with creation of angle groups file.

#Now process angles
anglenames=['theta_A', 'theta_B', 'phi_A', 'phi_B', 'phi_C']
num=1
for anglename in anglenames:
   num+=1
   #Need to decide if angle or dihedral
   #If it has a theta, type is angle
   if anglename.find('theta')>-1:
      type='angle'
   else:
      type='dihedral'
   #Now analyze
   os.system("echo \"%(num)s\" | g_angle -f %(infile)s -s %(tprname)s -n anglegps.ndx -ov %(outname)s_avg_%(anglename)s.xvg -type %(type)s" % vars())

#Done with angles. Do distance, also. Distance atom indices are 0 and 1.
os.system('echo \"0 \n 1\" | g_dist -f %(infile)s -s %(tprname)s -n anglegps.ndx -o %(outname)s_r_aA.xvg' % vars())


