#!/usr/bin/python

"""Add restraints to a topology file."""
"""Reads an anglegps.ndx index file to identify appropriate atoms for the restraints."""

import os
import re
import shutil
from optparse import OptionParser
parser=OptionParser()

def tail1(filename):
    """Returns the last line from the file named in the argument"""
    f=open(filename,"rt")
    line=""
    for l in f:
        line = l
    return line

def get_parms_from_gangle(prefix):
    """Gets the restraint parameters from the g_angle output files for the system, named
       [prefix]_avg_{phi_A,phi_B,phi_C,theta_A,theta_B,r_aA}.xvg
       Assumes only one line of data in the xvg file"""
    retvals = []
    for i in ["phi_A","phi_B","phi_C","theta_A","theta_B"]:
        rv = tail1(prefix+"_avg_"+i+".xvg").split()[1]
        retvals.append(rv)
    # special-case the distance because the naming convention is different
    rv = tail1(prefix+"_r_aA.xvg").split()[1]
    retvals.append(rv)
    return retvals

def dihedral_restraints(atomlists,phis,force_constants):
	""" Return a list corresponding to the dihedral restraint section """
	dih_type = 1
	dih_power = 2
	dphiA = 0
	kfacA = "0.00"
	dphiB = 0
	restr_sec = []
	restr_sec.append("[ dihedral_restraints ]\n")
	restr_sec.append(";  i    j    k    l type label power phiA dphiA  kfacA  phiB  dphiB   kfacB\n")
	restr_sec.append(";mixed state\n")
	for i in range(0,len(atomlists)):
		atoms = atomlists[i]
		phiA = phis[i]
		fcB = force_constants[i]
		restr_sec.append("%(atoms)s %(dih_type)4d %(i)4d %(dih_power)4d %(phiA)s %(dphiA)4d %(kfacA)s %(phiA)s %(dphiB)4d %(fcB)s\n"%vars());
	restr_sec.append("; zero state\n")
	for i in range(0,len(atomlists)):
		atoms = atomlists[i]
		phiA = phis[i]
		restr_sec.append(";%(atoms)s %(dih_type)4d %(i)4d %(dih_power)4d %(phiA)s %(dphiA)4d %(kfacA)s\n"%vars());
	restr_sec.append("; one state\n")
	for i in range(0,len(atomlists)):
		atoms = atomlists[i]
		phiA = phis[i]
		fcB = force_constants[i]
		restr_sec.append(";%(atoms)s %(dih_type)4d %(i)4d %(dih_power)4d %(phiA)s %(dphiB)4d %(fcB)s\n"%vars());
	restr_sec.append("\n")
	return restr_sec

def angle_restraints(atomlists,thetas,force_constants):
	""" Return a list corresponding to the dihedral restraint section """
	ang_type = 1
	mult = 1
	fcA = "0.00"
	restr_sec = []
	restr_sec.append("[ angle_restraints ]\n")
	restr_sec.append(";  i    j    k    l type theta0A fcA mult  theta0B    fcB  mult\n");
	restr_sec.append(";mixed state -- need to have MULT listed twice!!!\n")
	for i in range(0,len(atomlists)):
		atoms = atomlists[i]
		theta0A = thetas[i]
		fcB = force_constants[i]
		restr_sec.append("%(atoms)s %(ang_type)4d %(theta0A)s %(fcA)s %(mult)4d %(theta0A)s %(fcB)s %(mult)4d\n"%vars());
	
	restr_sec.append("; zero state\n")
	for i in range(0,len(atomlists)):
		atoms = atomlists[i]
		theta0A = thetas[i]
		fcB = force_constants[i]
		restr_sec.append(";%(atoms)s %(ang_type)4d %(theta0A)s %(fcA)s %(mult)4d\n"%vars());
		
	restr_sec.append("; one state\n")
	for i in range(0,len(atomlists)):
		atoms = atomlists[i]
		theta0A = thetas[i]
		fcB = force_constants[i]
		restr_sec.append(";%(atoms)s %(ang_type)4d %(theta0A)s %(fcB)s %(mult)4d\n"%vars());
	restr_sec.append("\n")
	return restr_sec

def distance_restraints(atomlists,rs,force_constants):
	""" Return a list corresponding to the dihedral restraint section """
	dis_type = 1
	dis_typeprime = 1
	r2 = "10.0"
	fcA = "0.0"
	restr_sec = []
	restr_sec.append("[ distance_restraints ]\n")
	restr_sec.append(";  i    j type label typeprime    r0A   r1A   r2A    fcA  r0B   r1B   r2B   fcB\n");
	restr_sec.append(";mixed state\n")
	for i in range(0,len(atomlists)):
		atoms = atomlists[i]
		r = rs[i]
		fcB = force_constants[i]
		restr_sec.append("%(atoms)s %(dis_type)4d %(i)4d %(dis_typeprime)4d %(r)s %(r)s %(r2)s %(fcA)s %(r)s %(r)s %(r2)s %(fcB)s\n"%vars());
	
	restr_sec.append("; zero state\n")
	for i in range(0,len(atomlists)):
		atoms = atomlists[i]
		r = rs[i]
		fcB = force_constants[i]
		restr_sec.append(";%(atoms)s %(dis_type)4d %(i)4d %(dis_typeprime)4d %(r)s %(r)s %(r2)s %(fcA)s\n"%vars());
		
	restr_sec.append("; one state\n")
	for i in range(0,len(atomlists)):
		atoms = atomlists[i]
		r = rs[i]
		fcB = force_constants[i]
		restr_sec.append(";%(atoms)s %(dis_type)4d %(i)4d %(dis_typeprime)4d %(r)s %(r)s %(r2)s %(fcB)s\n"%vars());
	restr_sec.append("\n")
	
	return restr_sec

####
#Read input
#####
parser.add_option("-f", "--file", dest="infile", help="Input angle groups file. Default: 'anglegps.ndx' in current directory.", default='anglegps.ndx', metavar="FILE")
parser.add_option("-n", "--name", dest="topname", help="Name used for topology files in which to put restraints. This will modifify the (name)_nochg and (name)_novdw files in the setup directory, and additionally create a (name)_restr.top file there from the restrained topology. Required.", metavar="FILE")
parser.add_option("-p", "--path", dest="path", help="Path to topology files to modify. Default: '../setup'.", default='../setup', metavar="PATH")

(options,args)=parser.parse_args()
infile=options.infile
pathprefix=options.path
#Error check
if not (os.path.isfile(infile)):
   parser.error("Cannot open input file %(infile)s." % vars())
if not options.topname:
   parser.error("Please enter a topology to modify using the -n option.")
tempname=os.path.join(pathprefix,options.topname+'.top')
if not (os.path.isfile(tempname)):
   parser.error("Cannot open topology files, i.e.  %(tempname)s." % vars())
topname=options.topname

#Create new topology file for restraining
shutil.copy(os.path.join(pathprefix,topname+'.top'),os.path.join(pathprefix,topname+'_restr.top'))
shutil.copy(os.path.join(pathprefix,topname+'.g96'),os.path.join(pathprefix,topname+'_restr.g96'))
#Set list of topology files to which to add restraints later.
topfiles=[os.path.join(pathprefix,topname+'_restr.top')]
#For convenience, back up these topologies as the existing names with _norst after them
# ih - this is already a copy, wtf?
#for file in topfiles:
#   shutil.copy(file,file[:-4]+'_norst.top')


#Get user input for preferred angles
# ih - get these from g_angle files
print get_parms_from_gangle(topname)
(phiA,phiB,phiC,thetaA,thetaB,raA) = get_parms_from_gangle(topname)


#Read atom numbers from input file
file=open(infile,'r')
inputlines=file.readlines()

#Currently assumes sections are in order: first atoms for distance, then thetas, then phis.
#Strip bracketed sections
atomlines=[]
p=re.compile(r'\[\s*\w+.*\]')
for line in inputlines:
   m=p.match(line)
   if not m:
     atomlines.append(line[:-1]) #Strip newlines


#Now modify the two angle lines to repeat the  atom numbers.
p=re.compile(r'(?P<a>\d+)\s(?P<b>\d+)\s(?P<c>\d+)')

anglenum=[2,3]
for num in anglenum:
  m=p.match(atomlines[num])
  #Replace line with repitition of one of the atoms; we want the angle from ab_cb (that's what g_angle computes).
  atomlines[num]=m.group('a')+' '+m.group('b')+' '+m.group('c')+' '+m.group('b')


#Error checking
FileError="Cannot open specified file."

#Modify topology files to add restraints
rmax=float(raA)+0.2
for top in topfiles:
   #Read topology file.
   if not os.path.isfile(top):
     print "Cannot open topology file %(top)s." % vars()
     raise FileError
   else:
     print "Adding restraints to topology file %(top)s." % vars() 
   file=open(top,'r')
   toplines=file.readlines()
   file.close()

   #Write topology file.
   file=open(top,'w')
   
   #Set up text to write:
   # hardcode the force constants
   # taken from dmobley's free_energy_setup example
   # fc_dihed = "fc_dihed"
   fc_dihed =   " 41.84   " 
   # fc_angle= "fc_angle"
   fc_angle =  "  41.84   "
   # fc_dist = "fc_dist"
   fc_dist  =  " 4184.0  "
   restrsec=[]
   
   #restrsec.append('[ dihedral_restraints ]\n')
   #restrsec.append('; i j k l            type     label  phi  dphi  kfac          power\n')
   # different defn of restraints for the free-energy code
   #restrsec.append(';  i    j    k    l type label power phiA dphiA  kfacA  phiB  dphiB   kfacB\n')
   #restrsec.append('%s   1        0      %s  0  %s     2\n' % (atomlines[4],phiA,fc_dihed))
   #restrsec.append('%s   1        1      %s  0  %s     2\n' % (atomlines[5],phiB,fc_dihed))
   #restrsec.append('%s   1        2      %s  0  %s     2\n\n' % (atomlines[6],phiC,fc_dihed))
   #restrsec.append('[ angle_restraints ]\n; i j k l            type    theta0     fc             mult\n')
   #restrsec.append('%s  1      %s       %s     1\n' % (atomlines[2],thetaA,fc_angle))
   #restrsec.append('%s  1      %s       %s     1\n\n' % (atomlines[3],thetaB,fc_angle))
   #restrsec.append('[ distance_restraints ]\n; i j       type    label typeprime  r0    r1   r2     fc\n')
   #restrsec.append('%s %s   1    0   1  %s  %s  %s   %s\n\n' % (atomlines[0],atomlines[1],raA,raA,rmax,fc_dist))
   
   dih_rest = dihedral_restraints([atomlines[4],atomlines[5],atomlines[6]],[phiA,phiB,phiC],[fc_dihed,fc_dihed,fc_dihed])
   ang_rest = angle_restraints([atomlines[2],atomlines[3]],[thetaA,thetaB],[fc_angle,fc_angle])
   dis_rest = distance_restraints(["%s %s"%(atomlines[0],atomlines[1])],[raA],[fc_dist])
   map(restrsec.append,dih_rest)
   map(restrsec.append,ang_rest)
   map(restrsec.append,dis_rest)

   #Restraints section needs to immediately precede the include of the water topology, so read and write
   #lines sequentially, inserting the restraints when we get to the appropriate spot.
   # ih - that looks like a dirty lie, seems to me like it needs to be before the top file switches to another molecule
   #      for now assume the protein/ligand combo is the first molecule in the file
   p=re.compile(r'^\[ moleculetype \]')
   restraints_inserted = False
   rst=re.compile(r'.*dihedral_restraints.*')
   line=0
   linenum=len(toplines)
   # skip the first [ moleculetype ] line
   while line<linenum:
           file.write(toplines[line])
	   line+=1
	   if p.match(toplines[line-1]):
		   break

   while line<linenum:
        #Check for old restraints section
        rstrmatch=rst.match(toplines[line])
        if rstrmatch:
           print "Found old restraint section; removing it to use new restraints."
           m=p.match(toplines[line])
           while not m:
               line+=1
               m=p.match(toplines[line])
        m=p.match(toplines[line])
        if m and not restraints_inserted:
           print "Found match at line %d; writing restraint section." % line
           file.writelines(restrsec)
           file.write(toplines[line])
	   restraints_inserted = True
        else:
           file.write(toplines[line])
        line+=1
           

   #Done writing, close file.
   file.close() 

#Done looping over files.

#Compute and output restraint energy.
springconst=10.0
kB=1.381*6.02214/4184.0
T=300
raA=float(raA)*10.0  #convert to angstroms
from math import pi
from math import sqrt
from math import sin
from math import log
beta=1.0/(kB*T)
theta_A=float(thetaA)*pi/180.0
theta_B=float(thetaB)*pi/180.0
numer=8*(pi**2)*1660*sqrt(springconst**6)
denom=(raA**2.0)*sin(theta_A)*sin(theta_B)*((2.0*pi*kB*T)**3.0)
G_r=kB*T*log(numer/denom)
print "Energy of restraining from standard state is %(G_r)s." % vars()

file=open('result.txt','w')
file.write('%.2f' % round(G_r,2))
file.close()

