#!/usr/bin/python

from GROMACS import *

gro=GromacsStructure()
gro.load("../sv-M12P/equil0.gro")

seq = gro.atoms.sequence()

out=""
for res in seq:
	res = res.strip()
	if not "SOL" in res: 
		olc = tlc2olc[ res ]
		out += olc

print out
