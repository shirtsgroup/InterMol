#!/usr/bin/python

import modelPDB
import sys

# check have right number of inputs
if len(sys.argv) != 3:
        print "usage: ./driver.py pdb_filename sequence_filename"
        sys.exit()

pdb_file = sys.argv[1]
seq_file = sys.argv[2]
outfile = "target.pdb"

print "PDB File: " + pdb_file
print "Sequence File: " + seq_file

myModel = modelPDB.ModelPDB()
myModel.makeModel(pdb_file, seq_file, outfile)


