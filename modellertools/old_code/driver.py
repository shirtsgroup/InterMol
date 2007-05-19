#!/usr/bin/python

import Model_PDB
import sys

# check have right number of inputs
if len(sys.argv) != 3:
        print "usage: ./driver.py pdb_filename sequence_filename"
        sys.exit()

pdb_file = sys.argv[1]
seq_file = sys.argv[2]

print "PDB File: " + pdb_file
print "Sequence File: " + seq_file

myModel = Model_PDB.Model_PDB()
myModel. make_model(pdb_file, seq_file)


