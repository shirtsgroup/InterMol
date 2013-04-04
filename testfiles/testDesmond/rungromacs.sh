#!/usr/bin/bash

# get the inport files
export PATH=/usr/local/gromacs/bin/:$PATH
echo $PATH
export GRO=../../3.gro
export TOP=../../3.top
export MDP=../../grompp.mdp
export TESTNAME=3
# run grompp
grompp -f $MDP -c $GRO -p $TOP -o $TESTNAME.tpr

# run mdp
mdrun -deffnm $TESTNAME -maxwarn 10

# run genergy
echo 1 2 3 4 5 6 7 8 9 10 11 12 13 0 | g_energy -f $TESTNAME -o $TESTNAME.xvg   
