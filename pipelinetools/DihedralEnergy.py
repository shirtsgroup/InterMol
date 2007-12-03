#=============================================================================================
# DihedralEnergy.py
#
#=============================================================================================
# AUTHORS
#
# Written by Vincent Voelz, Stanford University, 2007-07-06
#=============================================================================================
                     
#=============================================================================================
# IMPORTS

import sys, os
import random, math, string, tempfile 
import numpy

from mmtools.utilities import Units

# We need ZAM functions to modify dihedral angles, etc. 
sys.path.append('../zam')
import protein, proteinfunc, proteinconst, pdbtools, coords



#=============================================================================================
# CLASSES
#=============================================================================================


class DihedralEnergy:
    """An object for storing uncorrelated and residue-residue (NN) correlated Ramachandran dihedrals,
    with class methods to calculate the dihedral energy function a la Jha, Sosnick, et al. PNAS (2005)."""

    def __init__(self, useLib='native'):
        """Initialize the Dihedral Energy object.

        OPTIONS
        order		'single' - Do not consider nearest-neighbor (NN) effects
                        'double' - Consider NN (nearest-neighbors effects)
                        
        useLib          'coil'
                        'C_A'
                        'C_AB'
                        'C_ABT'
                        'native'
        """

        self.VERBOSE = True
        self.kB = 1.987 * Units.calorie                        # (...per mol per degree)  
        self.kT_ref = self.kB * 300. * Units.K / Units.kcal    # our reference is 300K
        
        self.useLib = useLib

        self.unlikely = 1e-8  # a pseudo-count baseline probability, so that log prob is defined

        self.singleProbs = self.readSingleProbs()   # Example: {'V': {1: 0.3462}}
        if self.VERBOSE:
            print 'self.singleProbs',self.singleProbs 

        self.doubleProbs = self.readDoubleProbs()   # Example: {('V','F'): {(1,4): 0.3462}}
        if self.VERBOSE:
            print 'self.doubleProbs',self.doubleProbs


    def readSingleProbs(self):
        """Read the single-residue basin probabilities from a default RCG *.par file."""

        parfile = os.path.join(os.environ['MMTOOLSPATH'],'pipelinetools/par/Monomers-prob-%s.par'%self.useLib)
        singleProbs = {}

        fin = open(parfile,'r')
        lines = fin.readlines()
        fin.close()

        for line in lines:
            fields = line.split()
            if len(fields) == 3:
                res = fields[0]
                basin = int(fields[1])
                prob = float(fields[2])
                if singleProbs.has_key(res):
                   singleProbs[res][basin] = prob  
                else:
                   singleProbs[res] = {basin:prob}

        return singleProbs


    def readDoubleProbs(self):
        """Read the double-residue basin probabilities from a default RCG *.par file."""

        parfile = os.path.join(os.environ['MMTOOLSPATH'],'pipelinetools/par/Dimers-prob-%s.par'%self.useLib)
        doubleProbs = {}

        fin = open(parfile,'r')
        lines = fin.readlines()
        fin.close()

        for line in lines:
            fields = line.split()
            if len(fields) == 5:
                res1 = fields[0]
                res2 = fields[1]
                basin1 = int(fields[2])
                basin2 = int(fields[3])
                prob = float(fields[4])
                if doubleProbs.has_key( (res1,res2) ):
                   doubleProbs[ (res1,res2) ][ (basin1, basin2) ] = prob
                else:
                   doubleProbs[ (res1,res2) ] = { (basin1, basin2): prob }

        return doubleProbs


    def single(self, res, basin):
        """Calculate the single-residue (no NN-interaction) energy as:

        U = -kT log(P(res, basin))

        ARGUMENTS
            res			single-letter residue code. e.g.  'G' or 'V'
            basin               basin number (int): 1, 2, 4, 5, or 6
        """
        if self.singleProbs.has_key(res):
            if self.singleProbs[res].has_key(basin):
                prob = self.singleProbs[res][basin]
            else:
                prob = self.unlikely
        else:
            prob = self.unlikely     

        return -1.0*self.kT_ref*math.log(prob)

    def double(self, res1, res2, basin1, basin2):
        """Calculate the double-residue (with NN-interaction) energy as

        U = -kT { log(P(res1, basin1)) + log(P(res2, basin2)) +log(P(res1,res2,basin1,basin2)/P(res1, basin1)P(res2, basin2) }
          = -kT { log(P(res1,res2,basin1,basin2) }

        ARGUMENTS
            res1                upstream (closest to N-term) single-letter residue (i) code. e.g.  'G' or 'V'
            res2                downstream single-letter residue (i+1) code. e.g.  'G' or 'V'
            basin1              corresponding basin number (int): 1, 2, 4, 5, or 6
            basin2              corresponding basin number (int): 1, 2, 4, 5, or 6

        """
        if self.doubleProbs.has_key( (res1,res2) ):
            if self.doubleProbs[ (res1,res2) ].has_key( (basin1, basin2) ):
                prob = self.doubleProbs[ (res1, res2) ][ (basin1, basin2) ]
            else:
                prob = self.unlikely
        else:
            prob = self.unlikely

        return -1.0*self.kT_ref*math.log(prob)


 
