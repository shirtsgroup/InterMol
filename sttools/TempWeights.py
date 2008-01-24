import os
import os.path
import shutil
import sys

from EnergyHistogram import *


class TempWeightUpdater(object):
    """An object to do the reweighting of different ST temperaure weights
    by equalizing the Acceptance Ratios across neighboring temps. """

    def __init__(self, tempList):
        """Initialize the TempWeightUpdater class with the list of temperatures."""

        self.myT = [float(x) for x in tempList]     # converteds to floats just in case    
        self.numT = len(self.myT)
    
    def getWeightsFromEnergyFiles(self, energyfiles):
        """Reads in a list of energy filenames, one per temperature, each containing
        a single-column list of energy values """
       
        allEnergy = []
        hists = []
        
        for i in range(0, numT):
          if Verbose: print "Reading T " + str(self.numT) + "..."
          allEnergy.append([])
          hists.append(EnergyHistogram())
          hists[i].set_temp(self.myT[i])
          hists[i].read_valuesFromFile(energyfiles[i])
        
        cutoff = 0.001
        bias = 0
        finalWeights = [0]
        for i in range(0, self.numT-1):
          Pij = hists[i].acceptRatio(hists[i+1])
          Pji = hists[i+1].acceptRatio(hists[i])
        
          dP = Pji - Pij - bias
        
          # use Newton's method to get weights
          while abs(dP) > cutoff:
            dP_prime = -Pji-Pij
            hists[i+1].g = hists[i+1].g - dP/dP_prime
        
            Pij = hists[i].acceptRatio(hists[i+1])
            Pji = hists[i+1].acceptRatio(hists[i])
            dP = Pji - Pij - bias
        
          finalWeights.append(hists[i+1].g)
          print "Pij " + str(i) + "->" + str(i+1) + " " + str(Pij)
          print "  Pji " + str(Pji)
        
        print "Final weights " + str(finalWeights)
        return finalWeights
        
    
