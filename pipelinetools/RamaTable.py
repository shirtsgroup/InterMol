#=============================================================================================
# RamaTable.py
#
# Written by Vincent Voelz, Stanford University, 2007-07-06
                     
#=============================================================================================
# IMPORTS

import sys, os
import random, math, string, tempfile 
import numpy

from optparse import OptionParser    # For parsing of command line arguments

from mmtools.utilities import Units

# We need ZAM functions to modify dihedral angles, etc. 
sys.path.append('/Users/vincentvoelz/scripts/mmtools/zam')
import protein, proteinfunc, proteinconst, pdbtools, coords



#=============================================================================================
# CLASSES
#=============================================================================================

class RamaTable:
    """An object for storing and querying Ramachandran Dihedral Angle distributions"""
 
    def __init__(self, title='Ramachandran Table'):
        """Initialize the RamaTable object"""

        self.title = title
        
        self.basins = [1,2,4,5,6]

        self.RamaProb = proteinconst.RamaProb
        self.RamaNPhi = self.RamaProb.shape[0]
        self.RamaNPsi = self.RamaProb.shape[1]
        self.RamaDPhi = 360. / self.RamaNPhi
        self.RamaDPsi = 360. / self.RamaNPsi
        self.PhiBins = numpy.arange(-180.,180., self.RamaDPhi)
        self.PsiBins = numpy.arange(-180.,180., self.RamaDPsi)

        # initialize the (normalized) tables for each structural basin 
        self.RamaBasin1Prob = None
        self.RamaBasin2Prob = None 
        self.RamaBasin4Prob = None 
        self.RamaBasin5Prob = None 
        self.RamaBasin6Prob = None 
        self.buildBasins() 

    def info(self):
        """Return text with info about this RamaTable object."""

        txt = 'Information for %s:\n'%self.title
        txt = txt + 'self.RamaProb = \n'
        txt = txt + self.table_repr(self.RamaProb)
        txt = txt + 'self.RamaBasin1Prob = \n'
        txt = txt + self.table_repr(self.RamaBasin1Prob)
        txt = txt + 'self.RamaBasin2Prob = \n'
        txt = txt + self.table_repr(self.RamaBasin2Prob)
        txt = txt + 'self.RamaBasin4Prob = \n'
        txt = txt + self.table_repr(self.RamaBasin4Prob)
        txt = txt + 'self.RamaBasin5Prob = \n'
        txt = txt + self.table_repr(self.RamaBasin5Prob)
        txt = txt + 'self.RamaBasin6Prob = \n'
        txt = txt + self.table_repr(self.RamaBasin6Prob)
        return txt 

    def table_repr(self, table, formatstr='%-4.3f'):
        """Return a space-delimited repr string for a table array."""

        txt = ''
        for i in range(0,table.shape[0]):
            for j in range(0,table.shape[1]):
                txt = txt + formatstr%table[i,j]+' '
            txt = txt + '\n'
        return txt


    def basinFromDihedral(self, phi, psi):
        """Given a pair of phi, psi dihedral angles, return the corresponding Jha/Sosnick basin (1,2,4,5, or 6).
        basin=1           Beta: (-140 +/- 40, -205 +/- 105)
        basin=2           PP2:  (-50 +/- 50, -205 +/- 105)
        basin=4           Helix-R: (-90 +/- 90, -25 +/- 75)
        basin=5           Helix-L (90 +/- 90, 25 +/- 75)
        basin=6           Gamma: (90 +/- 90, -155 +/- 105)
        """

        # transpose phi and psi back to the unit circle [-180., 180.)
        if phi < -180.:  phi=phi+360.
        if phi >= 180.:  phi=phi-360.
        if psi < -180.:  psi=psi+360.       
        if psi >= 180.:  psi=psi-360.

        # basin 1: Beta: (-140 +/- 40, -205 +/- 105)
        # two possible straddles
        if (phi >= -180.) & (phi < -100.):
          if (psi >= -180.) & (psi < -100.):
              return 1
          if (psi >= 50.) & (psi < 180.):
              return 1

        # basin 2: PP2:  (-50 +/- 50, -205 +/- 105)
        # two possible straddles
        if (phi >= -100.) & (phi < 0.):
          if (psi >= -180.) & (psi < -100.):
              return 2
          if (psi >= 50.) & (psi < 180.):
              return 2

        # basin 4: Helix-R: (-90 +/- 90, -25 +/- 75)
        if (phi >= -180.) & (phi < 0.):
          if (psi >= -100.) & (psi < 50.):
              return 4

        # basin 5: Helix-L: (90 +/- 90, 25 +/- 75)
        if (phi >= 0.) & (phi < 180.):
          if (psi >= -50.) & (psi < 100.):
              return 4

        # basin 6: Gamma: (90 +/- 90, -155 +/- 105)
        # gamma is "everything else"
        return 6

        
    def buildBasins(self):
        """
        basin=1           Beta: (-140 +/- 40, -205 +/- 105)
        basin=2           PP2:  (-50 +/- 50, -205 +/- 105)
        basin=4           Helix-R: (-90 +/- 90, -25 +/- 75)
        basin=5           Helix-L (90 +/- 90, 25 +/- 75)
        basin=6           Gamma: (90 +/- 90, -155 +/- 105)
        """

        self.RamaBasin1Prob = self.buildBasinRama(-180., -100., -310., -100.)
        self.RamaBasin2Prob = self.buildBasinRama(-100., 0., -310., -100.)
        self.RamaBasin4Prob = self.buildBasinRama(-180., 0., -100., 50.)
        self.RamaBasin5Prob = self.buildBasinRama( 0., 180., -50., 100.)
        self.RamaBasin6Prob = self.buildBasinRama( 0., 180., -260., -50.)


    def buildBasinRama(self, minphi, maxphi, minpsi, maxpsi):
        """Build a subset RamaTable from self.RamaProb"""
   
        # catch straddling phi regions
        phi_regions = []
        if (minphi < -180.0) & (maxphi >= -180.0):
            phi_regions.append( (-180.0, maxphi) )
            phi_regions.append( (minphi+360.0, 180.0) )
        elif (minphi <= 180.0) & (maxphi > 180.0):
            phi_regions.append( (minphi, 180.0) )
            phi_regions.append( (-180.0, maxphi-360.0) )
        else:
            phi_regions.append( (minphi, maxphi) )
 
        # catch straddling psi regions
        psi_regions = []
        if (minpsi < -180.0) & (maxpsi >= -180.0):
            psi_regions.append( (-180.0, maxpsi) )
            psi_regions.append( (minpsi+360.0, 180.0) )
        elif (minpsi <= 180.0) & (maxpsi > 180.0):
            psi_regions.append( (minpsi, 180.0) )
            psi_regions.append( (-180.0, maxpsi-360.0) )
        else:
            psi_regions.append( (minpsi, maxpsi) )

        RamaProb = numpy.zeros( self.RamaProb.shape )
        # for each rectangular block of phi-psi, copy the probabilities into the temp table
        for phi_region in phi_regions:
          for psi_region in psi_regions:

              print 'Filling phi_region',phi_region
              print 'Filling psi_region',psi_region

              # walk through each phi,psi bin and assign add the appropriate probability
              for i in range(0,self.RamaNPhi):
                for j in range(0, self.RamaNPsi):
                 
                  factor = 1.0
                  # PHI factor is the fraction of probabiily that comes from a edge-split bin 
                  if (self.PhiBins[i] <= phi_region[0]) & (self.PhiBins[i]+self.RamaDPhi > phi_region[0]):
                      factor = factor * (self.PhiBins[i]+self.RamaDPhi - phi_region[0])/self.RamaDPhi  
                  elif (self.PhiBins[i] > phi_region[0]) & (self.PhiBins[i]+self.RamaDPhi <= phi_region[1]):
                      factor = factor * 1.
                  elif (self.PhiBins[i] <= phi_region[1]) & (self.PhiBins[i]+self.RamaDPhi > phi_region[1]):
                      factor = factor * (phi_region[1]-self.PhiBins[i])/self.RamaDPhi
                  else:
                      factor = 0.0
                  
                  # PSI factor is the fraction of probabiily that comes from a edge-split bin
                  if factor > 0.0:
                    if (self.PsiBins[j] <= psi_region[0]) & (self.PsiBins[j]+self.RamaDPsi > psi_region[0]):
                      factor = factor * (self.PsiBins[j]+self.RamaDPsi - psi_region[0])/self.RamaDPsi
                    elif (self.PsiBins[j] > psi_region[0]) & (self.PsiBins[j]+self.RamaDPsi <= psi_region[1]):
                      factor = factor * 1.
                    elif (self.PsiBins[j] <= psi_region[1]) & (self.PsiBins[j]+self.RamaDPsi > psi_region[1]):
                      factor = factor * (psi_region[1]-self.PsiBins[j])/self.RamaDPsi
                    else:
                      factor = 0.0

                  RamaProb[i,j] = RamaProb[i,j] + self.RamaProb[i,j]*factor

        return (RamaProb + 1e-8) / RamaProb.sum()    # make sure cells are non-zero (so log prob is defined)


    def getRandomBasin(self, basin=None):
        """Returns a random basin from [1,2,4,5,6], each with equal probability."""
        
        return self.basins[ random.randint(0,len(self.basins)-1) ]

    
    def getRandomDihedral(self, basin=None):
        """Returns a random (phi, psi) angle and its probability."

    OPTIONAL
        basin             to specify a subset of Jha/Sosnick basins to draw from

            basin = 1, 2, 4, 5 or 6 or...
            basin = 'Beta', 'PP2', 'Helix-R', '', ''

        basin=1           Beta: (-140 +/- 40, -205 +/- 105)
        basin=2           PP2:  (-50 +/- 50, -205 +/- 105)
        basin=4           Helix-R: (-90 +/- 90, -25 +/- 75)
        basin=5           Helix-L (90 +/- 90, 25 +/- 75)
        basin=6           Gamma: (90 +/- 90, -155 +/- 105)
        """

        if basin==1:  RamaProb = self.RamaBasin1Prob
        elif basin==2:  RamaProb = self.RamaBasin2Prob
        elif basin==4:  RamaProb = self.RamaBasin4Prob
        elif basin==5:  RamaProb = self.RamaBasin5Prob
        elif basin==6:  RamaProb = self.RamaBasin6Prob
        else:
            RamaProb = self.RamaProb

        # return a random phi,psi bin, with probability
        r = random.random()
        CumProb = 0.
        for i in range(self.RamaNPhi):
          for j in range(self.RamaNPsi):
            CumProb += RamaProb[i,j]
            if r < CumProb:
              Phi = -180. + self.RamaDPhi * (i + random.random()) 
              Psi = -180. + self.RamaDPsi * (j + random.random()) 
              return Phi, Psi, RamaProb[i,j]
            
    def GetRamaProb(Phi, Psi, basin=None):
        """Returns the probability associated with a Phi, Psi angle."""
        if basin==1:  RamaProb = self.RamaBasin1Prob
        elif basin==2:  RamaProb = self.RamaBasin2Prob
        elif basin==4:  RamaProb = self.RamaBasin4Prob
        elif basin==5:  RamaProb = self.RamaBasin5Prob
        elif basin==6:  RamaProb = self.RamaBasin6Prob
        else:
            RamaProb = self.RamaProb

        PhiInd = int(clip((Phi + 180.) / RamaDPhi, 0, RamaNPhi - 1))
        PsiInd = int(clip((Psi + 180.) / RamaDPsi, 0, RamaNPsi - 1))
        return RamaProb[PhiInd, PsiInd]
 
