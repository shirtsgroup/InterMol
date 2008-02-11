# EnergyHistogram.py
#
# Last updated:  Dec 6, 2007
# 12/06/2007 - Added reading from flat one-column energy files


import sys, os, math, numpy

### Definitions ###



def acceptanceRatioTempSwap(hist1, hist2):
    """Return the average acceptance ratio between transitions from T1->T2, as
    < min(1, exp(-D)) >, where D = (beta2 - beta1)*E1 - (g2 - g1)"""

    # map the energies on the metropolis values min(1, exp(-D))
    metropolis = numpy.zeros( hist1.bins.shape )
#    metropolis = numpy.zeros( hist1.values.shape )
    for i in range(0, hist1.bins.shape[0]):
#    for i in range(0, hist1.values.shape[0]):
#       D = (hist2.beta - hist1.beta)*hist1.values[i] - (hist2.g - hist1.g)
       D = (hist2.beta - hist1.beta)*(hist1.bins[i]+hist1.binwidth/2) - (hist2.g - hist1.g)
       # print 'D[%d]'%i,D
       if D < 0:
           metropolis[i] = 1
       else:
           metropolis[i] = min(1., numpy.exp(-1.*D))

    # print 'metropolis', metropolis
    accept = hist1.expectation(metropolis) 
#    accept = sum(metropolis)/metropolis.shape[0]
    return accept 

### Classes ###
	
class Histogram(object):
	
    def __init__(self, nbins=100):
        """Initialize the EnergyHistogram object."""
		
	self.nbins = nbins
	self.values = []	    # holds a list of the energy values
	self.counts = None      # holds the histogram counts
	self.ncounts = None      # holds the normalized histogram counts
	self.bins = None        # holds the lower edges of the bins
	self.binwidth = None    # holds the width of the bin

    def read_values(self, values):
	"""Read in energies as a list."""

	self.values = numpy.array(values)
	[self.counts, self.bins] = numpy.histogram(values, bins=self.nbins)
	self.ncounts = self.counts/float(sum(self.counts))
	if self.nbins > 1:
            self.binwidth = self.bins[1]-self.bins[0]

    def read_valuesFromFile(self, filename, col=0, comment_chars=['#','@']):
	"""Read in energies from a one-column (or the first column of a) flat text file."""

	fin = open(filename,'r')
	lines = fin.readlines()
        fin.close()
	
	values = []
	for line in lines:
	  if comment_chars.count( line[0] ) == 0:
	    fields = line.strip().split()
	    if len(fields) > 0:
		values.append(float(fields[col]))
              	    
	self.read_values(values)
	

    def rebin(self, minval, maxval, nbins=100):
        """Re-bins the histogram according a new range."""
        self.nbins = nbins
        [self.counts, self.bins] = numpy.histogram(self.values, bins=self.nbins, range=(minval,maxval) )
        self.ncounts = self.counts/float(sum(self.counts))
        if self.nbins > 1:
            self.binwidth = self.bins[1]-self.bins[0]
 

    def expectation(self, observable):
        """Returns the expectation values of the numpy.array observable.  This
        must have the same shape (same number of bins) as self.values"""

        return numpy.dot(observable,self.ncounts)


    def mean(self, usehist=False):
        """Returns the mean of the values.  The usehist=True will compute
        the mean from the histogram bins, which in theory is faster for large data sets."""      
        if usehist:
            bincenters = self.bins + self.binwidth/2.
            return numpy.dot(bincenters,self.ncounts)
        else:
            return numpy.mean(self.values)

    def var(self, usehist=False):
        """Returns the mean of the values.  The usehist=True will compute
        the mean from the histogram bins, which in theory is faster for large data sets."""      
        if usehist:
            bincenters = self.bins + self.binwidth/2.
            return numpy.dot(bincenters*bincenters, self.ncounts) - (self.mean(usehist=usehist))**2
        else:
            return numpy.var(self.values)

    def write_dat(self, outfile):
        """Write a two-column (tab-delimited) bins, counts to dat file."""

        fout = open(outfile,'w')
        fout.write('E\tcounts\n')
        for i in range(0,len(self.bins)):
	        fout.write('%f\t%f\n'%(self.bins[i],self.ncounts[i]))
        fout.close()
        return


			
    def printASCII(self, nrows=20, normed=False):
        """Prints out an ASCII normalized plot of the histogram."""
	
	shading=""".,:;=@&%$*##"""
        	
        if normed:
            maxcount = float(max(self.ncounts))
            counts = self.ncounts
        else:
            maxcount = float(max(self.counts))
            counts = self.counts

        vspacing = maxcount/float(nrows)   # the spacing of each row in the printout
        print 'vspacing =',vspacing

        outstr = ''
	for j in numpy.arange(maxcount, 0., -1*vspacing):
	    # label the vertical axis
            outstr += '%-8s |'%('%-2.2f'%abs(j))
            for i in range(0, len(self.bins)):
	        if ((j - counts[i]) > 0):
                    if (j - counts[i]) < vspacing:
                        shading_index = int(10*(1.0 - (j-counts[i])/vspacing)) 
                        outstr += shading[shading_index]
                    else:
                        outstr += ' '
		else:
		    outstr += '#'
	    outstr += '\n'

        # print lower axis
        sofar = '%-8s *'%' ' + '%-8s'%('%f'%self.bins[0] )
        outstr += sofar
        for i in range(len(sofar)-10,len(self.bins)):
            outstr += '-'
        outstr += '%-8f'%self.bins[-1]
        print outstr
				
		
		
class EnergyHistogram(Histogram):
    """A histogram object, with special features for energy histograms."""

    def __init__(self, nbins=100, useRosetta_kB=False):
        """Initialize the class, with extras for EnergyHistograms."""

        # The energy files have energies defined in terms of kJ/mol
        self.kB = 8.314472 * 0.001    # in kJ/( mol*K )
        if useRosetta_kB:
            self.kB = 1  # for Rosetta
        self.T = 300.	# in K
        self.beta = 1./(self.kB*self.T)
        self.g = 0.     # unitless weight g_T for simulated tempering  

        self.nbins = nbins
        self.values = []            # holds a list/array of the energy values
        self.wvalues = []           # holds a weighted list of the energy values {E_i - g}
        self.counts = None      # holds the histogram counts
        self.ncounts = None      # holds the normalized histogram counts
        self.bins = None        # holds the lower edges of the bins
        self.binwidth = None    # holds the width of the bin


    def set_temp(self, T):
        """Set the temperature to T (K)."""
        self.T = T
        self.beta = 1./(self.kB*self.T)

    def set_weight(self, g):
        self.g = g

    def weight(self, g=0.):
        """Keep a set of weighted energies {E_i - g} in self.wvalues, and
        recompute the histograms using these energies"""

        self.g = g
        self.wvalues = array(values) - g
        [self.counts, self.bins] = numpy.histogram(self.wvalues, bins=self.nbins)
        self.ncounts = self.counts/float(sum(self.counts))
        if self.nbins > 1:
            self.binwidth = self.bins[1]-self.bins[0]

    def rebin(self, minval, maxval, nbins=100, weighted=True):
        """Re-bins the histogram according a new range."""
        self.nbins = nbins
        if weighted:
            [self.counts, self.bins] = numpy.histogram(self.wvalues, bins=self.nbins, range=(minval,maxval) )
        else:
            [self.counts, self.bins] = numpy.histogram(self.values, bins=self.nbins, range=(minval,maxval) )
        self.ncounts = self.counts/float(sum(self.counts))
        if self.nbins > 1:
            self.binwidth = self.bins[1]-self.bins[0]

     
    def read_xvgfile(self, xvgfile, col_index=1, ignore_nframes=False):
	"""Read in a Gromacs energy.xvg file, which looks like the lines below.
           It's assumed that the Potential energy will be in the second column, but
           you can specify other columns by the col_index parameter.

# This file was created by /home/server.171.64.65.111/programs/gromacs314/x86_64-unknown-linux-gnu/bin/g_energy314
# which is part of G R O M A C S:
# God Rules Over Mankind, Animals, Cosmos and Such
# All this happened at: Sun Oct 14 21:14:12 2007
#
@    title "Gromacs Energies"
@    xaxis  label "Time (ps)"
@    yaxis  label "E (kJ mol\S-1\N)"
@TYPE xy
@ view 0.15, 0.15, 0.75, 0.85
@ legend on
@ legend box on
@ legend loctype view
@ legend 0.78, 0.8
@ legend length 2
@ s0 legend "Potential"
@ s1 legend "Kinetic En."
@ s2 legend "Total Energy"
    0.000000  -166791.453125  48837.871094  -117953.578125
    0.500000  -229651.343750  53067.929688  -176583.406250
    1.000000  -235979.015625  47171.546875  -188807.468750
    1.500000  -237780.250000  48555.308594  -189224.937500

"""
        # read lines in the xvg file
        fin = open(xvgfile,'r')
        lines = fin.readlines()
        fin.close()

        # parse the potential energy column data
        keepers = []
        for line in lines:
            if len(line) > 0:
                if (line[0] != '#') & (line[0] != '@'):
                    keepers.append(line)

        if ignore_nframes:
            for i in range(0, ignore_nframes):
                keepers.pop(0)

        for line in keepers:
            fields = line.strip().split()
            self.values.append(float(fields[col_index]))

        # create the histograms
        self.read_values( self.values )

    def acceptRatio(self, hist2):
        metropolis = numpy.zeros( self.bins.shape )
        for i in range(0, self.bins.shape[0]):
            D = (hist2.beta - self.beta)*(self.bins[i]+self.binwidth/2) - (hist2.g - self.g)
#            print "d " + str(D)
            if D < 0:
              metropolis[i] = 1
            else:
              metropolis[i] = min(1., numpy.exp(-1.*D))
#            print "d " + str(D) + " met " + str(metropolis[i])
#            if i == 40:
#              sys.exit(0)

        accept = self.expectation(metropolis)
        return accept

