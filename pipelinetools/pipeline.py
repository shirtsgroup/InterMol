# pipeline.py
#
# Tools for running proteins through the pipeline!
# Vincent Voelz
# Last modified: June 28, 2007 


# GLOBAL IMPORTS
import mmtools.mccetools.mcce as mcce

from mmtools.gromacstools.System import *

import os
import os.path
import sys


# DEBUG flag: if true then save outputs
DEBUG = True


#----------------------------------
# CLASSES
#----------------------------------

class PipelineProtein:

    def __init__(self, pdb=None, seq=None, salt='NaCl', saltconc=0.100, pH=7.0, boxProtocol='small'):
        """Initialize a default data in protein object."""

        self.pdbfile  = pdb
        self.seqfile  = seq
        self.salt     = salt         # Molar concentration of salt (ionic strength)
        self.saltconc = saltconc
        self.pH       = pH
        self.boxProtocol = boxProtocol

        # to be added along the way
        self.modelPDBout = None   # PDB filename to write after the MODELLER phase    
        self.mccePDBout  = None   # PDB filename to write after the MCCE phase
        self.basename = None      # Basename for output files


    def setup(self, outdir):
        """set up filenames and working/output directory."""

        # create outdir if we need to
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        # determine base name
        self.basename = self.getBasename(self.seqfile)
        self.modelPDBout = os.path.join(outdir,(self.basename+"_model.pdb"))
        self.mccePDBout  = os.path.join(outdir,(self.basename+"_mcce.pdb"))

    def getBasename(self, filename):
        """Gets the basename of any file names 'basename.suffix'."""
        extInd = filename.rindex(".")
        return filename[0:extInd]
         
    def print_info(self):
        """Prints info about the PipelineProtein object"""

        print 'pdbfile', self.pdbfile
        print 'seqfile', self.seqfile
        print 'salt', self.salt
        print 'saltconc', self.saltconc 
        print 'pH', self.pH 
        print 'boxProtocol', self.boxProtocol 
        print 'modelPDBout', self.modelPDBout 
        print 'mccePDBout', self.mccePDBout 


