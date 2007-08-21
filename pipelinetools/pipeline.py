# pipeline.py
#
# Tools for running proteins through the pipeline!
# Vincent Voelz
# Last modified: June 28, 2007 


# GLOBAL IMPORTS
import mmtools.modellertools.modelPDB as modelPDB
import mmtools.mccetools.mcce as mcce
from mmtools.gromacstools.System import *
from mmtools import pdbtools

import os
import os.path
import sys


# DEBUG flag: if true then save outputs
DEBUG = True

#-----------------------------------
# FUNCTIONS
#-----------------------------------


def thread_model(pdbTemplate, sequence, outPdbFile, captermini=False):
    """Uses Modeller to thread a sequence onto a template structure.
    If captermini=True, will cap the N- and C- termnini with ACE and NH2, respectively. """
    

    # use either a sequence file or a sequence string
    if os.path.exists(sequence):
        sequenceFilename = sequence
    else:
        tmpdir = tempfile.mkdtemp()
        sequenceFilename = os.path.join(tmpdir,'sequence')
        fseq = open(sequenceFilename,'w')
        fseq.write(sequence)
        fseq.close()

    myModel = modelPDB.ModelPDB()
    myModel.makeModel(pdbTemplate, sequenceFilename, outPdbFile, captermini=captermini)
    
    # if capping, rebuild the PDB with ACE and NME residue names 
    if (captermini):
        pdbtools.rebuildPatchedTermini(outPdbFile)
        
    # cleanup
    if not os.path.exists(sequence):
        os.remove(sequenceFilename)
        os.rmdir(tmpdir)


def protonate(inpdbfile, outpdbfile, pH):
    """Generate a protonated PDB file from a template PDB using MCCE."""

    thisdir = os.path.abspath(os.curdir)
    mcceOut = os.path.abspath(outpdbfile)
    prmFile = '../mccetools/prmfiles/run.prm.quick'
    prmFile = os.path.abspath(prmFile)
    mcce.protonatePDB(inpdbfile, mcceOut, pH, os.environ['MCCE_LOCATION'], cleanup=True, prmfile=prmFile, labeledPDBOnly=False)
    os.chdir(thisdir)


def build_gmx(protein, outdir, forcefield='ffamber99p'):
    """Builds a gromacs project according to the specs in the pipelineProtein object.  Assumes
    the protein.pdbfile is in the right (AMBER) format."""

    baseName = protein.getBasename( protein.pdbfile )
    gromacsOut = baseName + "_final"
    g = System(protein.pdbfile, finalOutputDir=outdir, useff=forcefield)
    g.setup.setSaltConditions(protein.salt, protein.saltconc)
    if protein.boxProtocol == 'small':
        g.setup.set_boxType = 'octahedron'
        g.setup.set_boxSoluteDistance(1.5)   # periodic box margin distance, in nanometers
    elif protein.boxProtocol == 'big':
        g.setup.set_boxType = 'octahedron'
        g.setup.setUseAbsBoxSize(True)
        g.setup.setAbsBoxSize('8.0')   # periodic box absolute size, in nanometers (string)

    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
    g.prepare(outname=gromacsOut, outdir=outdir, verbose=True, cleanup=False, debug=DEBUG, protocol='racecar2', checkForFatalErrors=True)



def shoveit(protein, outdir, forcefield='ffamber99p', captermini=False):
    """Shoves a pipelineProtein object through the entire MODELLER -> MCCE --> gromacs pipeline.

    For capping the termini with ACE and NH2, the Modeller program has a nice way of doing
    residue 'patches', so ithis will be used for adding the caps.  WE need to do the following:

    1. Thread the sequence via Modeller, WITHOUT capping --> PDB
    2. Calculate the protonation state via MCCE --> return a non-renamed PDB, and renaming rules
    3. Use Modeller to cap the ends
    4. Apply the renaming rules for protonation states that came from MCCE  

    VV: THIS HASN'T BEEN IMPLEMENTED YET -- I'm working on it!!!!!! 8/21/07
    """

    thisdir = os.curdir

    print 'before:'
    protein.print_info()

    protein.setup(outdir)

    print 'after:'
    protein.print_info()

    # Build a PDB model from the pdbTemplate using MODELLER
    thread_model(protein.pdbfile, protein.seqfile, protein.modelPDBout)

    # run mcce
    if (1):
        mcceOut = os.path.abspath(protein.mccePDBout)
        prmFile = '../../mccetools/prmfiles/run.prm.quick'
        prmFile = os.path.abspath(prmFile)
        mcce.protonatePDB(protein.modelPDBout, mcceOut, protein.pH, os.environ['MCCE_LOCATION'], cleanup=False, prmfile=prmFile, labeledPDBOnly=False)

    # run gromacs setup
    if (1):
        gromacsOut = protein.basename + "_final.pdb"
        # g = system.GromacsSystem(mcceOut, useff=forcefield)    # the old gromacstools way
        g = System(mcceOut, finalOutputDir=outdir, useff=forcefield)
        g.setup.setSaltConditions(protein.salt, protein.saltconc)
        if protein.boxProtocol == 'small':
            g.setup.set_boxType = 'octahedron'
            g.setup.set_boxSoluteDistance(1.5)   # periodic box margin distance, in nanometers
        elif protein.boxProtocol == 'big':
            g.setup.set_boxType = 'octahedron'
            g.setup.setUseAbsBoxSize(True)
            g.setup.setAbsBoxSize('7.0')   # periodic box absolute size, in nanometers (string)

        thisOutDir = os.path.join(thisdir, protein.basename)
        print 'Writing equilibration directory to',thisOutDir,'...'
        if os.path.exists(thisOutDir) == False:
            os.mkdir(thisOutDir)
        os.chdir(thisOutDir)
        g.prepare(verbose=True, cleanup=False, debug=DEBUG, protocol='racecar2', checkForFatalErrors=True)
     

    if(1):
        # cleanup
        if not DEBUG:
            if os.path.exists(protein.modelPDBout):
                os.remove(protein.modelPDBout)
            if os.path.exists(protein.mccePDBout):
                os.remove(protein.mccePDBout)

    os.chdir(thisdir)




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


