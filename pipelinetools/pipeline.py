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


def thread_model(pdbTemplate, sequence, outPdbFile, chain="_", captermini=False):
    """Uses Modeller to thread a sequence onto a template structure.
    If captermini=True, will cap the N- and C- termnini with ACE and NH2, respectively. """
    
    print 'DEBUG:::'
    print 'pdbTemplate', pdbTemplate
    print 'outPdbFile', outPdbFile 
    print 'os.curdir', os.curdir
    print 'os.listdir', os.listdir('./')

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
    myModel.makeModel(pdbTemplate, sequenceFilename, outPdbFile, chain="_", captermini=captermini)
    
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


def build_gmx(protein, outdir, forcefield='ffamber99p', protocol='racecar2'):
    """Builds a gromacs project according to the specs in the pipelineProtein object.  Assumes
    the protein.pdbfile is in the right (AMBER) format."""


    if os.path.exists(outdir) == False:
        os.mkdir(outdir)

    baseName = protein.getBasename( protein.pdbfile )
    gromacsOut = baseName + "_final"
    g = System(protein.pdbfile, finalOutputDir=outdir, finalOutputName=gromacsOut, useff=forcefield)
    g.setup.setSaltConditions(protein.salt, protein.saltconc)
    if protein.boxProtocol == 'small':
        g.setup.set_boxType = 'octahedron'
        g.setup.set_boxSoluteDistance(1.5)   # periodic box margin distance, in nanometers
    elif protein.boxProtocol == 'big':
        g.setup.set_boxType = 'octahedron'
        g.setup.setUseAbsBoxSize(True)
        g.setup.setAbsBoxSize('8.0')   # periodic box absolute size, in nanometers (string)

    g.prepare(outname=gromacsOut, outdir=outdir, verbose=True, cleanup=False, debug=DEBUG, protocol='racecar2', checkForFatalErrors=True)



def shoveit(protein, outdir, forcefield='ffamber99p', protocol='racecar2', captermini=False, debug=False, verbose=True, cleanup=False, implicitOptions=None, useTable=None, MultiChain=False, SkipMCCE=False):
    """Shoves a pipelineProtein object through the entire MODELLER -> MCCE --> gromacs pipeline.

    For capping the termini with ACE and NH2:
    
    The Modeller program has a nice way of doing
    residue 'patches', so this will be used for adding the caps.  The program will do the following:

    1. Thread the sequence via Modeller, which will patch the termini caps--> PDB
    2. Calculate the protonation state via MCCE (which will know how to handle ACE and NH2 residues because
       of the patches ace.tpl and nh2.tpl in param04 and param08 --> PDB with ffamber-named residues
    3. feed through gromacstools in normal fashion

    OPTIONS

    MultiChain     If set to True, will thread a multichain model.  NOTE: the sequence of the 
                   mutiple chains (in the protein object) must be demarcated by '/'.  Example:
                   AKEEFWVY/AWEKKLELEQVID

    """

    thisdir = '%s'%os.path.abspath(os.curdir)
    protein.setup(outdir)    # Set up modelPDBout filenames and such using the basename of the outdir
    protein.print_info()

    # Build a PDB model from the pdbTemplate using MODELLER
    if MultiChain:
        thread_model(protein.pdbfile, protein.seqfile, protein.modelPDBout, captermini=captermini)
    else:
        thread_model(protein.pdbfile, protein.seqfile, protein.modelPDBout, captermini=captermini)
    

    # run mcce
    mcceOut = os.path.abspath(protein.mccePDBout)
    if SkipMCCE:
        # reformat the MODELLER PDB with the right names
        fin = open(protein.modelPDBout,'r')
        modelPDBlines = fin.readlines()
        fin.close()
        npdb = mcce.rename.nest_pdb(modelPDBlines)
        outlines = mcce.rename.unnest_pdb(mcce.rename.rename_MODELLER_termini(npdb))
        fout = open(mcceOut,'w')
        fout.writelines(outlines)
        fout.close() 
        print 'Writing PDB to', mcceOut
    else:
        prmFile = '../../mccetools/prmfiles/run.prm.quick'
        prmFile = os.path.abspath(prmFile)
        if (captermini):
            mcce.protonatePDB(protein.modelPDBout, mcceOut, protein.pH, os.environ['MCCE_LOCATION'], cleanup=cleanup, prmfile=prmFile, renameTermini=False)
        else:
            mcce.protonatePDB(protein.modelPDBout, mcceOut, protein.pH, os.environ['MCCE_LOCATION'], cleanup=cleanup, prmfile=prmFile)
       
    # run gromacs setup
    if (1):
        gromacsOut = protein.basename + "_final.pdb"
        # g = system.GromacsSystem(mcceOut, useff=forcefield)    # the old gromacstools way
        g = System(mcceOut, finalOutputDir=outdir, finalOutputName=protein.basename, useff=forcefield)
        g.setup.setSaltConditions(protein.salt, protein.saltconc)
        if protein.cosolvent != None:
            g.setup.setCosolventConditions(protein.cosolvent, protein.cosolventconc)
 
        if protein.boxProtocol == 'small':
            g.setup.set_boxType = 'octahedron'
            g.setup.set_boxSoluteDistance(1.5)   # periodic box margin distance, in nanometers
        elif protein.boxProtocol == 'big':
            g.setup.set_boxType = 'octahedron'
            g.setup.setUseAbsBoxSize(True)
            g.setup.setAbsBoxSize('7.0')   # periodic box absolute size, in nanometers (string)
        elif protein.boxProtocol[0:3] == 'unf':
            g.setup.set_boxType = 'octahedron'
            g.setup.setUseAbsBoxSize(True)
            if protein.boxProtocol == 'unf100':
                g.setup.setAbsBoxSize('10.0') 
            if protein.boxProtocol == 'unf110':
                g.setup.setAbsBoxSize('11.0')
            if protein.boxProtocol == 'unf120':
                g.setup.setAbsBoxSize('12.0')

        thisOutDir = os.path.join(thisdir, protein.basename)
        print 'Writing equilibration directory to',thisOutDir,'...'
        if os.path.exists(thisOutDir) == False:
            os.mkdir(thisOutDir)
        os.chdir(thisOutDir)
        g.prepare(verbose=verbose, cleanup=cleanup, debug=debug, protocol=protocol, checkForFatalErrors=True, implicitOptions=implicitOptions, useTable=useTable)
     

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

    def __init__(self, pdb=None, seq=None, salt='NaCl', saltconc=0.100, pH=7.0, boxProtocol='small',
                       cosolvent = None, cosolventconc = 0.0):
        """Initialize a default data in protein object."""

        self.pdbfile  = pdb
        self.seqfile  = seq
        self.salt     = salt         # Molar concentration of salt (ionic strength)
        self.saltconc = saltconc
        self.pH       = pH
        self.boxProtocol = boxProtocol
        self.boxProtocolOptions = {'small' : '15A margin around the protein. Octahedron',
                                   'big'   : '70A absolute boxsize.  Octahedron',
                                   'unf100': '100A absolute boxsize.  Octahedron',
                                   'unf110': '110A absolute boxsize.  Octahedron',
                                   'unf120': '120A absolute boxsize.  Octahedron' }
        self.cosolvent = cosolvent
        self.cosolventconc = cosolventconc
        self.cosolventOptions = { None : 'no cosolvent',
                                 'gnd' : 'guanidinium' }

        # to be added along the way
        self.modelPDBout = None   # PDB filename to write after the MODELLER phase    gencoil(sequence, outdir, nconf=10, temperature=300.0, dimension=None, useNN='single'):
        self.mccePDBout  = None   # PDB filename to write after the MCCE phase
        self.basename = None      # Basename for output files


    def setup(self, outdir):
        """set up filenames and working/output directory."""

        # create outdir if we need to
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        
        # determine base name from the sequence file
        fields = self.seqfile.split('/')
        if fields[-1] == '':
            fields.pop()
        self.basename = fields[-1].replace('.seq','')

        self.modelPDBout = os.path.join(outdir,(self.basename+"_model.pdb"))
        self.mccePDBout  = os.path.join(outdir,(self.basename+"_mcce.pdb"))


    def getBasename(self, filename):
        """Gets the basename of any file names 'basename.suffix'."""
        extInd = filename.rindex(".")
        return filename[0:extInd]
         
    def print_info(self):
        """Prints info about the PipelineProtein object"""
        print self.info()

    def info(self):
        """Returns a string with information about this PipelineProtein object"""
        s = ''
        s = s + '%-16s = %s\n'%('pdbfile', self.pdbfile)
        s = s + '%-16s = %s\n'%('seqfile', self.seqfile)
        s = s + '%-16s = %s\n'%('salt', self.salt)
        s = s + '%-16s = %s\n'%('saltconc', self.saltconc)
        s = s + '%-16s = %s\n'%('pH', self.pH)
        s = s + '%-16s = %s\n'%('boxProtocol', self.boxProtocol)
        s = s + '%-16s = %s\n'%('modelPDB_out', self.modelPDBout)
        s = s + '%-16s = %s\n'%('mccePDBout', self.mccePDBout)
        return s
