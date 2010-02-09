
# ----------------------------------------------------------------------
# IMPORTS
# ----------------------------------------------------------------------

import sys, os, string, commands, glob, random
import tempfile

from atomSelect import *
from Filenames import *
from IonCalculator import *
from Setup import *
from Status import *
from GromacsData import * 
from MdpFile import *
from mmtools.pdbtools import *


# ----------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# CLASSES
# ----------------------------------------------------------------------

class PrepWorker(Object):
    """An object to store info and coordinate a preparation run."""
  
    def __init__(self, infile, finalOutputDir=None, finalOutputName=None, workdir=None, verbose=True): 
        """Initialize the class
    
        REQUIRED INPUTS
        infile		the stating *.pdb, *.gro, etc. file to initialize the system 
    
        OPTIONAL INPUTS
        finalOutputDir      directory to save final *.trp, *.gro files for an MD simulation. Defaults to './out'
        finalOutputName      basename (*) to name the final *.trp, *.gro files for an MD simulation. Defaults to 'out'
        workdir		directory to perform gmx preprocessing.  If not specified, work
                            will be done in a unique temporary directory.
        """
    
        self.verbose = verbose
        
        # finalOutDir
        # the name of the directory to write the final prepared gmx files
        if finalOutputDir == None:
            self.finalOutputDir = os.path.abspath(os.path.join(os.curdir,'out'))
        else:
            self.finalOutputDir = os.path.abspath(finalOutputDir)
    
        # finalOutDir
        # the root name of the files written in the final prepared gmx files    
        if finalOutputName == None:        
            self.finalOutputName = 'out'    
        else:
            self.finalOutputName = finalOutputName
    
        # workdir
        # The working directory for the simulation preparation
        if workdir == None:
            self.workdir = tempfile.mkdtemp();
        else:
            self.workdir = workdir
    
        # files
        # An object for keeping track of a series of filenames
        self.files = PrepFilenames(infile, self.workdir)
	
	# status
        self.status = PrepStatus()

        
        
        #############################3
        # Start Building the GROMACS project
        
        print "GROMACS files %s.tpr and %s.gro, etc. will be written to directory %s..." % (self.finalOutputName, self.finalOutputName, self.finalOutputDir)
    
        print "Building GROMACS project in temporary directory %s..." % self.workdir
        output = commands.getoutput('cp %s %s'%(self.files.infile, os.path.join(self.workdir, self.files.infile_basename) ))
        print output
    
        # show a listing of all the current gmx files 
        if (self.verbose):
          self.files.show()    
      
    def mdrun(self, mdpfile):
	"""Run the commands to grompp and mdrun, and shift the current self.files."""
	pass
      
    def cleanupWorkingDir(self):
        """Clean up the working directory."""
        
        print
        print 'Cleaning up temporary directory %s...'%self.workdir
        curdir = os.getcwd()
        os.chdir(self.workdir)
        for i in os.listdir(self.workdir):
            thisdir = os.path.join(self.workdir,i) 
            if os.path.isdir(thisdir):
                for j in os.listdir(thisdir):
                   os.unlink(os.path.join(thisdir,j))
                os.rmdir(thisdir) 
            else:
                print 'removing', i, '...'
                os.unlink(thisdir)
        os.chdir(curdir)
        os.rmdir(self.workdir)

    def rungmx(self, cmd, mockrun=False, checkForFatalErrors=True):
        """Execute a gromacs executable on the command line, keeping track of the command and output for the logfile"""
      
        # origdir = os.curdir
        # os.chdir( os.environ['GMXPATH'] )
        # NOTE:  DON'T do this!!!  The gmx command should be the absolute path, and we should
        # stay in the current directory.  Other wise, posre.itp gets written to the GMXPATH directory, etc.
        
        print '>>',cmd	
        if mockrun==False:
          output=commands.getoutput(cmd)
          if self.verbose: print output
          self.status.loglines.append( (cmd, output) )    
          
          if checkForFatalErrors:def 
            if output.count('Fatal error') > 0:
              print '*** There was a fatal error in the following command: ***'
              print '*********************************************************'
              print cmd
              print '*********************************************************'
              print 'Exiting...'
              sys.exit(1)
              
        # os.chdir( origdir )
        return output

    def copyFilesToOutputDir(self, useTable=None):
    
    	# copy the necessary files to start an MD run back to the original curdir
	print
	print 'Copying the *.gro, *.tpr, and production script to ',self.finalOutputDir,'...'

        # create the final output directory if it doesn't yet exist
        if not os.path.exists(self.finalOutputDir):
            os.mkdir(self.finalOutputDir)
	
	### the GRO file
	out_grofile = os.path.join(self.finalOutputDir, self.finalOutputName+'.gro')
	copycmd = 'cp %s %s'%(self.files.grofile, out_grofile)
	if (self.verbose):
	  print copycmd 
	cmdout = commands.getoutput( copycmd )
	if self.verbose==True:
	    print cmdout
	    
	### the TPR file
	out_tprfile = os.path.join(self.finalOutputDir, self.finalOutputName+'.tpr')
	copycmd = 'cp %s %s'%(self.files.tprfile, out_tprfile)
	if (self.verbose): print copycmd 
	cmdout = commands.getoutput( copycmd )
	if (self.verbose): print cmdout

	### the TOP file
	out_topfile = os.path.join(self.finalOutputDir, self.finalOutputName+'.top')
	copycmd = 'cp %s %s'%(self.files.topfile, out_topfile)
	if (self.verbose): print copycmd 
	cmdout = commands.getoutput( copycmd )
	if (self.verbose): print cmdout
	    
	### the NDX file
	out_ndxfile = os.path.join(self.finalOutputDir, self.finalOutputName+'.ndx')
	copycmd = 'cp %s %s'%(self.files.ndxfile, out_ndxfile)
	if (self.verbose): print copycmd 
	cmdout = commands.getoutput( copycmd )
	if (self.verbose): print cmdout

        ### the MDP file, if there is one... ###
        out_mdpfile = os.path.join(self.finalOutputDir, self.finalOutputName+'.mdp')
        if os.path.exists(self.files.mdpfile):
            copycmd = 'cp %s %s'%(self.files.mdpfile, out_mdpfile)
            if (self.verbose): print copycmd
            cmdout = commands.getoutput( copycmd )
            if (self.verbose): print cmdout

        ### the 'table.xvg' file, if there is one... ###
	if useTable != None:
            out_xvgfile = os.path.join(self.finalOutputDir, 'table.xvg')
            if os.path.exists(self.files.xvgfile):
                copycmd = 'cp %s %s'%(self.files.xvgfile, out_xvgfile)
                if (self.verbose): print copycmd
                cmdout = commands.getoutput( copycmd )
		if (self.verbose): print cmdout


	### the mdrun script
	equilibrate = '%s/mdrun -v -s %s -c %s '%( os.environ['GMXPATH'], out_tprfile, out_grofile )
	equilscript = os.path.join(self.finalOutputDir, 'production')
	fout = open(equilscript,'w')
	fout.write(equilibrate+'\n')  def copyFilesToOutputDir(self, useTable=None):
    
    	# copy the necessary files to start an MD run back to the original curdir
	print
	print 'Copying the *.gro, *.tpr, and production script to ',self.finalOutputDir,'...'

        # create the final output directory if it doesn't yet exist
        if not os.path.exists(self.finalOutputDir):
            os.mkdir(self.finalOutputDir)
	
	### the GRO file
	out_grofile = os.path.join(self.finalOutputDir, self.finalOutputName+'.gro')
	copycmd = 'cp %s %s'%(self.files.grofile, out_grofile)
	if (self.verbose):
	  print copycmd 
	cmdout = commands.getoutput( copycmd )
	if self.verbose==True:
	    print cmdout
	    
	### the TPR file
	out_tprfile = os.path.join(self.finalOutputDir, self.finalOutputName+'.tpr')
	copycmd = 'cp %s %s'%(self.files.tprfile, out_tprfile)
	if (self.verbose): print copycmd 
	cmdout = commands.getoutput( copycmd )
	if (self.verbose): print cmdout

	### the TOP file
	out_topfile = os.path.join(self.finalOutputDir, self.finalOutputName+'.top')
	copycmd = 'cp %s %s'%(self.files.topfile, out_topfile)
	if (self.verbose): print copycmd 
	cmdout = commands.getoutput( copycmd )
	if (self.verbose): print cmdout
	    
	### the NDX file
	out_ndxfile = os.path.join(self.finalOutputDir, self.finalOutputName+'.ndx')
	copycmd = 'cp %s %s'%(self.files.ndxfile, out_ndxfile)
	if (self.verbose): print copycmd 
	cmdout = commands.getoutput( copycmd )
	if (self.verbose): print cmdout

        ### the MDP file, if there is one... ###
        out_mdpfile = os.path.join(self.finalOutputDir, self.finalOutputName+'.mdp')
        if os.path.exists(self.files.mdpfile):
            copycmd = 'cp %s %s'%(self.files.mdpfile, out_mdpfile)
            if (self.verbose): print copycmd
            cmdout = commands.getoutput( copycmd )
            if (self.verbose): print cmdout

        ### the 'table.xvg' file, if there is one... ###
	if useTable != None:
            out_xvgfile = os.path.join(self.finalOutputDir, 'table.xvg')
            if os.path.exists(self.files.xvgfile):
                copycmd = 'cp %s %s'%(self.files.xvgfile, out_xvgfile)
                if (self.verbose): print copycmd
                cmdout = commands.getoutput( copycmd )
		if (self.verbose): print cmdout


	### the mdrun script
	equilibrate = '%s/mdrun -v -s %s -c %s '%( os.environ['GMXPATH'], out_tprfile, out_grofile )
	equilscript = os.path.join(self.finalOutputDir, 'production')
	fout = open(equilscript,'w')
	fout.write(equilibrate+'\n')
	fout.close()
	os.chmod(equilscript, 0775)
	if self.verbose==True:
	    print cmdout
	fout.close()
	os.chmod(equilscript, 0775)
	if self.verbose==True:
	    print cmdout    