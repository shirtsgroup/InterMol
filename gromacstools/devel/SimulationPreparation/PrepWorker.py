
# ----------------------------------------------------------------------
# IMPORTS
# ----------------------------------------------------------------------

import sys, os, string, commands, glob, random
import tempfile

from PrepStatus import *
from PrepFilenames import *


# ----------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# CLASSES
# ----------------------------------------------------------------------

class PrepWorker():
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
	
	self.checkForFatalErrors = True   # these will be set by the prepare method
	self.mockrun = False
	
	# A set of paths to look in for mdpfiles 
	self.mdpPaths = []
	self.mdpPaths.append( os.path.join( os.environ['MMTOOLSPATH'], 'gromacstools/devel/SimulationPreparation/mdp') )

        
        #####################################
        # Start Building the GROMACS project
        
        print "GROMACS files %s.tpr and %s.gro, etc. will be written to directory %s..." % (self.finalOutputName, self.finalOutputName, self.finalOutputDir)
    
        print "Building GROMACS project in temporary directory %s..." % self.workdir
        output = commands.getoutput('cp %s %s'%(self.files.infile, os.path.join(self.workdir, self.files.infile_basename) ))
        print output
    
        # show a listing of all the current gmx files 
        if (self.verbose):
          self.files.show()    
	  
	# The PrepWorker needs to know where the directory of mdpfile is

	
    def prepare(self, systemSetup, cleanup=True, verbose=False, debug=False, checkForFatalErrors=True, mockrun=False):
	"""
	REQUIRED ARGUMENTS
	system  = a SystemSetup() object containing the system info: ion concentrations, forcefields, etc. 
    
	OPTIONAL ARGUMENTS
	cleanup             if True, delete the temporary pre-processing directory, and erase temporary files (Default: True)
	verbose             if True, print extra, more detailed progress statements (Default: False)
	debug               if True, turn on print statements for debugging (Default: False)    
	checkForFatalErrors Exits the program if any of the gmx produce "Fatal error" in the standard output (Default: True)	    
	mockrun             if True, prints out the gromacs commands, but does does execture them (Default: False)
				    
	 
	OUTPUT
      
	The following files will be written to the outdir:
	  
	    out.tpr		GROMACS tpr file
	    out.gro		GROMACS gro (structure) file
	    equilibrate*	executable script to call mdrun 
	    history.log         logfile containing all the commands used to prepare the simualtion 
	    output.log          logfile containg the output of the pre-processing programs
	"""	
	self.system = systemSetup
	
        # Overwrite the set of paths to look in for mdpfile -- instead use the user-specified list
	self.mdpPaths = systemSetup.mdpPaths 
	 
    
	self.checkForFatalErrors = checkForFatalErrors
	self.mockrun = mockrun
	
	# remember our original location before we descend into the working directory
        OriginalLocationDir = os.path.abspath( os.curdir )
	
	# From now on, do work in the working directory
	os.chdir( self.workdir )
	print 'Preparing simulation in directory %s...'%self.workdir
	

	
	### Put your preparation steps HERE!  ### 
	
	
        # Finalize the preparation: 1) copy files to output directory, 2) write logfile, and 3) clean up temp files
	self.finalizePreparation(cleanTempFiles=cleanup)   # this is defined in the parent class
	

	os.chdir( OriginalLocationDir )	### change back out of the working directory to our original location
	return

	
	
	  
    def finalizePreparation(self, cleanTempFiles=False): 
        """Finalize the preparation by:
	    1) copying files to the output directoy
	    2) writing log files
	    3) doing cleanup
	"""

        # copy the final files to the final output directory
	self.copyFilesToOutputDir(useTable=None)
	    
        # write logfiles
	self.status.writeLog( os.path.join(self.finalOutputDir,'output.log') )
        self.status.writeHistory( os.path.join(self.finalOutputDir,'history.log') )
	    
	    
        # Clean up the working directory
        if (cleanTempFiles):
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

	
	  
      
    def mdrun(self, mdpfile, useIndexFile=False):
	"""Run the commands to grompp and mdrun, and shift the current self.files.
	NOTE: This function assumes that the desired product of the mdrun is the output conformation."""
	
	# First, build the tpr using grompp
	self.buildTpr(mdpfile, useIndexFile=useIndexFile)
	
	# Then, run the mdrun job
	myjob = '%s/mdrun -v -s %s -c %s '%(os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
	self.rungmx( myjob )
	self.files.increment_gro()    # must increment filename for any new gmx file 
	

    def buildTpr(self, mdpfile, useIndexFile=False):
	"""Run the commands to grompp, but DON'T perform the mdrun.
	NOTE: This function assumes that the desired product of the mdrun is the output conformation."""
	
	# Find the mdpfile in our list of possible paths
	myMdpfile = self.findMdpfile(mdpfile)
	
	# make a tpr file with grompp
	if not useIndexFile:
	    grompp = '%s/grompp -f %s -c %s -o %s -p %s'%(os.environ['GMXPATH'], myMdpfile, self.files.grofile, self.files.next_tpr(), self.files.topfile)
	else:
	    grompp = '%s/grompp -f %s -c %s -o %s -p %s -n %s'%(os.environ['GMXPATH'], myMdpfile, self.files.grofile, self.files.next_tpr(), self.files.topfile, self.files.ndxfile)
	self.rungmx( grompp )	
	self.files.increment_tpr()    # must increment filename for any new gmx file 


	

    def findMdpfile(self, mdpfile):
	"""Look in mdpPaths to try and find the mdpfile.
	RETURNS
	The path of the mdpfile.  (An exception will be raised if it can't find it.)
	"""
	
	myMdpFile = None
	for path in self.mdpPaths:
	    myMdpFile = os.path.join(path,mdpfile)
	    if os.path.exists(myMdpFile):
		break
        if myMdpFile == None:
	    print 'Cannot find file the mdpfile %s in the set of paths %r'%(mdpfile, self.mdpPaths)
	    raise MdpFileNotFound
	
	return os.path.abspath(myMdpFile)
	

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
	

    def rungmx(self, cmd):
        """Execute a gromacs executable on the command line, keeping track of the command and output for the logfile"""
      
        # origdir = os.curdir
        # os.chdir( os.environ['GMXPATH'] )
        # NOTE:  DON'T do this!!!  The gmx command should be the absolute path, and we should
        # stay in the current directory.  Other wise, posre.itp gets written to the GMXPATH directory, etc.
        
        print '>>',cmd	
        if self.mockrun==False:
          output=commands.getoutput(cmd)
          if self.verbose: print output
          self.status.loglines.append( (cmd, output) )    
          
          if self.checkForFatalErrors:
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
	fout.write(equilibrate+'\n')
	

class MdpFileNotFound(Exception):
    pass