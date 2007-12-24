# gromacstools/System.py
#
# REQUIREMENTS
#
#    1) GROMACS ( http://www.gromacs.org/ ) must be installed on your machine, 
#       and the following Environment variables must be defined (in your .bash_profile, e.g.):
# 
#        GMXLIB             the pathname of the gmx parameter files.  
#        GMXPATH            the pathname of the Gromacs executables  (This also needs to be in your PATH)
#        GMXPATH_PS3        the pathname of the Gromacs executables used for building PlayStation3 WUs
#        MMTOOLSPATH        the pathname of the mmtools library
#
#    2) The AMBER ports for GROMACS ( http://folding.stanford.edu/ffamber/ ) by Sorin, Park
#       and Pande must be installed, with the *.itp, *.rtp, etc., files in the $GMXLIB directory (.../share/top )
#
#    3) If you need to simulate any proteins with a norleucine (NLE) residue, you will need to insert the
#       following lines in the *.hdb and *.rtp files for the ffamber gmx port forcefields:
#
#           ffamber/norleucine.hdb
#           ffamber/norleucine.rtp
#
# ---------------------------------------------------------------------
# AUTHORS
#
# Written by Vincent Voelz, Stanford University (c) 2007
# ----------------------------------------------------------------------
# MODIFICATION HISTORY
#
# 10/01/07 VAV - Added racecar3 protocol
#                fixed the -H14 flag in pdb2gmx
#                in make_ndx() created (if protocol=='racecar3') atomgroup "bath" which is the solvent and the ions, to be controlled
# 09/02/07 VAV - Added make_PS3_visualization_groups(), and helper function find_residue_groups(), which is shared now by make_ndx() and make_PS3...
# 08/28/07 VAV - Fixed make_ndx() in the case that there are NO solvent residues (e.g.)
# 08/27/07 VAV - Added implicitSolvationSteps for PS3: GB, AGBNP, etc.
# 07/08/07 VAV - RE-fixed the make_ndx() function to include ions!!!!
# 06/26/07 GRB - added code to allow use of absolute box size

# ----------------------------------------------------------------------
# TO DO
# * Right now user-tables DO NOT work with non-implicit preparations!!!  Tis is easy to fix -- we just need to write a modified *.mdp file and table.xvg
#   evbery time we need to use a mdp file
# * write a better MDP 'protocol' object?  right now protocols are 'cases'
#
# *** MOST IMPORTANT ***
# * Implement Exceptions!!!  One failed protein should cripple the rest of the pipeline.  
# *** IMPORTANT BUT LESS SO ***
# * make_ndx() doesn't look to see all the Ion types available.  Fix for robustness 
# * PROPER documentation of all functions according to style guide
# * Incorporate Units.py
# * a better way of keeping track of gromacs preprocessed files!!!
# * --- set() and get() for mdp attribyutes 
# * an elegant way of indexing subgroups
# * Support for GROMACS 3.3.1
# * Support for more salts in GromacsSystemSetup
# * Make ioncalc a more useful object
# ---- compile and keep all possible ion pairs, parms 
#
# FUTURE
# * Full python wrapper for gromacs functions? bundled API?

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


# ----------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# CLASSES
# ----------------------------------------------------------------------


class System:

  def __init__(self, infile, finalOutputDir=None, finalOutputName=None, workdir=None, useff='ffamber03',  version='3.1', verbose=True): 
    """A class to initialze, minimize, and equilibrate Gromacs simulations.

    REQUIRED INPUTS
    infile		the stating *.pdb, *.gro, etc. file to initialize the system 

    OPTIONAL INPUTS
    finalOutputDir      directory to save final *.trp, *.gro files for an MD simulation. Defaults to './out'
    finalOutputName      basename (*) to name the final *.trp, *.gro files for an MD simulation. Defaults to 'out'
    workdir		directory to perform gmx preprocessing.  If not specified, work
                        will be done in a unique temporary directory.
    ff			string to specify forcefield (default is 'ffamber99p').  'ffamber03' also supported
    version             gmx version -- supported: '3.1.4' ('3.3.1' COMING SOON - not that diffferent)
    
"""

    # initialize status object
    self.status = Status()

    # initialize setup object
    self.setup = Setup()
    
    # finalOutDir is the name of the directory to write the final prepared gmx files
    if finalOutputDir == None:
        self.finalOutputDir = os.path.abspath(os.path.join(os.curdir,'out'))
    else:
        self.finalOutputDir = os.path.abspath(finalOutputDir)

    # finalOutDir - the name of the directory to write the final prepared gmx files    
    if finalOutputName == None:        
        self.finalOutputName = 'out'    
    else:
        self.finalOutputName = finalOutputName

    print "GROMACS files %s.tpr and %s.gro, etc. will be written to directory %s..." % (self.finalOutputName, self.finalOutputName, self.finalOutputDir)

    # process workdir
    if workdir == None:
        self.workdir = tempfile.mkdtemp();
    else:
        self.workdir = workdir

    # initialize files object  
    self.files = Filenames(infile, self.workdir)

    # Start Building the GROMACS project
    print "Building GROMACS project in temporary directory %s..." % self.workdir
    output = commands.getoutput('cp %s %s'%(self.files.infile, os.path.join(self.workdir, self.files.infile_basename) ))
    print output

    
    # process other options
    self.useff = useff
    self.version = version
    self.verbose = verbose
    self.protocol = 'default'    # keyword for Preparation protocol 

    # find available forcefields with this GROMACS installation
    [self.forcefields, self.forcefieldCodes] = self.getForcefields()
    
    # check if user-specifiedforcefield is available
    print 'Trying to use forcefield', self.useff, '...'
    if self.forcefields.has_key(useff):
        print 'Using', self.useff, ' with key',self.forcefieldCodes[self.useff],':',self.forcefields[self.useff]
    else:
	print 'Cannot find forcefield', self.useff,'\n'
	self.printForcefields()
	print '\nExiting....'
        sys.exit(1)
  	
       
    # to keep track of the total charge in the system
    # NOTE: this will be calculated from the pdb2gmx output text
    self.totalChargeBeforeIons = 0.0    
    
    # initialize ionCalculator
    ### This is needed for calculating numbers of counterions from salt concentrations
    self.ioncalc = IonCalculator(self.useff)


    # show a listing of all the current gmx files 
    if (verbose):
      self.files.show()

 
  def printTwoColumnStrings(self,s1,s2):
    """Prints two strings in a fixed-width column."""
    print '%-36s %s'%(s1,s2) 
    return


  def getForcefields(self):
    """Get the forcefields and codes from FF.dat; produce error if this is not found."""

    try:
      fin = open(os.path.join(self.files.GMXLIB,'FF.dat'),'r')
    except IOError:
      print "File FF.dat cannot be found in GMXLIB=%s   ..... \nExiting....\n"%self.files.GMXLIB
      sys.exit(1)
    forcefields = {}
    forcefieldCodes = {}
    lines = fin.readlines()
    numff = int(lines.pop(0)) 
    for i in range(0,numff):
        fields = lines[i].split()
	key = fields.pop(0)
	forcefields[ key ] = string.joinfields(fields)
	forcefieldCodes[ key ] = str(i)
    return forcefields, forcefieldCodes
  
   
  def printForcefields(self):
    """Prints the available forcefields."""

    print 'Available forcefields in $GMXLIB =',os.path.join(self.files.GMXLIB,'FF.dat')
    for key,value in self.forcefields.items():
      print '%-16s%s'%(key,value)
    
    

  def prepare(self, cleanup=True, verbose=False, debug=False, protocol=None, checkForFatalErrors=True, mockrun=False, implicitOptions=None, useTable=None):
    
    """Prepare an equilibration simulation for this Gromacs system, namely:
    
    1.  pdb2gmx
    2.  edit topology file (if necessary) to use tip3p water
    3.  solvate the box
    4.  minimize the box
    5.  add ions to the solvated box
    6.  minimize once again
    7.  prepare a molecular dynamics simulation for equilibration...

    OPTIONAL ARGUMENTS
    
    ARG                 DEFAULT     DESCRIPTION
    cleanup             True        if True, delete the temporary pre-processing directory, and erase temporary files 
    verbose             False       if True, print extra, more detailed progress statements 
    debug               False       if True, turn on print statements for debugging    
    protocol            'racecar2'  Currently 'quicktest', 'racecar2', 'implicitPS3'  and 'default' supported:
				    * 'default':   minimize.mdp -> equilibrate.mdp
				    * 'quicktest':   quicktest_min1.mdp -> quicktest_equil1.mdp -> quicktest_equil2.mdp -> quicktest_grompp.mdp
				    * 'racecar2':  racecar2_min1.mdp -> racecar2_equil1.mdp -> racecar2_equil2.mdp -> racecar2_grompp.mdp
                                    * 'racecar3':  racecar33_min1.mdp -> racecar3_equil1.mdp -> racecar3_equil2.mdp -> racecar3_grompp.mdp

    checkForFatalErrors	True        Exits the program if any of the gmx produce "Fatal error" in the standard output	    
    mockrun             False       if True, prints out the gromacs commands, but does does execture them  
    implicitOptions     None        can be a MdpFile.ImplicitOptions() object used to pass in PS3 implicit setup information   
    useTable            None        can be one of the TableMaker.DielectricTable derived classes.  Supported as of 9/2007:
				            TableMaker.Switch_Shift(), TableMaker.ReactionField(), TableMaker.SF3()
			        
     
    OUTPUT
  
    The following files will be written to the outdir:
      
        out.tpr		GROMACS tpr file
        out.gro		GROMACS gro (structure) file
        equilibrate*	executable script to call mdrun 
        history.log     logfile containing all the commands used to prepare the simualtion 
        output.log      logfile containg the output of the pre-processing programs

"""

    # remember our original location before we descend into the working directory
    OriginalLocationDir = os.path.abspath( os.curdir )
    
    # initialize input parameters
    self.cleanup = cleanup          # Erase the temporary files when done
    self.verbose = verbose          # Print extra, more detailed progress statements 
    self.debug   = debug            # Turns on print statements for debugging
    if protocol == None:            # Current only 'default' and 'racecar2' are supported:
      self.protocol = 'default'     #     'default':   minimize.mdp -> equilibrate.mdp -> simulate.mdp 
    else:                           #     'racecar2':  racecar2_min1.mdp -> racecar2_equil1.mdp -> racecar2_equil2.mdp -> racecar2_grompp.mdp  
      self.protocol = protocol      #     'racecar3':  racecar3_min1.mdp -> racecar3_equil1.mdp -> racecar3_equil2.mdp -> racecar3_grompp.mdp  
                                    #     'quicktest':  quicktest_min1.mdp -> quicktest_equil1.mdp -> quicktest_equil2.mdp -> quicktest_grompp.mdp 
                                    #     'implicitPS3'
    self.checkForFatalErrors = checkForFatalErrors
    self.mockrun = mockrun
    
    
    # Choose a set of *.mdp files according to the specified protocol
    if self.protocol=='racecar2':
        self.files.set_mdpMinimization('racecar2_min1.mdp')
	self.files.set_mdpEquilibration('racecar2_equil1.mdp')    # equil2 comes later....    
        self.files.set_mdpSimulation('racecar2_grompp.mdp')  

    if self.protocol=='racecar3':
        self.files.set_mdpMinimization('racecar3_min1.mdp')
	self.files.set_mdpEquilibration('racecar3_equil1.mdp')    # equil2 comes later....    
        self.files.set_mdpSimulation('racecar3_grompp.mdp')  

    if self.protocol=='quicktest':
        self.files.set_mdpMinimization('quicktest_min1.mdp')
        self.files.set_mdpEquilibration('quicktest_equil1.mdp')    # equil2 comes later....    
        self.files.set_mdpSimulation('quicktest_grompp.mdp')

    if self.protocol=='implicitPS3_AGBNP':
        self.files.set_mdpMinimization('implicitPS3_AGBNP_min1.mdp')
        self.files.set_mdpEquilibration('implicitPS3_AGBNP_equil1.mdp')
        self.files.set_mdpSimulation('implicitPS3_AGBNP_grompp.mdp')
	os.environ['GMXPATH'] = os.environ['GMXPATH_PS3']

    if self.protocol=='implicitPS3_GB':
        self.files.set_mdpMinimization('implicitPS3_GB_min1.mdp')
        self.files.set_mdpEquilibration('implicitPS3_GB_equil1.mdp')
        self.files.set_mdpSimulation('implicitPS3_GB_grompp.mdp')
        os.environ['GMXPATH'] = os.environ['GMXPATH_PS3']

    if self.protocol=='implicitPS3_StillGB':
        self.files.set_mdpMinimization('implicitPS3_StillGB_min1.mdp')
        self.files.set_mdpEquilibration('implicitPS3_StillGB_equil1.mdp')
        self.files.set_mdpSimulation('implicitPS3_StillGB_grompp.mdp')
	os.environ['GMXPATH'] = os.environ['GMXPATH_PS3']

    if self.protocol=='implicitPS3_Sigmoidal':
        self.files.set_mdpMinimization('implicitPS3_Sigmoidal_min1.mdp')
        self.files.set_mdpEquilibration('implicitPS3_Sigmoidal_equil1.mdp')
        self.files.set_mdpSimulation('implicitPS3_Sigmoidal_grompp.mdp')
        os.environ['GMXPATH'] = os.environ['GMXPATH_PS3']
	
	
    # If there is a user-defined coulomb table, 	

    # From now on, do work in the working directory
    os.chdir( self.workdir )
    print 'Preparing simulation in directory %s...'%self.workdir
    		
    # prepare the simultion 	
    if self.version == '3.1':	

        # pdb2gmx, and useTIP3P
        pdb2gmx = 'echo %s | %s/pdb2gmx -f %s -o %s -p %s -ignh -H14'%(self.forcefieldCodes[self.useff], os.environ['GMXPATH'], self.files.infile, self.files.grofile, self.files.topfile)
	self.rungmx( pdb2gmx, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors)
        self.useTIP3P( self.files.topfile )
	
	pdb2gmxlines = self.status.loglines[-1][1].split('\n')
	for line in pdb2gmxlines:
	  if line[0:12] == 'Total charge':
	    self.totalChargeBeforeIons = float(line.split()[2])
	    
        # make a (periodic boundary conditions) box
        if self.setup.useAbsBoxSize:
	    editconf = '%s/editconf -bt %s -f %s -o %s -box %s'%(os.environ['GMXPATH'], self.setup.boxType, self.files.grofile, self.files.next_gro(), self.setup.absBoxSize )
        else:
       	    editconf = '%s/editconf -bt %s -f %s -o %s -d %s'%(os.environ['GMXPATH'], self.setup.boxType, self.files.grofile, self.files.next_gro(), self.setup.boxSoluteDistance )
	self.rungmx( editconf, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
	self.files.increment_gro()    # must increment filename for any new gmx file 

        if self.protocol.count('implicitPS3') > 0:
            self.implicitSolvationPreparationSteps(implicitOptions=implicitOptions, useTable=useTable)
        else:

            if self.setup.cosolvent == 'gnd':
                # solvate the box with water and guanidinium 
	        self.Cosolvate()
                editconf = '%s/genbox -cp %s -cs ffamber_tip3p.gro -o %s -p %s'%(os.environ['GMXPATH'], self.files.grofile, self.files.next_gro(), self.files.topfile)
	        self.rungmx( editconf, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
	        self.files.increment_gro()    # must increment filename for any new gmx file 
            else:
                # solvate the box with plain water 
                editconf = '%s/genbox -cp %s -cs ffamber_tip3p.gro -o %s -p %s'%(os.environ['GMXPATH'], self.files.grofile, self.files.next_gro(), self.files.topfile)
                self.rungmx( editconf, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
                self.files.increment_gro()    # must increment filename for any new gmx file 

            # minimize the system
            ### make a tpr file with grompp
            grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.mdpfile_Minimization, self.files.grofile, self.files.tprfile, self.files.topfile)
	    self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )	
            ### run minimization
            minimize = '%s/mdrun -v -s %s -c %s '%(os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
	    self.rungmx( minimize, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
	    self.files.increment_gro()    # must increment filename for any new gmx file 
	
            # do the rest of the preparation steps
            self.postSolvationPreparationSteps()

        # copy the final files to the final output directory
	self.copyFilesToOutputDir(useTable=useTable)
	    
        # write logfiles
	self.status.writeLog( os.path.join(self.finalOutputDir,'output.log') )
        self.status.writeHistory( os.path.join(self.finalOutputDir,'history.log') )
	    
	    
        # Clean up the working directory
        if (cleanup):
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

	
    # elif self.version == '3.3':	
      	
        # pdb2gmx, and useTIP3P
	# pdb2gmx = 'echo %s | pdb2gmx -f %s -o %s -p %s -water tip3p -ignh'%(self.forcefieldCodes[self.useff], self.infile, self.grofile_prep, self.topfile_prep)
	# self.rungmx( pdb2gmx, mockrun=self.mockrun  )
	
        # make a cubic (periodic boundary conditions) box
	### FILL THIS IN! ###
	
        # solvate the box
	# editconf = 'genbox -cp %s -cs ffamber_tip3p.gro -o %s -p %s'%(self.grofile_box, self.grofile_solvated, self.topfile_prep)
	# self.rungmx( editconf, mockrun=self.mockrun  )

        # minimize the system
        ### FILL IN!
	

    else:
        print 'This version of gromacs (',self.version,') is not supported.'
	sys.exit(1)

    os.chdir( OriginalLocationDir )	### change back out of the working directory to our original location
    return

  
  def prepareFromSolvated(self, cleanup=True, verbose=False, debug=False, protocol=None, checkForFatalErrors=True, mockrun=False):
    """Preare simulation like prepare(), but start from a Solvated PDB or *.gro self.infile"""    

    self.cleanup = cleanup          # Erase the temporary files when done
    self.verbose = verbose          # Print extra, more detailed progress statements 
    self.debug   = debug            # Turns on print statements for debugging
    if protocol == None:            # Current only 'default', 'racecar2', and 'racecar3' are supported:
      self.protocol = 'default'     #     'default':   minimize.mdp -> equilibrate.mdp -> simulate.mdp 
    else:                           #     'racecar2':  racecar2_min1.mdp -> racecar2_equil1.mdp
      self.protocol = protocol      #     'quicktest': quicktest_min1.mdp -> quicktest_equil1.mdp -> quicktest_equil2.mdp -> quicktest_grompp.mdp
    self.checkForFatalErrors = checkForFatalErrors
    self.mockrun = mockrun
    
    if self.protocol=='racecar2':
        self.files.set_mdpMinimization('racecar2_min1.mdp')
	self.files.set_mdpEquilibration('racecar2_equil1.mdp')    # equil2 comes later....    
        self.files.set_mdpSimulation('racecar2_grompp.mdp')  

    if self.protocol=='racecar3':
        self.files.set_mdpMinimization('racecar3_min1.mdp')
	self.files.set_mdpEquilibration('racecar3_equil1.mdp')    # equil2 comes later....    
        self.files.set_mdpSimulation('racecar3_grompp.mdp')  

    elif self.protocol=='quicktest':
        self.files.set_mdpMinimization('quicktest_min1.mdp')
	self.files.set_mdpEquilibration('quicktest_equil1.mdp')    # equil2 comes later....    
        self.files.set_mdpSimulation('quicktest_grompp.mdp')  

    os.chdir( self.workdir )
    print 'Preparing simulation in directory %s...'%self.workdir

    # convert the solvated PDB or gro to a gromacs *.gro and *.top 
    # pdb2gmx, and useTIP3P
    pdb2gmx = 'echo %s | %s/pdb2gmx -f %s -o %s -p %s -ignh -H14'%(self.forcefieldCodes[self.useff], os.environ['GMXPATH'], self.files.infile, self.files.grofile, self.files.topfile)
    self.rungmx( pdb2gmx, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors)
    self.useTIP3P( self.files.topfile )
    pdb2gmxlines = self.status.loglines[-1][1].split('\n')
    for line in pdb2gmxlines:
        # for a system with water the last line  of the pdb2gmx output should read:
        # 'Total charge in system 3.000 e' 
        if line[0:22] == 'Total charge in system':
            self.totalChargeBeforeIons = float(line.split()[4]) 

    # minimize the system        ### make a tpr file with grompp
    grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.mdpfile_Minimization, self.files.grofile, self.files.tprfile, self.files.topfile)
    self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
    ### run minimization
    minimize = '%s/mdrun -v -s %s -c %s '%( os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
    self.rungmx( minimize, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
    self.files.increment_gro()    # must increment filename for any new gmx file 


    # do all the post-solvation steps
    self.postSolvationPreparationSteps()  


  def Cosolvate(self):
      """Use the genbox feature to co-solvate one solvent molecule at a time....."""

      from mmtools.pdbtools import *      # need for countSolventResiduesFromGrofile()

      # Solvate the box with water to see how many will fit in the box 
      editconf = '%s/genbox -cp %s -cs ffamber_tip3p.gro -o %s -p %s'%(os.environ['GMXPATH'], self.files.grofile, self.files.next_gro(), self.files.topfile)
      self.rungmx( editconf, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
      [ nClres, nNares, nsolres ] = countSolventResiduesFromGrofile(self.files.next_gro())
      # don't update the file increments, we just need the current file to count the water 
      #self.files.increment_gro()    # must increment filename for any new gmx file 
      print '**** nClres, nNares, nsolres ****', nClres, nNares, nsolres

      water_moles_per_L = 55.55
      cosolv_fraction = self.setup.cosolventconc/(self.setup.cosolventconc + water_moles_per_L) 
      ncosolv = int(cosolv_fraction*float(nsolres))

      # let's try to conserve the total mass, in the case the cosolvent is much larger than the solvent (water)   
      water_mw = 18.0
      total_mass = float(nsolres*water_mw)
      water_grofile = os.path.join(os.environ['GMXLIB'],'ffamber_tip3p_single.gro')
      if self.setup.cosolvent == 'gnd':
          cosolv_mw = 60.07
          cosolv_grofile = os.path.join(os.environ['GMXLIB'],'GND.gro')
      else:
          print 'Non gnd cosolvents not currently supported.'
          sys.exit(1)
      unscaled_mass  = ncosolv*cosolv_mw + (nsolres - ncosolv)*water_mw
      scaling_factor = total_mass/unscaled_mass

      ncosolv = int(scaling_factor*float(ncosolv))
      nwater = int(scaling_factor*float(nsolres))

      print 'Adding: ncosolv =', ncosolv, '    nwater =', nwater
      print '*******************************'

      waters_added = 0
      cosolv_added = 0
      while (ncosolv + nwater) > 0:     
        r = random.random()
        print 'random number r =', r
        print 'cosolv_fraction', cosolv_fraction
        #if ((r < cosolv_fraction) | nwater == 0) & (ncosolv > 0):
        if (r < cosolv_fraction):
          cosolv_added += 1
          # add cosolv molecule 
          print '### Adding cosolvent %d of %d total needed ###'%(cosolv_added,(ncosolv + nwater))
          editconf = '%s/genbox -cp %s -nmol %d -try 100 -seed %d -ci %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.grofile, 1, (ncosolv+nwater), cosolv_grofile, self.files.grofile, self.files.topfile)
          self.rungmx( editconf )
          ncosolv -= 1
        else: 
          # add water -- THIS SUCKS.  The 3.1.4 genbox only works when inserting ONE molecule, so let's get loop over inserting one solvengt molecule at a time...    
          waters_added += 1
          print '### Adding water %d of %d total needed ###'%(waters_added,(ncosolv + nwater))
          editconf = '%s/genbox -cp %s -nmol %d -try 100 -seed %d -ci %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.grofile, 1, (ncosolv+nwater), water_grofile, self.files.grofile, self.files.topfile)
          self.rungmx( editconf )
          nwater -= 1

        # round up those damn '#*' backup files and delete them!
        if (ncosolv+nwater) % 100 == 0:
            fields = self.files.grofile.split('/')
            print 'fields', fields
            fields.pop()
            wkdirglob = string.join(fields,'/') + '/#*'
            print 'wkdirglob'
            damnedBackupFiles = glob.glob(  wkdirglob )
            print 'damnedBackupFiles',damnedBackupFiles
            for dbf in damnedBackupFiles:
                os.unlink(dbf)


    
  def prepareResolvatedFromTemplate(self, grofile, grotemplate, growater, protocol='racecar3', checkForFatalErrors=True, mockrun=False, cleanup=False):
      """
      REQUIRED ARGUMENTS
      grofile                  a *.gro file with only the protein residues, and the correct periodic box dimensions
                               The residue names are allowed to clash with the protonation state, which will be reset using 'pdb2gmx -ignh' anyway
      grotemplate              a template grofile containing the correct numbers of ions and solvent molecules   
      growater                 a *.gro file for a single molecule
      
      """

      from mmtools.pdbtools import *

      self.checkForFatalErrors = checkForFatalErrors
      self.mockrun = mockrun
      self.debug = False

      # get the grofile in the right protonation state
      pdb2gmx = 'echo %s | %s/pdb2gmx -f %s -o %s -p %s -ignh -H14'%(self.forcefieldCodes[self.useff], os.environ['GMXPATH'], grofile, self.files.grofile, self.files.topfile)
      self.rungmx( pdb2gmx )
      self.useTIP3P( self.files.topfile )
      
      # center the protein in the box
      editconf = '%s/editconf -f %s -c -o %s'%(os.environ['GMXPATH'], self.files.grofile, self.files.next_gro()  )
      self.rungmx( editconf, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
      self.files.increment_gro()    # must increment filename for any new gmx file 
      
      # fill the coil grofile with the right number of ion and SOL residues from the grofile template
      [ nClres, nNares, nsolres ] = countSolventResiduesFromGrofile(grotemplate)
      totalsolventres = (nsolres + nClres + nNares)
            
      # solvate the box 
      editconf = '%s/genbox -cp %s -cs ffamber_tip3p.gro -o %s -p %s'%(os.environ['GMXPATH'], self.files.grofile, self.files.next_gro(), self.files.topfile)
      self.rungmx( editconf, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
      self.files.increment_gro()    # must increment filename for any new gmx file 
      [ nClres_sofar, nNares_sofar, nsolres_sofar ] = countSolventResiduesFromGrofile(self.files.grofile)
      totalsolventres_sofar = (nsolres_sofar + nClres_sofar + nNares_sofar)

      print 'totalsolventres', totalsolventres
      print 'totalsolventres_sofar', totalsolventres_sofar
      
      nsol_needed = (totalsolventres - totalsolventres_sofar)
      
      if nsol_needed > 0:
	  # add water -- THIS SUCKS.  The 3.1.4 genbox only works when inserting ONE molecule, so let's get loop over inserting one solvengt molecule at a time...    
	  for nsol in range(0, nsol_needed):
	      print '### Adding water %d of %d total needed ###'%(nsol+totalsolventres_sofar+1,totalsolventres)
	      editconf = '%s/genbox -cp %s -nmol %d -try 100 -seed %d -ci %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.grofile, 1, nsol, growater, self.files.grofile, self.files.topfile)
	      self.rungmx( editconf )
	      # round up the '#*' backup files and delete them!
	      if nsol % 100 == 0:
		fields = self.files.grofile.split('/')
		print 'fields', fields
		fields.pop()
		wkdirglob = string.join(fields,'/') + '/#*'
		print 'wkdirglob'
		damnedBackupFiles = glob.glob(  wkdirglob )
		print 'damnedBackupFiles',damnedBackupFiles
		for dbf in damnedBackupFiles:
		    os.unlink(dbf)

      elif (nsol_needed < 0):
	  print '### We need to get RID of %d waters ###'%(-1*nsol_needed)
	  g = GromacsStructureFromGrofile(self.files.grofile)
	  [ nClres_sofar, nNares_sofar, nsolres_sofar ] = countSolventResiduesFromGromacsStructure(g)
	  while nsolres_sofar > totalsolventres:
	    g.removeatom(-1)
	    print '### Now %d water residues - we want %d total ###'%(nsolres_sofar,totalsolventres)
	    
	  g.write( self.files.next_gro() )
	  self.files.increment_gro()
	  
      else:
	  print 'Perfect number of waters!'
	  
      # minimize the system
      ### make a tpr file with grompp
      grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.mdpfile_Minimization, self.files.grofile, self.files.tprfile, self.files.topfile)
      self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )	
      ### run minimization
      minimize = '%s/mdrun -v -s %s -c %s '%(os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
      self.rungmx( minimize, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
      self.files.increment_gro()    # must increment filename for any new gmx file 


      # add ions
      # create a make_ndx *.ndx file and 
      self.make_ndx(self.files.grofile, self.files.ndxfile, groups=['protein', 'sol'] )
      # the solvent groupu is guaranteed to be group "1"
      solgroup = 1
      genion = 'echo %d | %s/genion -s %s -n  %s -o %s -pname %s -np %d -pq %d -nname %s -nn %d -nq %d -g genion.log'%(solgroup, os.environ['GMXPATH'], self.files.tprfile,  self.files.ndxfile, self.files.next_gro(), self.setup.positiveIonName, nNares, self.setup.positiveIonCharge, self.setup.negativeIonName, nClres, self.setup.negativeIonCharge)
      self.rungmx( genion, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
      self.files.increment_gro()    # must increment filename for any new gmx file 

      ### generate a new topolgy file for the ion-ated grofile  
      pdb2gmx = 'echo %s | %s/pdb2gmx -f %s -o %s -p %s -ignh -H14'%(self.forcefieldCodes[self.useff], os.environ['GMXPATH'], self.files.grofile, self.files.next_gro(), self.files.next_top())
      self.rungmx( pdb2gmx, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
      self.files.increment_gro()    # must increment filename for any new gmx file 
      self.files.increment_top()    # must increment filename for any new gmx file 
      self.useTIP3P( self.files.topfile )
	      
      if (self.protocol == 'racecar2') | (self.protocol == 'racecar3') | (self.protocol == 'quicktest'):
  
	  # run minimimzation once more with the new ions
	  ### make a tpr file with grompp
	  grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.mdpfile_Minimization, self.files.grofile, self.files.next_tpr(), self.files.topfile)
	  self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	  self.files.increment_tpr()    # must increment filename for any new gmx file 
	  ### run minimization
	  minimize = '%s/mdrun -v -s %s -c %s '%( os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
	  self.rungmx( minimize, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	  self.files.increment_gro()    # must increment filename for any new gmx file 
  
	  # run equilibration phase 1: just water with frozen protein.  
	  ### make a tpr file with grompp
	  grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.mdpfile_Equilibration, self.files.grofile, self.files.next_tpr(), self.files.topfile)
	  self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	  self.files.increment_tpr()    # must increment filename for any new gmx file 
	  ### run equilibration
	  equilibrate = '%s/mdrun -v -s %s -c %s '%(os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
	  self.rungmx( equilibrate, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	  self.files.increment_gro()    # must increment filename for any new gmx file 
	  
	  # Make an index for the separately temperature groups "protein" and "sol"
	  self.make_ndx( self.files.grofile, self.files.ndxfile )
	  
	  # run equilibration phase 2: protein and sol.  
	  if (self.protocol == 'racecar2'):
	      self.files.set_mdpEquilibration('racecar2_equil2.mdp')
	  elif (self.protocol == 'racecar3'):
	      self.files.set_mdpEquilibration('racecar3_equil2.mdp')
	  elif (self.protocol == 'quicktest'): 
	      self.files.set_mdpEquilibration('quicktest_equil2.mdp')
	  else:
	      raise "Protocol Error!"
  
	  ### make a tpr file with grompp
	  grompp = '%s/grompp -f %s -c %s -o %s -p %s -n %s'%(os.environ['GMXPATH'], self.files.mdpfile_Equilibration, self.files.grofile, self.files.next_tpr(), self.files.topfile, self.files.ndxfile)
	  self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	  self.files.increment_tpr()    # must increment filename for any new gmx file 
	  ### run equilibrate
	  equilibrate = '%s/mdrun -v -s %s -c %s '%( os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
	  self.rungmx( equilibrate, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	  self.files.increment_gro()    # must increment filename for any new gmx file 
  
	  # Just to be safe, make a new index for the separate temperature groups "protein" and "sol"
	  self.make_ndx( self.files.grofile, self.files.next_ndx() )
	  self.files.increment_ndx()    # must increment filename for any new gmx file 
		  
	  # setup files for simulation 
	  ### make a tpr file with grompp
	  grompp = '%s/grompp -f %s -c %s -o %s -p %s -n %s'%(os.environ['GMXPATH'], self.files.mdpfile_Simulation, self.files.grofile, self.files.next_tpr(), self.files.topfile, self.files.ndxfile)
	  self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	  self.files.increment_tpr()    # must increment filename for any new gmx file 
		  
  
      elif self.protocol == 'default':
	
	  # run minimimzation once more with the new ions
	  ### make a tpr file with grompp
	  grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.mdpfile_Minimization, self.files.grofile, self.files.next_tpr(), self.files.topfile)
	  self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	  self.files.increment_tpr()    # must increment filename for any new gmx file 
	  ### run minimization
	  minimize = '%s/mdrun -v -s %s -c %s '%( os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
	  self.rungmx( minimize, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	  self.files.increment_gro()    # must increment filename for any new gmx file 
	  
	  # setup files for equilibration 
	  ### make a tpr file with grompp
	  grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.mdpfile_Equilibration, self.files.grofile, self.files.next_tpr(), self.files.topfile)
	  self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	  self.files.increment_tpr()    # must increment filename for any new gmx file 
		  
  
      else:
	
	  print 'The preparation protocol %s is not supported! Exiting....'%self.protocol
	  sys.exit(1)
  

	  

      # copy the final files to the final output directory
      self.copyFilesToOutputDir()
	  
      # write logfiles
      self.status.writeLog( os.path.join(self.finalOutputDir,'output.log') )
      self.status.writeHistory( os.path.join(self.finalOutputDir,'history.log') )
	  
	  
      # Clean up the working directory
      if (cleanup):
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


      return


    


  def getSolventGroupFromGenion(self):
      """Gets the number of the 'SOL' group from the text output of genion."""
     
      [np, nn, nwaters] = self.counterions()
      genion = 'echo q | genion -s %s -o %s -pname %s -np %d -pq %d -nname %s -nn %d -nq %d -g genion.log'%(self.files.tprfile, self.files.next_gro(), self.setup.positiveIonName, np, self.setup.positiveIonCharge, self.setup.negativeIonName, nn, self.setup.negativeIonCharge)
      genion_lines = commands.getoutput(genion).split('\n')
      groupstr = None
      for line in genion_lines:
          if line[0:5] == 'Group':
              fields = line.replace(')',' ').replace('(',' ').split()
              if fields[2] == 'SOL':
                  groupstr = fields[1]      
      if groupstr == None:
          print 'Cannot find atomgroup number for SOL from genion output:'
          for line in genion_lines:
             print '>>>>', line
          raise IOError
      return int(groupstr) 


  def implicitSolvationPreparationSteps(self, implicitOptions=None, useTable=None):
      """Do the steps needed for implicit solvation runs"""
      
      # Find the shortest box vector
      grostruct = GromacsStructure()
      grostruct.load(self.files.grofile)
      boxstring = grostruct.getboxstring()
      boxvector = []
      for str in boxstring.strip().split():
          num = float(str)
	  if num > 0.0:
	    if len(boxvector) < 3:
	      boxvector.append(num)
      shortDimension = min(boxvector)
      cutoffLength = shortDimension/2.0 - 0.05
      
      ##### MINIMIZATION:  #####
      # Make rlist, rcoulomb, and rvdw smaller than half the box length.
      # NOTE:  Edgar warns: Update: Looks like ns_type simple forbids twin range cutoffs as well, so make sure rlist=rcoulomb=rvdw.
      min1_mdpfile = MdpFile(self.files.mdpfile_Minimization)
      min1_mdpfile.setParameter('rlist','%3.2f'%(cutoffLength))
      min1_mdpfile.setParameter('rcoulomb','%3.2f'%(cutoffLength))
      min1_mdpfile.setParameter('rvdw','%3.2f'%(cutoffLength))
      if useTable != None:
	  # if we have a periodic box, override the shift and cutoff parameters of the table to conform to our "box"???
          if min1_mdpfile.params['pbc'] == 'xyz':
	      useTable.rc = cutoffLength + 0.5
	      if useTable.name == 'Switch_Shift':
		    useTable.rs = cutoffLength
	  min1_mdpfile.use_table(useTable)
	  min1_mdpfile.write_table(self.files.xvgfile)
      if implicitOptions != None:
	  # All StillGB, AGBNP, EdGB etc. should  be turned off for minimization
	  implicitOptions.use_GB = False
	  implicitOptions.use_StillGB = False
	  implicitOptions.use_AGBNP = False
          # Otheriwse, make sure all the userints and userreals are what we want
          min1_mdpfile.setParametersFromDict( implicitOptions.mdpParams ) 
      min1_mdpfile.write(self.files.mdpfile)
      
      # minimize the system            ### make a tpr file with grompp
      grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.mdpfile, self.files.grofile, self.files.tprfile, self.files.topfile)
      self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
      ### run minimization
      minimize = '%s/mdrun -v -s %s -c %s '%(os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
      self.rungmx( minimize, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
      self.files.increment_gro()    # must increment filename for any new gmx file 
      


      # #### EQUILIBRATION #### #
      # Make rlist, rcoulomb, and rvdw smaller than half the box length
      equil1_mdpfile = MdpFile(self.files.mdpfile_Equilibration)
      equil1_mdpfile.setParameter('rlist','%3.2f'%(cutoffLength))
      equil1_mdpfile.setParameter('rcoulomb','%3.2f'%(cutoffLength))
      equil1_mdpfile.setParameter('rvdw','%3.2f'%(cutoffLength))
      if useTable != None:
	  # if we have a periodic box, override the shift and cutoff parameters of the table to conform to our "box"???
          if equil1_mdpfile.params['pbc'] == 'xyz':
	      useTable.rc = cutoffLength + 0.5
	      if useTable.name == 'Switch_Shift':
		    useTable.rs = cutoffLength
	  equil1_mdpfile.use_table(useTable)
	  equil1_mdpfile.write_table(self.files.xvgfile)
      # turn back on the implicit solvation for equilibration and production
      if implicitOptions != None:
          if self.protocol.count('implicitPS3_GB') > 0:
              implicitOptions.use_GB = True
          if self.protocol.count('implicitPS3_StillGB') > 0:
              implicitOptions.use_StillGB = True
          if self.protocol.count('implicitPS3_AGBNP') > 0:
              implicitOptions.use_AGBNP = True
          # Otheriwse, make sure all the userints and userreals are what we want
      equil1_mdpfile.setParametersFromDict( implicitOptions.mdpParams )  # make sure all the userints and userreals are what we want
      self.files.increment_mdp()
      equil1_mdpfile.write(self.files.mdpfile)
      
      # setup files for equilibration 
      ### make a tpr file with grompp
      grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.mdpfile, self.files.grofile, self.files.next_tpr(), self.files.topfile)
      self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
      self.files.increment_tpr()    # must increment filename for any new gmx file 
      ### run equilibrate
      equilibrate = '%s/mdrun -v -s %s -c %s '%( os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
      self.rungmx( equilibrate, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
      self.files.increment_gro()    # must increment filename for any new gmx file 


      # #### PRODUCTION ####
      # FOR PRODUCTION:  Make rlist, rcoulomb, and rvdw smaller than half the box length

      # Create an MdpFile object -- we'll write out the mdpfile at the end.
      prod1_mdpfile = MdpFile(self.files.mdpfile_Simulation)
      prod1_mdpfile.setParameter('rlist','%3.2f'%(cutoffLength))
      prod1_mdpfile.setParameter('rcoulomb','%3.2f'%(cutoffLength))
      prod1_mdpfile.setParameter('rvdw','%3.2f'%(cutoffLength))
      if useTable != None:
	  # if we have a periodic box, override the shift and cutoff parameters of the table to conform to our "box"???
          if prod1_mdpfile.params['pbc'] == 'xyz':
	      useTable.rc = cutoffLength + 0.5
	      if useTable.name == 'Switch_Shift':
		    useTable.rs = cutoffLength
	  prod1_mdpfile.use_table(useTable)
	  prod1_mdpfile.write_table(self.files.xvgfile)
      # Create an index with the separate temperature groups "protein" and "sol" 
      self.make_ndx( self.files.grofile, self.files.next_ndx() )
      self.files.increment_ndx()    # must increment filename for any new gmx file 

      # build the mdpParams for the extra PS3 Visualization groups inside the implicitOptions() object.
      # NOTE:  This will also create a new, incremented index file with the visualiztion groups defined
      if implicitOptions != None:
          self.make_PS3_visualization_groups( self.files.grofile, self.files.next_ndx(), implicitOptions)
          self.files.increment_ndx()    # must increment filename for any new gmx file 
  	  prod1_mdpfile.setParametersFromDict( implicitOptions.mdpParams )  # make sure all the userints and userreals are what we want
      self.files.increment_mdp()
      prod1_mdpfile.write(self.files.mdpfile)
      

      # setup files for simulation 
      ### make a tpr file with grompp
      grompp = '%s/grompp -f %s -c %s -o %s -p %s -n %s'%(os.environ['GMXPATH'], self.files.mdpfile, self.files.grofile, self.files.next_tpr(), self.files.topfile, self.files.ndxfile)
      self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
      self.files.increment_tpr()    # must increment filename for any new gmx file 


  def postSolvationPreparationSteps(self):
    """Perform the preparation steps (addig ions,  minimization equil, and production prep) after the box is solvated."""

    # calculate how many ions to add to the box
    [np, nn, nwaters] = self.counterions()	 

    # get the solvent atomgroup number from genion
    solgroup = self.getSolventGroupFromGenion()

    # add ions	   
    genion = 'echo %d | %s/genion -s %s -o %s -pname %s -np %d -pq %d -nname %s -nn %d -nq %d -g genion.log'%(solgroup, os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro(), self.setup.positiveIonName, np, self.setup.positiveIonCharge, self.setup.negativeIonName, nn, self.setup.negativeIonCharge)
    self.rungmx( genion, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
    self.files.increment_gro()    # must increment filename for any new gmx file 

    ### generate a new topolgy file for the ion-ated grofile  
    pdb2gmx = 'echo %s | %s/pdb2gmx -f %s -o %s -p %s -ignh -H14'%(self.forcefieldCodes[self.useff], os.environ['GMXPATH'], self.files.grofile, self.files.next_gro(), self.files.next_top())
    self.rungmx( pdb2gmx, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
    self.files.increment_gro()    # must increment filename for any new gmx file 
    self.files.increment_top()    # must increment filename for any new gmx file 
    self.useTIP3P( self.files.topfile )
	    
    if (self.protocol == 'racecar2') | (self.protocol == 'racecar3') | (self.protocol == 'quicktest'):

	# run minimimzation once more with the new ions
	### make a tpr file with grompp
	grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.mdpfile_Minimization, self.files.grofile, self.files.next_tpr(), self.files.topfile)
	self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	self.files.increment_tpr()    # must increment filename for any new gmx file 
	### run minimization
	minimize = '%s/mdrun -v -s %s -c %s '%( os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
	self.rungmx( minimize, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	self.files.increment_gro()    # must increment filename for any new gmx file 

	# run equilibration phase 1: just water with frozen protein.  
	### make a tpr file with grompp
	grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.mdpfile_Equilibration, self.files.grofile, self.files.next_tpr(), self.files.topfile)
	self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	self.files.increment_tpr()    # must increment filename for any new gmx file 
	### run equilibration
	equilibrate = '%s/mdrun -v -s %s -c %s '%(os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
	self.rungmx( equilibrate, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	self.files.increment_gro()    # must increment filename for any new gmx file 
	
	# Make an index for the separately temperature groups "protein" and "sol"
	self.make_ndx( self.files.grofile, self.files.ndxfile )
	
	# run equilibration phase 2: protein and sol.  
	if (self.protocol == 'racecar2'):
            self.files.set_mdpEquilibration('racecar2_equil2.mdp')
	elif (self.protocol == 'racecar3'):
            self.files.set_mdpEquilibration('racecar3_equil2.mdp')
        elif (self.protocol == 'quicktest'): 
            self.files.set_mdpEquilibration('quicktest_equil2.mdp')
        else:
            raise "Protocol Error!"

	### make a tpr file with grompp
	grompp = '%s/grompp -f %s -c %s -o %s -p %s -n %s'%(os.environ['GMXPATH'], self.files.mdpfile_Equilibration, self.files.grofile, self.files.next_tpr(), self.files.topfile, self.files.ndxfile)
	self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	self.files.increment_tpr()    # must increment filename for any new gmx file 
	### run equilibrate
	equilibrate = '%s/mdrun -v -s %s -c %s '%( os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
	self.rungmx( equilibrate, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	self.files.increment_gro()    # must increment filename for any new gmx file 

	# Just to be safe, make a new index for the separate temperature groups "protein" and "sol"
	self.make_ndx( self.files.grofile, self.files.next_ndx() )
	self.files.increment_ndx()    # must increment filename for any new gmx file 
		
	# setup files for simulation 
	### make a tpr file with grompp
	grompp = '%s/grompp -f %s -c %s -o %s -p %s -n %s'%(os.environ['GMXPATH'], self.files.mdpfile_Simulation, self.files.grofile, self.files.next_tpr(), self.files.topfile, self.files.ndxfile)
	self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	self.files.increment_tpr()    # must increment filename for any new gmx file 
		

    elif self.protocol == 'default':
      
	# run minimimzation once more with the new ions
	### make a tpr file with grompp
	grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.mdpfile_Minimization, self.files.grofile, self.files.next_tpr(), self.files.topfile)
	self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	self.files.increment_tpr()    # must increment filename for any new gmx file 
	### run minimization
	minimize = '%s/mdrun -v -s %s -c %s '%( os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
	self.rungmx( minimize, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	self.files.increment_gro()    # must increment filename for any new gmx file 
	
	# setup files for equilibration 
	### make a tpr file with grompp
	grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.mdpfile_Equilibration, self.files.grofile, self.files.next_tpr(), self.files.topfile)
	self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors  )
	self.files.increment_tpr()    # must increment filename for any new gmx file 
		

    else:
      
	print 'The preparation protocol %s is not supported! Exiting....'%self.protocol
	sys.exit(1)




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
	fout.close()
	os.chmod(equilscript, 0775)
	if self.verbose==True:
	    print cmdout




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
      
      if checkForFatalErrors:
	if output.count('Fatal error') > 0:
	  print '*** There was a fatal error in the following command: ***'
	  print '*********************************************************'
	  print cmd
	  print '*********************************************************'
          print 'Exiting...'
	  sys.exit(1)
	  
    # os.chdir( origdir )
    return output
  

  def useTIP3P(self, topfile):
    """Change a line in the topology file to use ffamber_tip3p."""

    print '>> Using TIP3P water - changing topology file %s...'%topfile
    fin = open( topfile, 'r')
    lines = fin.readlines()
    fin.close()
	
    fout = open( topfile, 'w')
    for line in lines:
      if line[0:8] == '#include':
        if line.count('#include "%s.itp"'%self.useff) > 0:
            fout.write('#include "%s.itp"\n#include "ffamber_tip3p.itp"\n'%self.useff)
        elif line.count('#include "flexspc.itp"') > 0:
           continue
        elif line.count('#include "spc.itp"') > 0:
           continue
        else:
            ifile = line.split()[1].replace('"','')
            if os.path.isfile(ifile):
              finclude = open( ifile , 'r')
            else:
              finclude = open(os.path.join(self.files.GMXLIB,ifile),'r')
            includelines =  finclude.readlines()
            for iline in includelines:
                fout.write( iline )
            finclude.close()
      else:
        fout.write(line)
    fout.close()


  def counterions(self):
    
    # Find out the number of water molecules there are in the box
    fin = open(self.files.topfile, 'r')
    lines = fin.readlines()
    fin.close()
    foundit = False
    nwaters = None
    while ( (len(lines) > 0) & (foundit==False)):
       fields = lines.pop().split()
       print 'fields', fields       
       if (fields[0] == 'SOL') or (fields[0] == 'HOH'):
	   nwaters = int(fields[1])
	   foundit = True

    # Calculate the number of ions that should be in solution based on the concentration
    numberfraction = self.ioncalc.molarityToFraction( self.setup.saltconc )
    print
    print 'Ratio of salt molecules to solvent molecules:', numberfraction
   
    ######################################################################################## 
    # Let f be the fraction of salt molecules in a solution of total solvent molecules N.
    # When the ions dissociate, there are a total of M*N*f (cations+anions), where
    #
    #     M = (self.setup.positiveStoichiometry + self.setup.negativeStoichiometry)
    #
    # When we add cations an anions with the genion* program, each *ion* displaces a
    # solvent molecule! How many total salt molecules x should we introduce to maintain
    # the correct concentration f?  x must satisfy the equation:
    #
    #     x/(N-Mx) = f
    #
    # which has the solution:
    #
    #     x = fN/(1+Mf)
    #
    #######################################################################################

    M = (self.setup.positiveStoichiometry + self.setup.negativeStoichiometry)
    # counterions for every salt molecule

    nsaltmolecules = int(numberfraction*float(nwaters)/(1.+float(M)*numberfraction))
    print 'Number of salt molecules to add:', nsaltmolecules
     
    # Balance the charge with the appropriate number of counterions
    # Our philosophy here will adhere to simple rules:
    #     1) always neutralize the charge in the box
    #     2) always decrease the total numbers of ions, if possible
    #     3) there will be divalent anions, only -1 charges.

    # Here's the formula:

    # Start with the number of counterions predicted from salt concentration
    np = self.setup.positiveStoichiometry*nsaltmolecules
    nn = self.setup.negativeStoichiometry*nsaltmolecules

    print
    print "Finding the best numbers of addded positive and negative ions..."

    chargesOk=False 
    total_ioncharge = (nn*self.setup.negativeIonCharge + np*self.setup.positiveIonCharge)                   
    while (chargesOk==False):
        if (total_ioncharge + self.totalChargeBeforeIons) > 0:
            # are there enough pos ions to remove?
            if np > 0:
                np=np-1
            else:
                nn=nn+1
        elif (total_ioncharge + self.totalChargeBeforeIons) < 0:
            # are there enough neg ions to remove?
            if nn > 0:
                nn=nn-1
            else:
                np=np+1
        else:
            chargesOk=True
	
        total_ioncharge = (nn*self.setup.negativeIonCharge + np*self.setup.positiveIonCharge)                   
	if total_ioncharge == self.totalChargeBeforeIons:
	    chargesOk=True

    print '...Done.'
    print 'Total charge of system before added ions:', self.totalChargeBeforeIons
    print 'Total charge of added ions:', total_ioncharge
    print '    ',nn,self.setup.negativeIonName,np,self.setup.positiveIonName
    
    return [np, nn, (nwaters-np-nn)]


  def make_ndx(self, grofile, ndxfile, groups=['protein','ions','sol', 'bath'], implicitOptions=None):
    """Make an *.ndx file with atom groups corresponding to a list of commmon-sense names
    given in the groups=[] list.  Supported group names so far:
    
        protein            All protein atoms
	sol                All solvent atoms.  (protein + sol + ions) shouuld equal (system) !
        ions               All ions
	bath               (sol + ions)
	
    Note: This subroutine assumes it's being called from the current working directory
    """
 
    [protein_residues, ion_residues, solvent_residues] = self.find_residue_groups(grofile)
    bath_residues = ion_residues + solvent_residues

    if self.debug:
      print 'protein_residues', protein_residues
      print 'ion_residues', ion_residues
      print 'solvent_residues', solvent_residues
      print 'bath_residues', bath_residues
    
    # write a atomSelect-style selections file:
    ### selections file consists of blank lines, lines beginning with these keywords, and      
    ###  titles. Keywords: 'all', 'not', 'atomname', 'residues'
    ###  titles begin and end with "["/"]"
    ###  example: 
    ###  [ helix1-CA ]
    ###  residues 4-10
    ###  only CA
    ###
    ### this little file will produce an index group called 'helix1-CA' which consists of atoms
    ### in resiudes 4 through 10 which are called 'CA'
  
    selectionsFile = 'selections'
    fout = open(selectionsFile, 'w')
    
    if groups.count('protein') > 0:
      if len(protein_residues) > 0:
        fout.write( '[ protein ]\nresidues ' )
	for res in protein_residues[0:-1]:
	    fout.write('%d, '%res)  
	fout.write('%d\n\n'%protein_residues[-1])
	
    if groups.count('sol') > 0:     
      if len(solvent_residues) > 0:
        fout.write( '[ sol ]\nresidues ' )
	for res in solvent_residues[0:-1]:
	    fout.write('%d, '%res)  
	fout.write('%d\n\n'%solvent_residues[-1])
	
    if groups.count('ions') > 0:     
      if len(ion_residues) > 0:
        fout.write( '[ ions ]\nresidues ' )
	for res in ion_residues[0:-1]:
	    fout.write('%d, '%res)  
	fout.write('%d\n\n'%ion_residues[-1])

    if groups.count('bath') > 0:     
      if len(bath_residues) > 0:
        fout.write( '[ bath ]\nresidues ' )
	for res in bath_residues[0:-1]:
	    fout.write('%d, '%res)  
	fout.write('%d\n\n'%bath_residues[-1])

        
    fout.close()
   
    indexGroups = getSelections( selectionsFile, grofile )
    if indexGroups == None:
      # TO DO: *should* raise an exception
      print 'There was a problem with selections index atoms with atomSelect.py!'
      print 'Check the selections file: %s',os.path.abspath(selectionsFile)
      print 'Exiting...'
      sys.exit(1)
      
    fndx = open(ndxfile,'w')
    for grp in indexGroups:
      fndx.write( grp.getIndexGroupText() )



  def find_residue_groups(self, grofile):
    """From a grofile, read through the residue names and return three lists of numbers:
      [protein_residues, ion_residues, solvent_residues]."""

    # First, scan through the grofile for the residue numbers that correspond to the protein and solvent
    fin = open(grofile, 'r')
    lines = fin.readlines()
    lines.pop(0)    # skip title
    lines.pop(0)    # skip number of atoms
    lines.pop()     # skip box size lines
    protein_residues = []
    solvent_residues = []
    ion_residues = []
    res = -1
    lastRes = -1
    for line in lines:    # these lines should all have atoms on them
      # for first line, get res number, after that just increment
      if res == -1:
        res = int(line[0:5])
        lastRes = res

      # if current res# from reading line != lastRes then move to next res
      if lastRes != int(line[0:5]):
        res += 1
        lastRes = int(line[0:5])

      if (line[5:8] == 'SOL') or (line[5:8] == 'HOH'):
          if solvent_residues.count(res) == 0:
            solvent_residues.append(res)
      elif (line[5:8] == 'Na ') or (line[5:8] == 'Cl '):
          if ion_residues.count(res) == 0:
            ion_residues.append(res)
      else:
          if protein_residues.count(res) == 0:
            protein_residues.append(res)

    return [protein_residues, ion_residues, solvent_residues]

    
  def make_PS3_visualization_groups(self, grofile, ndxfile, implicitOptions, groups=['protein','ions','sol']):     
    """This routine does two things:
       
    1) Add user1_grps and user2_grps atom groups to the index file.  This
       is for the visualiztion options on the PS3.  For details, see make_ndx().
    2) Adds user1_grps = 'sticks'  (e.g.) and user2_grps = 'space-filled' (e.g.)
       to the mdp file

    """

    [protein_residues, ion_residues, solvent_residues] = self.find_residue_groups(grofile)

    selectionsFile = 'selections'
    fout = open(selectionsFile, 'w')
 
    # 0) Add protein, ion and solvent groups as usual....
    if groups.count('protein') > 0:
      if len(protein_residues) > 0:
        fout.write( '[ protein ]\nresidues ' )
        for res in protein_residues[0:-1]:
            fout.write('%d, '%res)
        fout.write('%d\n\n'%protein_residues[-1])

    if groups.count('sol') > 0:
      if len(solvent_residues) > 0:
        fout.write( '[ sol ]\nresidues ' )
        for res in solvent_residues[0:-1]:
            fout.write('%d, '%res)
        fout.write('%d\n\n'%solvent_residues[-1])

    if groups.count('ions') > 0:
      if len(ion_residues) > 0:
        fout.write( '[ ions ]\nresidues ' )
        for res in ion_residues[0:-1]:
            fout.write('%d, '%res)
        fout.write('%d\n\n'%ion_residues[-1])
 

    # 1) Add PS3 Advanced visualization groups to the index file -- for use only with implicitOptions
    if implicitOptions.user1_residues != None:
        fout.write( '[ %s ]\nresidues '%implicitOptions.user1_grps )
        for res in implicitOptions.user1_residues[0:-1]:
            fout.write('%d, '%res)
        fout.write('%d\n\n'%implicitOptions.user1_residues[-1])
    if implicitOptions.user2_residues != None:
        fout.write( '[ %s ]\nresidues '%implicitOptions.user2_grps )
        for res in implicitOptions.user2_residues[0:-1]:
            fout.write('%d, '%res)
        fout.write('%d\n\n'%implicitOptions.user2_residues[-1])

    fout.close()

    indexGroups = getSelections( selectionsFile, grofile )
    if indexGroups == None:
      # TO DO: *should* raise an exception
      print 'There was a problem with selections index atoms with atomSelect.py!'
      print 'Check the selections file: %s',os.path.abspath(selectionsFile)
      print 'Exiting...'
      sys.exit(1)

    fndx = open(ndxfile,'w')
    for grp in indexGroups:
      fndx.write( grp.getIndexGroupText() )

    # 2) Add the user1_grps and user2_grps to the mdp file
    print 'in make_PS3...  before implicitOptions.buildMdpParams:  implicitOptions.mdpParams =',implicitOptions.mdpParams
    implicitOptions.buildMdpParams(PS3_Visualization=True)

    return
 
