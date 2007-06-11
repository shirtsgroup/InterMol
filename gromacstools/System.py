# gromacstools/system.py
#
# REQUIREMENTS
#
#    1) GROMACS ( http://www.gromacs.org/ ) must be installed on your machine, 
#       and the following Environment variables must be defined (in your .bash_profile, e.g.):
# 
#        GMXLIB             the pathname of the gmx parameter files.  
#        GMXPATH            the pathname of the Gromacs exectuables  (This also needs to be in your PATH)
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



# ----------------------------------------------------------------------
# TODO:
# *** MOST IMPORTANT ***
# * Implement Exceptions!!!  One failed protein should cripple the rest of the pipeline.  
# *** IMPORTANT BUT LESS SO ***
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

import sys, os, string, commands
import tempfile

from atomSelect import *
from Filenames import *
from IonCalculator import *
from Setup import *
from Status import *


# ----------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# CLASSES
# ----------------------------------------------------------------------


class System:

  def __init__(self, infile, outdir=None, workdir=None, useff='amber99p',  version='3.1', verbose=True): 
    """A class to initialze, minimize, and equilibrate Gromacs simulations.

    REQUIRED INPUTS
    infile		the stating *.pdb, *.gro, etc. file to initialize the system 

    OPTIONAL INPUTS
    outdir		directory to save final *.trp, *.gro files for an MD simulation 
    workdir		directory to perform gmx preprocessing.  If not specified, work
                        will be done in a unique temporary directory.
    ff			string to specify forcefield (default is 'amber99p') 
    version             gmx version -- supported: '3.1.4' ('3.3.1' COMING SOON)
    
"""


    # initialize status object
    self.status = Status()

    # initialize setup object
    self.setup = Setup()
    
    # process outdir
    if outdir == None:
        self.outdir = os.curdir 
    else:
        self.outdir = outdir
        print "GROMACS *.tpr and *.gro for simulation will be written to directory %s..." % self.outdir
    
    # process workdir
    if workdir == None:
        self.workdir = tempfile.mkdtemp();
    else:
        self.workdir = workdir

    # initializee files object  
    self.files = Filenames(infile, self.workdir)

    # Start Building the GROMACS project
    print "Building GROMACS project in temporary directory %s..." % self.workdir
    output = commands.getoutput('cp %s %s'%(self.files.infile, os.path.join(self.workdir, self.files.infile_basename) ))
    print output

    
    # process other options
    self.useff = useff
    self.version = version
    self.verbose = verbose


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
    
    

  def prepare(self, outname=None, outdir=None, cleanup=True, verbose=False, debug=False, protocol=None, checkForFatalErrors=True, mockrun=False):
    
    """Prepare an equilibration simulation for this Gromacs system, namely:
    
    1.  pdb2gmx
    2.  edit topology file (if necessary) to use tip3p water
    3.  solvate the box
    4.  minimize the box
    5.  add ions to the solvated box
    6.  minimize once again
    7.  prepare a molecular dynamics simulation for equilibration...

    OPTIONAL ARGUMENTS
    
    ARG                 DEFAULT    DESCRIPTION
    outname             'out'      user-specified name for <outname>.tpr and .gro. Default: "out"
    outdir              './'       save output files to directory outdir.  If not specified, outdir is the current directory
    cleanup             True       if True, delete the temporary pre-processing directory, and erase temporary files 
    verbose             False      if True, print extra, more detailed progress statements 
    debug               False      if True, turn on print statements for debugging    
    protocol            'default'  Currently only 'default' and 'racecar2' are supported:
				   * 'default':   minimize.mdp -> equilibrate.mdp -> simulate.mdp    
				   * 'racecar2':  racecar2_min1.mdp -> racecar2_equil1.mdp -> racecar2_equil2.mdp -> racecar2_grompp.mdp
    checkForFatalErrors	True       Exits the program if any of the gmx produce "Fatal error" in the standard output	    
			        
    
    OUTPUT
  
    The following files will be written to the outdir:
      
        out.tpr		GROMACS tpr file
        out.gro		GROMACS gro (structure) file
        equilibrate*	executable script to call mdrun 
        history.log     logfile containing all the commands used to prepare the simualtion 
        output.log      logfile containg the output of the pre-processing programs

"""

    # PARSE INPUT ARGUMENTS
    
    # outdir
    ### do work in the specified working directory
    if outdir == None:
      cwd = os.path.abspath(os.curdir)  ### remember our original location
      print 'cwd', cwd
    else:
      cwd = os.path.abspath(outdir)
    
    # outname    
    ### specify default outfile name as 'out'
    if outname == None:
      outname = 'out'

    self.cleanup = cleanup                             # Erase the temporary files when done
    self.verbose = verbose                             # Print extra, more detailed progress statements 
    self.debug   = debug                               # Turns on print statements for debugging
    if protocol == None:                               # Current only 'default' and 'racecar2' are supported:
      self.protocol = 'default'                        #     'default':   minimize.mdp -> equilibrate.mdp -> simulate.mdp 
                                                       #     'racecar2':  racecar2_min1.mdp -> racecar2_equil1.mdp -> racecar2_equil2.mdp -> racecar2_grompp.mdp  
    self.checkForFatalErrors = checkForFatalErrors
    self.mockrun = mockrun
    
    

    os.chdir( self.workdir )
    print 'Preparing simulation in directory %s...'%self.workdir
    		
    # prepare the simualtion 	
    if self.version == '3.1':	

        # pdb2gmx, and useTIP3P
        pdb2gmx = 'echo %s | pdb2gmx -f %s -o %s -p %s -ignh'%(self.forcefieldCodes[self.useff], self.files.infile, self.files.grofile, self.files.topfile)
	self.rungmx( pdb2gmx, mockrun=self.mockrun)
        self.useTIP3P( self.files.topfile )
	
	pdb2gmxlines = self.status.loglines[-1][1].split('\n')
	for line in pdb2gmxlines:
	  if line[0:12] == 'Total charge':
	    self.totalChargeBeforeIons = float(line.split()[2])
	    
        # make a cubic (periodic boundary conditions) box
	editconf = 'editconf -bt %s -f %s -o %s -d %s'%(self.setup.boxType, self.files.grofile, self.files.next_gro(), self.setup.boxSoluteDistance )
	self.rungmx( editconf, mockrun=self.mockrun  )
	self.files.increment_gro()    # must increment filename for any new gmx file 

        # solvate the box 
	editconf = 'genbox -cp %s -cs ffamber_tip3p.gro -o %s -p %s'%(self.files.grofile, self.files.next_gro(), self.files.topfile)
	self.rungmx( editconf, mockrun=self.mockrun )
	self.files.increment_gro()    # must increment filename for any new gmx file 

        # minimize the system
        ### make a tpr file with grompp
        grompp = 'grompp -f %s -c %s -o %s -p %s '%(self.files.mdpfile_Minimization, self.files.grofile, self.files.tprfile, self.files.topfile)
	self.rungmx( grompp, mockrun=self.mockrun )	
        ### run minimization
        minimize = 'mdrun -v -s %s -c %s '%( self.files.tprfile, self.files.next_gro() )
	self.rungmx( minimize, mockrun=self.mockrun )
	self.files.increment_gro()    # must increment filename for any new gmx file 
	
        # add ions to box
	### Run genion
	[np, nn, nwaters] = self.counterions()	  
	genion = 'echo 12 | genion -s %s -o %s -pname %s -np %d -pq %d -nname %s -nn %d -nq %d -g genion.log'%(self.files.tprfile, self.files.next_gro(), self.setup.positiveIonName, np, self.setup.positiveIonCharge, self.setup.negativeIonName, nn, self.setup.negativeIonCharge)
	self.rungmx( genion, mockrun=self.mockrun )
	self.files.increment_gro()    # must increment filename for any new gmx file 
        ### generate a new topolgy file for the ion-ated grofile  
        pdb2gmx = 'echo %s | pdb2gmx -f %s -o %s -p %s -ignh'%(self.forcefieldCodes[self.useff], self.files.grofile, self.files.next_gro(), self.files.next_top())
	self.rungmx( pdb2gmx, mockrun=self.mockrun )
	self.files.increment_gro()    # must increment filename for any new gmx file 
	self.files.increment_top()    # must increment filename for any new gmx file 
        self.useTIP3P( self.files.topfile )
		
	
	# run minimimzation once more with the new ions
        ### make a tpr file with grompp
        grompp = 'grompp -f %s -c %s -o %s -p %s '%(self.files.mdpfile_Minimization, self.files.grofile, self.files.next_tpr(), self.files.topfile)
	self.rungmx( grompp, mockrun=self.mockrun  )
	self.files.increment_tpr()    # must increment filename for any new gmx file 
        ### run minimization
        minimize = 'mdrun -v -s %s -c %s '%( self.files.tprfile, self.files.next_gro() )
	self.rungmx( minimize, mockrun=self.mockrun  )
	self.files.increment_gro()    # must increment filename for any new gmx file 
	
	# setup files for equilibration 
        ### make a tpr file with grompp
        grompp = 'grompp -f %s -c %s -o %s -p %s '%(self.files.mdpfile_Equilibration, self.files.grofile, self.files.next_tpr(), self.files.topfile)
	self.rungmx( grompp, mockrun=self.mockrun  )
	self.files.increment_tpr()    # must increment filename for any new gmx file 
 		
        # copy the necessary files to start an MD run back to the original curdir
	print
	print 'Copying the *.gro, *.tpr, and equilibration script to outdir=',outdir,'...'
	
	### the GRO file
	out_grofile = os.path.join(outdir,outname+'.gro')
	copycmd = 'cp %s %s'%(self.files.grofile, out_grofile)
	if (self.verbose):
	  print copycmd 
	cmdout = commands.getoutput( copycmd )
	if self.verbose==True:
	    print cmdout
	    
	### the TPR file
	out_tprfile = os.path.join(outdir,outname+'.tpr')
        copycmd = 'cp %s %s'%(self.files.tprfile, out_tprfile)
	if (self.verbose): print copycmd 
        cmdout = commands.getoutput( copycmd )
	if (self.verbose): print cmdout
	    
        ### the mdrun script
	equilibrate = 'mdrun -v -s %s -c %s '%( out_tprfile, out_grofile )
	equilscript = os.path.join(outdir,'equilibrate')
	fout = open(equilscript,'w')
	fout.write(equilibrate+'\n')
	fout.close()
	os.chmod(equilscript, 0775)
	if self.verbose==True:
	    print cmdout


        # write logfiles
	self.status.writeLog( os.path.join(outdir,'output.log') )
        self.status.writeHistory( os.path.join(outdir,'history.log') )
	    
	    
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

    os.chdir( cwd )	### change back to our original directory
    return
  
  
  
  

  def rungmx(self, cmd, mockrun=False):
    """Execute a gromacs executable on the command line, keeping track of the command and output for the logfile"""
  
    print '>>',cmd	
    if mockrun==False:
      output=commands.getoutput(cmd)
      if self.verbose: print output
      self.status.loglines.append( (cmd, output) )    
    
    return output
  

  def useTIP3P(self, topfile):
    """Change a line in the topology file to use ffamber_tip3p."""

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
            finclude = open( (line.split()[1]).replace('"','') , 'r')
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
       if fields[0] == 'SOL':
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


  def genbox(self, genboxOptions = {}):
    """generate a pbc box and Solvate with water.  genboxOptions"""
    return

  def genion(self, genionOptions={}):
    """generate ions in the box"""
    return

  def minimize(self, mdpfile, outdir=None, SetupOnly=False, customParms={}):
    """Minimize the system according to the specified mdpfile.   If outdir is specified,
    perform calculations/write files to outdir."""
    ### RIGHT NOW this is a dummy function... will be added later.
    return

  def equilibrate(self, outdir, mdpfile, SetupOnly=True, customParms={}):
    """Equilibrate the system."""                     
    ### RIGHT NOW this is a dummy function... will be added later.
    return

  def simulate(self, outdir, mdpfile, SetupOnly=True, customParms={}):
    ### RIGHT NOW this is a dummy function... will be added later.
    """Simulate the system."""                     
    return