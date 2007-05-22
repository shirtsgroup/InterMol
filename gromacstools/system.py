# gromacstools/system.py
#
# REQUIREMENTS
#
#    1) GROMACS ( http://www.gromacs.org/ ) must be installed on your machine, 
#    and the following Environment variables must be defined (in your .bash_profile, e.g.):
# 
#        GMXLIB             the pathname of the gmx parameter files.  
#        GMXPATH            the pathname of the Gromacs exectuables  (This also needs to be in your PATH)
#        MMTOOLSPATH        the pathname of the mmtools library
#
#    2) The AMBER ports for GROMACS ( http://folding.stanford.edu/ffamber/ ) by Sorin, Park
#    and Pande must be installed, with the *.itp, *.rtp, etc., files in the $GMXLIB directory (.../share/top )   

# ----------------------------------------------------------------------
# TODO:
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
import ioncalc


# ----------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# CLASSES
# ----------------------------------------------------------------------


class GromacsSystem:

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

    # process infile name
    self.infile = infile
    self.infile_basename = os.path.basename(self.infile) 
    self.infile_suffix = self.infile.split('.').pop()    

    # process outdir
    if outdir == None:
        self.outdir = os.curdir 
    else:
        self.outdir = outdir
        print "GROMACS *.tpr and *.gro for simulation will be written to directory %s..." % self.outdir
    
    # process workdir
    if workdir == None:
        self.workdir = tempfile.mkdtemp();
	print "Building GROMACS project in temporary directory %s..." % self.workdir
    else:
        self.workdir = workdir
        print "Building GROMACS project in directory %s..." % self.workdir
    output = commands.getoutput('cp %s %s'%(self.infile, os.path.join(self.workdir,self.infile_basename) ))
    print output
 
    # process other options
    self.useff = useff
    self.version = version
    self.verbose = verbose

    # check to see if the necessary environment variables are defined   
    try:
        self.GMXPATH     = os.environ['GMXPATH']    
	self.GMXLIB      = os.environ['GMXLIB']
        self.MMTOOLSPATH = os.environ['MMTOOLSPATH']
    except KeyError:
	print """Cannot find one or more of the following shell environment variables

    GMXLIB             the pathname of the gmx parameter files.  
    GMXPATH            the pathname of the Gromacs exectuables
    MMTOOLSPATH        the pathname of the mmtools library\n"""
	
	for env in os.environ.keys():
	    print '%-16s\t\t%s'%(env,os.environ[env])
	print '-------------\nExiting....'
	sys.exit(1)

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
  	
    # initialize default mdpfiles
    self.mdpfile_Minimization  = os.path.join(self.MMTOOLSPATH,'gromacstools/mdp/minimize.mdp')
    self.mdpfile_Equilibration = os.path.join(self.MMTOOLSPATH,'gromacstools/mdp/equilibrate.mdp')
    self.mdpfile_Simulation    = os.path.join(self.MMTOOLSPATH,'gromacstools/mdp/simulate.mdp')
    
    # initialize status object
    self.status = GromacsSystemStatus()

    # initialize setup object
    self.setup = GromacsSystemSetup()
    
    # to keep track of the total charge in the system
    # NOTE: this will be calculated from the pdb2gmx output text
    self.totalChargeBeforeIons = 0.0    
    
    # initialize ionCalculator
    # This is needed for calculating numbers of counterions from salt concentrations
    import ioncalc
    self.ioncalc = ioncalc.ionCalculator(self.useff)


    # initialize GROMACS file names
    # NOTE: These names do not correspond to files yet, but will be used
    #       by self.prepare() as a consistent naming scheme when buiilding the tpr 
    if self.infile[-4:] == '.pdb':
        self.grofile_prep = os.path.join(self.workdir, self.infile_basename.replace('.pdb','.gro'))
    elif self.infile[-4:] == '.gro':
        self.grofile_prep = os.path.join(self.workdir, self.infile_basename) 
    else:
        print 'Error:  infile must be either *.pdb or *.gro file'
	sys.exit(1)
    self.topfile_prep = self.grofile_prep.replace('.gro','.top')
    
    self.grofile_box= self.grofile_prep.replace('.gro','_box.gro')
    self.grofile_solvated = self.grofile_prep.replace('.gro','_solvated.gro')
    self.tprfile_solvated = self.grofile_prep.replace('.gro','_solvated.tpr') 
    self.grofile_solvated_afterem = self.grofile_solvated.replace('.gro','_afterem.gro')

    self.grofile_ions = self.grofile_prep.replace('.gro','_ions.gro')
    self.grofile_ions2 = self.grofile_prep.replace('.gro','_ions2.gro')
    self.topfile_ions = self.topfile_prep.replace('.top','_ions.top')

    self.tprfile_ions = self.grofile_prep.replace('.gro','_ions.tpr')
    self.grofile_ions_afterem = self.grofile_ions.replace('.gro','_afterem.gro')
    
    self.tprfile_equil = self.tprfile_ions.replace('_ions.tpr','_equil.tpr')
    self.grofile_equil = self.tprfile_equil.replace('.tpr','.gro')
    
    
    if (verbose):
      self.printGromacsFiles()

 
  def printGromacsFiles(self):
    """A printout of the defined GROMACS filenames"""
    
    print
    print 'Files used in the GROMACS pre-processing:'
    self.printTwoColumnStrings('self.grofile_prep', self.grofile_prep)
    self.printTwoColumnStrings('self.topfile_prep', self.topfile_prep)

    self.printTwoColumnStrings('self.grofile_box', self.grofile_box)
    self.printTwoColumnStrings('self.grofile_solvated', self.grofile_solvated )
    self.printTwoColumnStrings('self.tprfile_solvated', self.tprfile_solvated )
    self.printTwoColumnStrings('self.grofile_solvated_afterem', self.grofile_solvated_afterem)

    self.printTwoColumnStrings('self.grofile_ions', self.grofile_ions ) 
    self.printTwoColumnStrings('self.topfile_ions', self.topfile_ions )


  def printTwoColumnStrings(self,s1,s2):
    """Prints two strings in a fixed-width column."""
    print '%-36s %s'%(s1,s2) 
    return



  def getForcefields(self):
    """Get the forcefields and codes from FF.dat; produce error if this is not found."""

    try:
      fin = open(os.path.join(self.GMXLIB,'FF.dat'),'r')
    except IOError:
      print "File FF.dat cannot be found in GMXLIB=%s   ..... \nExiting....\n"%self.GMXLIB
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

    print 'Available forcefields in $GMXLIB =',os.path.join(self.GMXLIB,'FF.dat')
    for key,value in self.forcefields.items():
      print '%-16s%s'%(key,value)
    
    

  def prepare(self, outname=None, outdir=None, cleanup=True, verbose=False):
    
    """Prepare an equilibration simulation for this Gromacs system, namely:
    
    1.  pdb2gmx
    2.  edit topology file (if necessary) to use tip3p water
    3.  solvate the box
    4.  minimize the box
    5.  add ions to the solvated box
    6.  minimize once again
    7.  prepare a molecular dynamics simulation for equilibration...

    OPTIONAL ARGUMENTS
    outname             user-specified name for <outname>.tpr and .gro. Default: "out"
    outdir              save output files to directory outdir
    cleanup             if cleanup=True, delete the temporary pre-processing directory        

    OUTPUT

    Outputs the following files in the current working directory:

        out.tpr		GROMACS tpr file
        out.gro		GROMACS gro (structure) file
        equilibrate*	executable script to call mdrun 
        history.log    logfile containing all the commands used to prepare the simualtion 
        output.log      logfile containg the output of the pre-processing programs

"""
    # do work in the specified working directory
    if outdir == None:
      cwd = os.path.abspath(os.curdir)  ### remember our original location
      print 'cwd', cwd
    else:
      cwd = outdir
      
    if outname == None:
      outname = 'out'
      
    self.verbose=verbose  

    os.chdir( self.workdir )
    print 'Preparing simulation in directory %s...'%self.workdir
    		
    # prepare the simualtion 	
    if self.version == '3.1':	

        # pdb2gmx, and useTIP3P
        pdb2gmx = 'echo %s | pdb2gmx -f %s -o %s -p %s -ignh'%(self.forcefieldCodes[self.useff], self.infile, self.grofile_prep, self.topfile_prep)
	self.rungmx( pdb2gmx )
        self.useTIP3P( self.topfile_prep )
	
	pdb2gmxlines = self.status.loglines[-1][1].split('\n')
	for line in pdb2gmxlines:
	  if line[0:12] == 'Total charge':
	    self.totalChargeBeforeIons = float(line.split()[2])
	    
        # make a cubic (periodic boundary conditions) box
	editconf = 'editconf -bt %s -f %s -o %s -d %s'%(self.setup.boxType, self.grofile_prep, self.grofile_box, self.setup.boxSoluteDistance )
	self.rungmx( editconf )

        # solvate the box 
	editconf = 'genbox -cp %s -cs ffamber_tip3p.gro -o %s -p %s'%(self.grofile_box, self.grofile_solvated, self.topfile_prep)
	self.rungmx( editconf )

        # minimize the system
        ### make a tpr file with grompp
        grompp = 'grompp -f %s -c %s -o %s -p %s -check14'%(self.mdpfile_Minimization, self.grofile_solvated, self.tprfile_solvated, self.topfile_prep)
	self.rungmx( grompp )	
        ### run minimization
        minimize = 'mdrun -v -s %s -c %s '%( self.tprfile_solvated, self.grofile_solvated_afterem )
	self.rungmx( minimize )
	
        # add ions to box
	### Run genion
	[np, nn, nwaters] = self.counterions()	  
	genion = 'echo 12 | genion -s %s -o %s -pname %s -np %d -pq %d -nname %s -nn %d -nq %d -g genion.log'%(self.tprfile_solvated, self.grofile_ions, self.setup.positiveIonName, np, self.setup.positiveIonCharge, self.setup.negativeIonName, nn, self.setup.negativeIonCharge)
	self.rungmx( genion )
        ### generate a new topolgy file for the ion-ated grofile  
        pdb2gmx = 'echo %s | pdb2gmx -f %s -o %s -p %s -ignh'%(self.forcefieldCodes[self.useff], self.grofile_ions, self.grofile_ions2, self.topfile_ions)
	self.rungmx( pdb2gmx )
        self.useTIP3P( self.topfile_ions )
		
	
	# run minimimzation once more with the new ions
        ### make a tpr file with grompp
        grompp = 'grompp -f %s -c %s -o %s -p %s -check14'%(self.mdpfile_Minimization, self.grofile_ions2, self.tprfile_ions, self.topfile_ions)
	self.rungmx( grompp )
        ### run minimization
        minimize = 'mdrun -v -s %s -c %s '%( self.tprfile_ions, self.grofile_ions_afterem )
	self.rungmx( minimize )
	
	# setup files for equilibration 
        ### make a tpr file with grompp
        grompp = 'grompp -f %s -c %s -o %s -p %s -check14'%(self.mdpfile_Equilibration, self.grofile_ions_afterem, self.tprfile_equil, self.topfile_ions)
	self.rungmx( grompp )
 		
        # copy the necessary files to start an MD run back to the original curdir
	
	### the GRO file
	out_grofile = os.path.join(outdir,outname+'.gro')
	cmdout = commands.getoutput( 'cp %s %s'%(self.grofile_ions_afterem, out_grofile) )
	if self.verbose==True:
	    print cmdout
	    
	### the TPR file
	out_tprfile = os.path.join(outdir,outname+'.tpr')
	cmdout = commands.getoutput( 'cp %s %s'%(self.tprfile_equil, out_tprfile) )
	if self.verbose==True:
	    print cmdout
	    
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
	# self.rungmx( pdb2gmx )
	
        # make a cubic (periodic boundary conditions) box
	### FILL THIS IN! ###
	
        # solvate the box
	# editconf = 'genbox -cp %s -cs ffamber_tip3p.gro -o %s -p %s'%(self.grofile_box, self.grofile_solvated, self.topfile_prep)
	# self.rungmx( editconf )

        # minimize the system
        ### FILL IN!
	

    else:
        print 'This version of gromacs (',self.version,') is not supported.'
	sys.exit(1)

    os.chdir( cwd )	### change back to our original directory
    return
  
  
  
  

  def rungmx(self, cmd):
    """Execute a gromacs executable on the command line, keeping track of the command and output for the logfile"""
  
    if self.verbose: print '>>',cmd	
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
    fin = open(self.topfile_prep, 'r')
    lines = fin.readlines()
    fin.close()
    foundit = False
    nwaters = None
    while ( (len(lines) > 0) & (foundit==False)):
       fields = lines.pop().split()
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


class GromacsSystemStatus:

  def __init__(self): 
    """A class to keep track of the status of of the simulation system."""

    # logffile lines
    self.loglines = []     # append (cmd, output)
        
    # Benchmarks
    self.groprepDone = False
    self.genboxDone = False
    self.genionDone = False
    self.minimizeDone = False
    self.equilibrateDone = False
    self.simulateDone = False
  
  
  def writeHistory(self, logfile):
    """Writes logfile lines to file in a standard file format"""
   
    fout = open(logfile,'w')
    for line in self.loglines:
      fout.write(line[0]+'\n')
    fout.close()

  def writeLog(self, logfile):
    """Writes logfile lines to file in a standard file format"""

    fout = open(logfile,'w')
    for line in self.loglines:
      fout.write(line[1]+'\n')
    fout.close()


    

  def __repr__(self): 
    """A class to keep trackk of the status of of the simulation system."""
    
    return "wowwiw!"

class GromacsSystemSetup:

  def __init__(self, infile=None): 
    """A class to keep track of parameters important to the system setup."""

    # editconf parms
    self.boxType = 'cubic'            #  -bt   enum   tric  Box type for -box and -d: tric, cubic, dodecahedron or octahedron
    self.boxSoluteDistance = '0.9'    #  -d   real      0  Distance between the solute and the box
    
    # salt conditions
    self.salt = 'NaCl'                #  Supported:  'NaCl' or 'SodiumChloride'
    self.saltconc = 0.050             #  Molar salt concentration
    self.positiveIonName = 'Na'
    self.negativeIonName = 'Cl'
    self.positiveStoichiometry = 1
    self.negativeStoichiometry = 1
   
    # setup from file?  THIS is still a DUMMY function....
    if infile != None:
        self.read(infile)

	
  def setSaltConditions(self, saltname, saltconcentration): 
    """Set ion names and stoichiometry for the chosen salt."""
    if saltname == 'NaCl' or saltname=='sodium chloride':
      self.salt = 'NaCl'                #  Supported:  'NaCl' or 'SodiumChloride'
      self.saltconc = saltconcentration            #  Molar salt concentration
      self.positiveIonName = 'Na'
      self.negativeIonName = 'Cl'
      self.positiveIonCharge = 1
      self.negativeIonCharge = -1
      self.positiveIonStoichiometry = 1
      self.negativIoneStoichiometry = 1

    else:
        print 'salt type', salt, 'not supported!  Try \'NaCl\'.  Exiting.'
	sys.exit(1)
    return None
  
  
  def read(self, filename): 
    """Read parameters in from a GromacsSystemSetup text file."""
    return None
  
  def write(self, filename):
    """Write parameters to a GromacsSystemSetup text file."""
    return None
  

  
  

