import sys, os, string, commands, glob, random
import tempfile

from PrepWorker import *
from SystemSetup import *

from mmtools.gromacstools.devel.Helpers.IndexFile import *
from mmtools.gromacstools.devel.Helpers.TopologyFile import *
from mmtools.gromacstools.devel.Helpers.MdpFile import *
from mmtools.gromacstools.devel.Scripts.calculateCounterionsToAdd import *
from mmtools.gromacstools.devel.Scripts.extractNumberSolventMoleculesFromTopology import *
from mmtools.gromacstools.devel.Structure.PDBStructure import *



class PrepWorkerExplicitSolvation(PrepWorker):
    """An object to store info and coordinate a preparation run.
    
    This class inherents methods from PrepWorker:
        def __init__()
        def mdrun(self, mdpfile)
	def buildTpr(self, mdpfile)
        def cleanupWorkingDir()
        def rungmx()
        def copyFilesToOutputDir()
    """
  

    def prepare(self, systemSetup, cleanup=True, verbose=False, debug=False, checkForFatalErrors=True, mockrun=False):
    
	"""Prepare an equilibration simulation for this Gromacs system, namely:
	
	1.  pdb2gmx
	2.  edit topology file (if necessary) to use tip3p water
	3.  solvate the box
	4.  minimize the box
	5.  add ions to the solvated box
	6.  minimize once again
	7.  prepare a molecular dynamics simulation for equilibration...
	
	REQUIRES ARGUMENTS
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


	##################################################
	# step 1.  Create a grofile from the input file (which can be either a *.pdb or *.gro)
	
	if self.system.gmx_version == '3.1.4':	
	    pdb2gmx = 'echo %s | %s/pdb2gmx -f %s -o %s -p %s -ignh -H14'%(self.system.forcefieldCodes[self.system.forcefield], os.environ['GMXPATH'], self.files.infile, self.files.grofile, self.files.topfile)
	elif self.system.gmx_version == '3.3':
	    pdb2gmx = '%s/pdb2gmx -ff %s -f %s -o %s -p %s -ignh'%(os.environ['GMXPATH'], self.system.forcefield, self.files.infile, self.files.grofile, self.files.topfile)
	else:
	    raise "Only system.gmx_version = '3.1.4' or '3.3' are supported"
	self.rungmx( pdb2gmx )
	
	
	##################################################
	# Step 2.  Since this is implicit solvent, remove any references to SPC or TIP3P water in the topology file
	
        self.useTIP3P()
	    
	    
	##################################################
        # Step 3. make a (periodic boundary conditions) box of an appropriate size
	
        if self.system.useAbsBoxSize:
	    editconf = '%s/editconf -bt %s -f %s -o %s -box %s'%(os.environ['GMXPATH'], self.system.boxType, self.files.grofile, self.files.next_gro(), self.system.absBoxSize )
        else:
       	    editconf = '%s/editconf -bt %s -f %s -o %s -d %s'%(os.environ['GMXPATH'], self.system.boxType, self.files.grofile, self.files.next_gro(), self.system.boxSoluteDistance )
	self.rungmx( editconf )
	self.files.increment_gro()    # PrepWorker must increment filename for any new gmx file 


	##################################################
        # Step 4.  Solvate the system
	
	if not self.system.cosolventName == None:
	    waterbox = self.system.cosolventWaterBoxFile
	else:
	    waterbox = self.system.WaterBoxFile
	editconf = '%s/genbox -cp %s -cs %s -o %s -p %s'%(os.environ['GMXPATH'], self.files.grofile, waterbox, self.files.next_gro(), self.files.topfile)
	self.rungmx( editconf )
	self.useTIP3P()
	self.files.increment_gro()    # PrepWorker must increment filename for any new gmx file 
	
	# Minimize
	self.mdrun('minimize.mdp')

	##################################################
        # Step 5. Add counterions to neutralize the charge

        # get the total charge in the system
	totalSystemCharge = self.calcTotalCharge()
	
        # calculate how many ions and water molecules *should* be in the box
	[numberPositiveIons, numberNegativeIons, numWaters] = calculateCounterionsToAdd(self.system, self.files.topfile, totalSystemCharge)
	
	# Add the ions to the grofile and topfile
	self.AddIons(numberPositiveIons, numberNegativeIons, numWaters)

 
	##################################################
        # Step 6. Minimize and Equilibrate

	# run minimimzation once more with the new ions
	self.mdrun('minimize.mdp')
	
	# Make a new index file for the simulation
	self.makeIndexFile(self.files.grofile, self.files.ndxfile)

	# run equilibration phase 1: just water with frozen protein.  
	self.mdrun('equilibrateWithFrozenProtein.mdp', useIndexFile=True)
			
	# run equilibration phase 2: the whole system can move  
	self.mdrun('equilibrate.mdp', useIndexFile=True)

	
	##################################################
        # Step 7.   Build a tpr for a production run
	self.buildTpr('production.mdp', useIndexFile=True)

	
	##################################################
        # Step 8. Finalize the preparation: 1) copy files to output directory, 2) write logfile, and 3) clean up temp files
	self.finalizePreparation(cleanTempFiles=cleanup)   # this is defined in the parent class
	

	os.chdir( OriginalLocationDir )	### change back out of the working directory to our original location
	return
	
	
    def useTIP3P(self):
        """Change the appropriate lines in the topology file to use ffamber_tip3p."""
    
        print '>> Using TIP3P water <<'
	top = TopologyFile(topfile=self.files.topfile)
	
	# print '>>> top.includes:', top.includes
	# print '>>> top.ifdefs:', top.ifdefs

	# Remove other water models
	i = 0
	TIP3PLine = False
	while i < len(top.includes):
	    if top.includes[i].count('"ffamber_tip3p.itp"') > 0:
		TIP3PLine = True
	    if top.includes[i].count('"spc.itp"') > 0:
		print '>>> REMOVING include:', top.includes.pop(i) # remove any instances of '#include "spc.itp"'
	    elif top.includes[i].count('"flexspc.itp"') > 0:
		print '>>> REMOVING include:', top.includes.pop(i) # remove any instances of '#include "flexspc.itp"'
	    else:
		i += 1
		
        # Remove any "spc.itp" mention from ifdefs too
	#
	# Recall that top.ifdefs are lists of ifdef groups of lines
	# Example:
	#     [ ['#ifdef POSRES\n', '#include "posre.itp"\n', '#endif\n'],
	#       ['#ifdef FLEX_SPC\n', '#include "flexspc.itp"\n', '#else\n', '#include "spc.itp"\n', '#endif\n'], .... ]
	i = 0
	while i < len(top.ifdefs):
	    spcMentioned = False
	    for line in top.ifdefs[i]:
		if line.count('spc.itp') > 0:
		    spcMentioned = True
		    break
	    if spcMentioned:
		top.ifdefs.pop(i)
	    else:
		i += 1
	
	# If there's no line for ffamber_tip3p.itp, add it
	if not TIP3PLine:
	    top.includes.append('#include "ffamber_tip3p.itp"\n')
	    
	# write the new topfile
	top.rebuild_lines()  # This must be done if there are any changes made the the data
	top.write(self.files.next_top())
	self.files.increment_top()
	
      

    def calcTotalCharge(self):
	"""Use a clever application of editconf to calculate the total charge in of the system."""
	
	# does a *.tpr file exist?
	tprFileExists = os.path.exists(self.files.tprfile)
	
	# build a tpr if one doesn't exist
	if not tprFileExists:
    	    mdpfile = self.findMdpfile('minimize.mdp')
	    grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], mdpfile, self.files.grofile, self.files.tprfile, self.files.topfile)
	    self.rungmx( grompp )		    
	    
	# create a *.pdb file where the B-factor column contains the charge of each atom
	pdbfile  = tempfile.mktemp()+'.pdb'
	editconf = '%s/editconf -f %s -o %s -mead'%(os.environ['GMXPATH'], self.files.tprfile, pdbfile)
	self.rungmx( editconf )
	
        # Load the pdbfile and add up the atomic charge values in the occupancy column 
	# Note -- because of the non-standard PDB file output by editconf, we are
	# parsing this by hand rather than use Structure.PDBStructure
	#
	# Parse pdbfile.
	# Example:                                                *chgs*
	# ATOM  13472  HW1 SOL  4153      49.750  47.060   0.920  0.4170  0.4000
	# ATOM  13473  HW2 SOL  4153      51.110  46.410   0.970  0.4170  0.4000
	# ATOM  13474  OW  SOL  4154      56.190  48.230  28.640 -0.8340  1.0500
	# ...

	fin = open(pdbfile,'r')
	pdblines = fin.readlines()
	fin.close()
	
	totalCharge = 0.0
	for line in pdblines:
	    lastTwoColumnFields = line[54:].split()
	    if len(lastTwoColumnFields) == 2:
	        totalCharge += float(lastTwoColumnFields[0])
	
	# remove temporary files we don't need anymore
	os.remove(pdbfile)
 
	# return the total charge
	return totalCharge


    def AddIons(self, np, nn, nwaters):
	"""Use the gromacs program genbox to replace water molecules in the current simulation box
	with ions.
	
	REQUIRED
	np        - the number of positive ions species 
	nn        - the number of negative ions species 
	nwaters   - the number of water molecules *currently* in the box  
	
	After the ions are added (and the water molecules replaced), there should be:
	np cations, nn anions, and (nwaters-nn-np) water molecules. 
	"""
    
        # If we're using the 3.1.4 version of gromacs, we first need to find the index of the solvent group
	# that genbox will use to add the ions.
	if self.system.gmx_version == '3.1.4':
	    
	    # Get the solvent group index
	    solventGroupIndex = self.getSolventGroupIndexFromGenion()
	    	
	    # Make a new grofile with the added ions	   
	    genion = 'echo %d | %s/genion -s %s -o %s -pname %s -np %d -pq %d -nname %s -nn %d -nq %d -g genion.log -n %s '%(solventGroupIndex, os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro(), self.system.positiveIonName, np, self.system.positiveIonCharge, self.system.negativeIonName, nn, self.system.negativeIonCharge, self.files.ndxfile)
	    self.rungmx( genion )
	    self.files.increment_gro()    # must increment filename for any new gmx file 
	    
	    
	elif self.system.gmx_version == '3.3':
	    
	    # make a new grofile with the added ions
	    genion = '%s/genion -s %s -o %s -pname %s -np %d -pq %d -nname %s -nn %d -nq %d -g genion.log'%(os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro(), self.system.positiveIonName, np, self.system.positiveIonCharge, self.system.negativeIonName, nn, self.system.negativeIonCharge)
	    self.rungmx( genion )
	    if (np+nn) > 0:  # genion 3.3 WILL NOT produce a new gro if there is nothing to be done
		self.files.increment_gro()    # must increment filename for any new gmx file 
	    
	# Generate a new index file for the ion-ated grofile  
	self.makeIndexFile(self.files.grofile, self.files.next_ndx())
	self.files.increment_ndx()    # must increment filename for any new gmx file 

	# Generate a new topology file that we will add ions to (and decrease the number of water molecules)
	top = TopologyFile(topfile=self.files.topfile)

	# Add ions to topology file, which will eventually be in blocks that look like so:
	# [ moleculetype ]
	# ; molname       nrexcl
	# K               1
	#
	# [ atoms ]
	# ; id    at-type res-nr  residuename     at-name  cgnr   charge    mass
	# 1       K       1       K               K        1       1        39.098

	ions = IonEnvironmentCalculation(self.system).getIons()  # a list of Ion() objects we will use to build ion molecules
	
	# Formatting needed for TopologyFileMolecule lines
	blockText = """[ moleculetype ]
; molname       nrexcl
%(molname)s               1

[ atoms ]
; id    atomtype        resnum  residuename    atoname       cgnr   charge      mass
1       %(atom_type)s  1       %(molname)s    %(molname)s   1       %(charge)d  %(mass)f """
	
	# Make a TopologyMolecule object for the negative ion
	negIonDict = {}
	negIonDict['molname']   = self.system.negativeIonName
	negIonDict['atom_type'] = ions[self.system.negativeIonName].atomType
	negIonDict['charge']    = ions[self.system.negativeIonName].charge
	negIonDict['mass']      = ions[self.system.negativeIonName].atomicMass
	print 'negIonDict', negIonDict
	lines = [s+'\n' for s in ((blockText%negIonDict).split('\n')) ]
	print 'lines', lines

	lines = (blockText%negIonDict).split('\n')
	#lines = [(s+'\n') for s in ((blockText%negIonDict).split('\n')) ]
	print 'lines', lines
	negIonMolecule = TopologyFileMolecule( lines )
	
	# Make a TopologyMolecule object for the positive ion
	posIonDict = {}
	posIonDict['molname']   = self.system.positiveIonName
	posIonDict['atom_type'] = ions[self.system.positiveIonName].atomType
	posIonDict['charge']    = ions[self.system.positiveIonName].charge
	posIonDict['mass']      = ions[self.system.positiveIonName].atomicMass
	lines = (blockText%posIonDict).split('\n')
	# lines = [(s+'\n') for s in () ]
	print 'lines', lines
	posIonMolecule = TopologyFileMolecule( lines )

	# Add the ion molecules to the topology file
	top.AddMolecule(negIonMolecule, nn)   # arguments: opologyFileMolecule object, NumberOfMolecules 
	top.AddMolecule(posIonMolecule, np)

	# The number of water molecules needs to be changed to reflect the current count
	numberWatersAfterGenion = nwaters-nn-np
	top.nmols['SOL'] = str(numberWatersAfterGenion)
	
	# Rebuild all the lines in the topfile and save it
	top.rebuild_lines()
	top.write(self.files.next_top())
	self.files.increment_top()    # PrepWorker must increment filename for any new gmx file 

   
    
    
	
    def makeIndexFile(self, grofile, ndxfile, protein=True, ions=True, sol=True, bath=True, system=True, debug=False):
    
	"""Reads in the grofile, and make an *.ndx file with atom groups corresponding to a list of commmon-sense names
	given in the groups=[] list.  
	
	REQUIRED
	grofile   - gromacs structure filename from which to compile the atom groups
	ndxfile   - filename to write the index file to 
	
	Supported group names so far:
	
	    protein            All protein atoms
	    sol                All solvent atoms.  (protein + sol + ions) shouuld equal (system) !
	    ions               All ions
	    bath               (sol + ions)
	    system             All atoms
	    
	Note: This subroutine assumes it's being called from the current working directory
	
	OUTPUT
	Writes an index file to the given filename. 
	"""
     
        ndx = IndexFile()
	
	all_indices = ndx.getAtomListFromGrofile(grofile)
	ion_indices = ndx.getAtomListFromGrofile(grofile, ResName=['Na','Cl','Na+','Cl-','Mg','Ca'])
	sol_indices = ndx.getAtomListFromGrofile(grofile, ResName=['SOL','HOH'])
	
	# The protein indices are considered to be all the residues besides ions and solvent
	protein_indices = [i for i in all_indices if (ion_indices.count(i) == 0 and sol_indices.count(i) == 0) ]
	
	# The bath is just the ions plus the solvent
	bath_indices = ion_indices + sol_indices
	
	# Create the index file groups
	ndx.addIndexGroup('protein', protein_indices)
	ndx.addIndexGroup('sol', sol_indices)
	ndx.addIndexGroup('ions', ion_indices)
	ndx.addIndexGroup('bath', bath_indices)
	ndx.addIndexGroup('system', all_indices)

	# write the ndxfile 
	ndx.writeIndexFile(ndxfile)	
	
    def getSolventGroupIndexFromGenion(self):
	"""Use the gromacs program "genion" to report the solvent group index."""
	
	# Do we have an indexfile made yet?
	indexFileExists = os.path.exists(self.files.ndxfile)
	if not indexFileExists:
	    # If not, make a standard indexfile for our project, with groups 'protein', 'ions', 'sol', 'bath', 'system'
	    self.makeIndexFile(self.files.grofile, self.files.ndxfile)
	    
	# We then run genion, with 'q' (quit) as an input argument, just so we can parse the output
	genion = 'echo q | %s/genion -s %s -o %s -pname %s -np %d -pq %d -nname %s -nn %d -nq %d -g genion.log -n %s'%(os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro(), self.system.positiveIonName, 1, self.system.positiveIonCharge, self.system.negativeIonName, 1, self.system.negativeIonCharge, self.files.ndxfile)
	genion_lines = commands.getoutput(genion).split('\n')
	
	# The output will look something like this:
	# ....
	# Select a continuous group of solvent molecules
	# Group     0 (     protein) has  1254 elements
	# Group     1 (      system) has  1254 elements
	# (etc...)

	# Go through each of these output lines to find the 'sol' Group idex 
	solventGroupIndex = None
	for line in genion_lines:
	    if line[0:5] == 'Group':
		fields = line.replace(')',' ').replace('(',' ').split()
		if fields[2] == 'sol':
		    solventGroupIndex = int(fields[1])
		    
	# If we can't find a solvent Group index, something has gon horribly awry :(		
	if solventGroupIndex == None:
	    print "Cannot find solvent Group index number for 'sol' from genion output:"
	    for line in genion_lines:
	       print '>>>>', line
	    raise GenionParsingError
	
	return solventGroupIndex	
    
	
# Exceptions
class GenionParsingError(Exception):
    print "There was an error parsing the information from the gmx program 'genion'."
        
