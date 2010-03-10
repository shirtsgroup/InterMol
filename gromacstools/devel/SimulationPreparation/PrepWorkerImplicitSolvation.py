
# ----------------------------------------------------------------------
# IMPORTS
# ----------------------------------------------------------------------

import sys, os, string, commands, glob, random
import tempfile

from PrepWorker import *
from SystemSetup import *

from mmtools.gromacstools.devel.Helpers.IndexFile import *


# ----------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# CLASSES
# ----------------------------------------------------------------------

class PrepWorkerImplicitSolvation(PrepWorker):
    """An object to store info and coordinate a preparation run.
    
    This class inherents methods from PrepWorker:
        def __init__()
        def mdrun(self, mdpfile)
        def cleanupWorkingDir()
        def rungmx()
        def copyFilesToOutputDir()
    """
  

    def prepare(self, cleanup=True, verbose=False, debug=False, checkForFatalErrors=True, mockrun=False):
    
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
    
	checkForFatalErrors True        Exits the program if any of the gmx produce "Fatal error" in the standard output	    
	mockrun             False       if True, prints out the gromacs commands, but does does execture them  
				    
	 
	OUTPUT
      
	The following files will be written to the outdir:
	  
	    out.tpr		GROMACS tpr file
	    out.gro		GROMACS gro (structure) file
	    equilibrate*	executable script to call mdrun 
	    history.log         logfile containing all the commands used to prepare the simualtion 
	    output.log          logfile containg the output of the pre-processing programs
	"""
	
	self.checkForFatalErrors = checkForFatalErrors
	self.mockrun = mockrun
	
	# remember our original location before we descend into the working directory
        OriginalLocationDir = os.path.abspath( os.curdir )
	
	# From now on, do work in the working directory
	os.chdir( self.workdir )
	print 'Preparing simulation in directory %s...'%self.workdir


	# Step 0.  Initialize the parameters for this simulation system
	system = SystemSetup()

	# step 1.  Create a grofile from the input file (which is either a *.pdb or *.gro)
	if system.version == '3.1.4':	
	    pdb2gmx = 'echo %s | %s/pdb2gmx -f %s -o %s -p %s -ignh -H14'%(system.forcefieldCodes[system.forcefield], os.environ['GMXPATH'], self.files.infile, self.files.grofile, self.files.topfile)
	elif self.version == '3.3':
	    pdb2gmx = '%s/pdb2gmx -ff %s -f %s -o %s -p %s -ignh'%(os.environ['GMXPATH'], system.forcefield, self.files.infile, self.files.grofile, self.files.topfile)
	self.rungmx( pdb2gmx )
	
	# Step 2.  Since this is implicit solvent, remove any references to SPC or TIP3P water in the topology file
        ### TopologyFile.useTIP3P( self.files.topfile )
	    
        # Step 3. make a (periodic boundary conditions) box of an appropriate size
        if system.useAbsBoxSize:
	    editconf = '%s/editconf -bt %s -f %s -o %s -box %s'%(os.environ['GMXPATH'], system.boxType, self.files.grofile, self.files.next_gro(), self.setup.absBoxSize )
        else:
       	    editconf = '%s/editconf -bt %s -f %s -o %s -d %s'%(os.environ['GMXPATH'], system.boxType, self.files.grofile, self.files.next_gro(), system.boxSoluteDistance )
	self.rungmx( editconf )
	self.files.increment_gro()    # PrepWorker must increment filename for any new gmx file 

        # step 4.  self.implicitSolvationPreparationSteps(implicitOptions=implicitOptions, useTable=useTable)

	# 4a. Load in an *.mdp file for the implicit prep
	# 4b. Find the shortest box vector
	# 4c. Use the is length to determine 

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
	
	
	
	    
    def make_ndx(self, grofile, ndxfile, protein=True, ions=True, sol=True, bath=True, system=True, debug=False):
    
	"""Make an *.ndx file with atom groups corresponding to a list of commmon-sense names
	given in the groups=[] list.  Supported group names so far:
	
	    protein            All protein atoms
	    sol                All solvent atoms.  (protein + sol + ions) shouuld equal (system) !
	    ions               All ions
	    bath               (sol + ions)
	    system             All atoms
	    
	Note: This subroutine assumes it's being called from the current working directory
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
    

