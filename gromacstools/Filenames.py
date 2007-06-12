# Modification history:
#
# VAV:  June 11, 2007:  Added set_mdpMinimization(), set_mdpEquilibration(), and set_mdpSimulation()


import sys, os, tempfile, string



class Filenames(object):
    
  def __init__(self, infile, workdir): 
    """A class to keep track of the most current gmx input and output files
    as they get puched through the pipeline.
    
    INPUT
    
    infile        the input *.pdb or *.gro file from which to build the simulation tpr
    workdir       the working (temporary) directory in which to do the processing    
    """

    # process infile name
    self.infile = infile
    self.infile_basename = os.path.basename(self.infile) 
    self.infile_suffix = self.infile.split('.').pop()
 
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
    
    # initialize default mdpfiles
    self.mdpfile_Minimization  = os.path.join(self.MMTOOLSPATH,'gromacstools/mdp/minimize.mdp')
    self.mdpfile_Equilibration = os.path.join(self.MMTOOLSPATH,'gromacstools/mdp/equilibrate.mdp')
    self.mdpfile_Simulation    = os.path.join(self.MMTOOLSPATH,'gromacstools/mdp/simulate.mdp')
    
    # process workdir
    self.workdir = workdir
      
    # FILES #
    
    self.grofile = ''
    self.pdbfile = ''
    self.topfile = ''
    self.tprfile = ''
    self.ndxfile = ''
    
    # initialize file names
    # NOTE: These names do not correspond to files yet, but will be used
    #       by self.prepare() as a consistent naming scheme when building the tpr 
    
    if self.infile[-4:] == '.pdb':     
        self.grofile = os.path.join(self.workdir, self.infile_basename.replace('.pdb','0000.gro') )
        self.pdbfile = os.path.join(self.workdir, self.infile_basename.replace('.pdb','0000.pdb') )
        self.topfile = os.path.join(self.workdir, self.infile_basename.replace('.pdb','0000.top') )
        self.tprfile = os.path.join(self.workdir, self.infile_basename.replace('.pdb','0000.tpr') )
        self.ndxfile = os.path.join(self.workdir, self.infile_basename.replace('.pdb','0000.ndx') )
    elif self.infile[-4:] == '.gro':
        self.grofile = os.path.join(self.workdir, self.infile_basename.replace('.gro','0000.gro') )
        self.pdbfile = os.path.join(self.workdir, self.infile_basename.replace('.gro','0000.pdb') )
        self.topfile = os.path.join(self.workdir, self.infile_basename.replace('.gro','0000.top') )
        self.tprfile = os.path.join(self.workdir, self.infile_basename.replace('.gro','0000.tpr') )
        self.ndxfile = os.path.join(self.workdir, self.infile_basename.replace('.gro','0000.ndx') )
    else:
        print 'Error:  infile must be either *.pdb or *.gro file'
	sys.exit(1)
	
    oldway = """
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
    
    """

  def increment(self, fileExtensions):
    """For each file extension in the list fileExtensions=[ 'gro', 'tpr', ...],
    increment the the fiilename from "filename0001.gro" --> "filename0002.gro", e.g."""
    
    for ext in fileExtensions:
      if ext=='gro':
	self.increment_gro()
      elif ext=='pdb':
	self.increment_pdb()
      elif ext=='top':
	self.increment_top()
      elif ext=='tpr':
	self.increment_tpr()
      elif ext=='ndx':
	self.increment_ndx()
      else:
	print 'File extension',ext,'not supported!  Exiting....'
	sys.exit(1)
	
  def increment_gro(self):
    """Increments the current gro file"""
    self.grofile = self.grofile[0:-8] + string.zfill( (int(self.grofile[-8:-4])+1), 4) + self.grofile[-4:]
        
  def increment_pdb(self):
    """Increments the current pdb file"""
    self.pdbfile = self.pdbfile[0:-8] + string.zfill( (int(self.pdbfile[-8:-4])+1), 4) + self.pdbfile[-4:]
    
  def increment_top(self):
    """Increments the current top file"""
    self.topfile = self.topfile[0:-8] + string.zfill( (int(self.topfile[-8:-4])+1), 4) + self.topfile[-4:]
    
  def increment_tpr(self):
    """Increments the current tpr file"""
    self.tprfile = self.tprfile[0:-8] + string.zfill( (int(self.tprfile[-8:-4])+1), 4) + self.tprfile[-4:]
    
  def increment_ndx(self):
    """Increments the current ndx file"""
    self.ndxfile = self.ndxfile[0:-8] + string.zfill( (int(self.ndxfile[-8:-4])+1), 4) + self.ndxfile[-4:]
    
  def next_gro(self):
    """Returns the next name of the grofile, without incrementing"""
    return self.grofile[0:-8] + string.zfill( (int(self.grofile[-8:-4])+1), 4) + self.grofile[-4:]
  
  def next_pdb(self):
    """Returns the next name of the pdbfile, without incrementing"""
    return self.pdbfile[0:-8] + string.zfill( (int(self.pdbfile[-8:-4])+1), 4) + self.pdbfile[-4:]
  
  def next_top(self):
    """Returns the next name of the topfile, without incrementing"""
    return self.topfile[0:-8] + string.zfill( (int(self.topfile[-8:-4])+1), 4) + self.topfile[-4:]
  
  def next_tpr(self):
    """Returns the next name of the tprfile, without incrementing"""
    return self.tprfile[0:-8] + string.zfill( (int(self.tprfile[-8:-4])+1), 4) + self.tprfile[-4:]
  
  def next_ndx(self):
    """Returns the next name of the ndxfile, without incrementing"""
    return self.ndxfile[0:-8] + string.zfill( (int(self.ndxfile[-8:-4])+1), 4) + self.ndxfile[-4:]

  
  def set_mdpMinimization(self, filename):
    """Sets the *.mdp filename in the gromacstools/mdp directory to use for minimzation."""
    self.mdpfile_Minimization  = os.path.join(self.MMTOOLSPATH,'gromacstools/mdp/%s'%filename)

  def set_mdpEquilibration(self, filename):
    """Sets the *.mdp filename in the gromacstools/mdp directory to use for equilibration."""
    self.mdpfile_Equilibration = os.path.join(self.MMTOOLSPATH,'gromacstools/mdp/%s'%filename)

  def set_mdpSimulation(self, filename):
    """Sets the *.mdp filename in the gromacstools/mdp directory to use for simulation."""
    self.mdpfile_Simulation    = os.path.join(self.MMTOOLSPATH,'gromacstools/mdp/%s'%filename)    

       
  def show(self):
    """A printout of the current GROMACS filenames"""
    
    print
    print 'Files currently being used in the GROMACS pre-processing:'
    print '%-20s %-20s'%('grofile',self.grofile)
    print '%-20s %-20s'%('pdbfile',self.pdbfile)
    print '%-20s %-20s'%('topfile',self.topfile)
    print '%-20s %-20s'%('tprfile',self.tprfile)
    print '%-20s %-20s'%('ndxfile',self.ndxfile)
    print '%-20s %-20s'%('mdpfile_Minimization',self.mdpfile_Minimization)
    print '%-20s %-20s'%('mdpfile_Equilibration',self.mdpfile_Equilibration)
    print '%-20s %-20s'%('mdpfile_Simulation',self.mdpfile_Simulation)

