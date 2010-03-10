# MdpFile.py
#
# HISTORY
# JDC: 1/29/08    Added capability to initialize from contents of an .mdp file.  Added random seed handling.
# VAV: 9/09/07    Added self.delParameter(), allowing keywords and their valuestrs to be removed from the mdp file
# VAV: 9/02/07    Added self.use_table() and TableMaker object 
# VAV: 8/27/07    File created
#
# TO DO:
# * Have default *.mdp template in case __init__() template file is not provided.

import sys, os, tempfile, string, copy
from TableMaker import *


class MdpFile(object):
    """A class to store and manipulate Gromacs mdp files.
  
    TODO
      Can we make use of the Units class to allow user to specify parameters with units attached and automatically convert to/from gromacs units?
    """
  
    # Class data
    MAXSEED = 2**30 # maximum random number seed
    
      
    def __init__(self, mdpfile = None, debug = False): 
        """Constructor to initialize mdp file container.
      
        OPTIONAL ARGUMENTS
          mdpfile        filename (string) or contents (list of strings) of a gromacs .mdp file to be read
        """
    
        # set debug flag
        self.debug = debug 
        self.mdpfile = mdpfile
                  
        # initialize internal storage
        self.params = {}    # {keyword:value} dictionary.
    
        # read .mdp file or process contents, if give
        if mdpfile:
            self.read(mdpfile)
        else:
            self.setDefaults()
        
        return
    
    def setDefaults(self):
        """Set default parameters for the mdp file."""
        
        self.params['title']         = 'Title'
        self.params['cpp']           = '/lib/cpp -traditional'
        self.params['include']       = ''
        self.params['define']        = ''
        
        self.params['integrator']    = 'md'
        self.params['tinit']         = '0'
        self.params['dt']            = '0.002'
        self.params['nsteps']        = '1000'
        self.params['comm_mode']     = 'linear'
        self.params['nstcomm']       = '1'
        self.params['comm_grps']     = 'system'
        
        self.params['bd-temp']       = '300'
        self.params['bd-fric']       = '0'
        self.params['ld-seed']       = '1993'

        self.params['emtol']         = '100'
        self.params['emstep']        = '0.01'
        self.params['niter']         = '20'
        self.params['fcstep']        = '0'
        self.params['nstcgsteep']    = '1000'
        
        self.params['nstxout']       = '100'
        self.params['nstvout']       = '0'
        self.params['nstfout']       = '0'
        self.params['nstlog']        = '100'
        self.params['nstenergy']     = '100'
        self.params['nstxtcout']     = '100'
        self.params['xtc-precision'] = '1000'
        self.params['xtc_grps']      = 'system'
        self.params['energygrps']    = 'system'
        
        
        # NEIGHBORSEARCHING PARAMETERS
        self.params['nstlist']       = '10'
        self.params['ns_type']       = 'grid'
        self.params['pbc']           = 'xyz'
        self.params['rlist']         = '0.9'
        self.params['domain-decomposition']         = 'no'
        
        # OPTIONS FOR ELECTROSTATICS AND VDW 
        self.params['coulombtype']   = 'PME'
        self.params['rcoulomb-switch']= '0.8'
        self.params['rcoulomb']      = '0.9'
        self.params['epsilon_r']     = '78'
        self.params['vdwtype']       = 'switch'
        self.params['rvdw_switch']   = '0.8'
        self.params['rvdw']          = '0.9'
        self.params['DispCorr']      = 'Ener'
        self.params['fourierspacing']= '0.08'
        self.params['fourier_nx']    = '0'
        self.params['fourier_ny']    = '0'
        self.params['fourier_nz']    = '0'
        self.params['pme_order']     = '6'
        self.params['ewald_rtol']    = '1e-06'
        self.params['ewald_geometry']= '3d'
        self.params['epsilon_surface'] = '0'
        self.params['optimize_fft']  = 'no'
        
        # OPTIONS FOR WEAK COUPLING ALGORITHMS = 
        self.params['tcoupl']        = 'nose-hoover'        
        self.params['tc_grps']       = 'system'        
        self.params['tau_t']         = '0.0109'
        self.params['ref_t']         = '300.0'
        self.params['Pcoupl']        = 'No'        
        self.params['Pcoupltype']    = 'Isotropic'        
        self.params['tau-p']         = '1.0'
        self.params['compressibility']   = '0'
        self.params['ref-p']         = '1.0'
        
        # SIMULATED ANNEALING CONTROL = 
        self.params['annealing']     = 'no'
        self.params['zero-temp_time']= '0'
  
        # GENERATE VELOCITIES FOR STARTUP RUN = 
        self.params['gen_vel']       = 'yes'
        self.params['gen_temp']      = '300'
        self.params['gen_seed']      = '-1'

        # OPTIONS FOR BONDS     = 
        self.params['constraints']   = 'hbonds'
        self.params['constraint-algorithm']   = 'lincs'
        self.params['unconstrained-start']   = 'no'
        self.params['Shake-SOR']    = 'no'
        self.params['shake-tol']     = '1e-04'
        self.params['lincs-order']   = '4'
        self.params['lincs-warnangle']= '30'
        self.params['morse']         = 'no'
  
        # ENERGY GROUP EXCLUSIONS 
        self.params['energygrp_excl']=''
          
        # NMR refinement stuff 
        self.params['disre']          = 'No'
        self.params['disre-weighting']= 'Conservative'
        self.params['disre-mixed']    = 'no'
        self.params['disre-fc']       = '1000'
        self.params['disre-tau']      = '0'
        self.params['nstdisreout']    = '100'
        self.params['orire']          = 'no'
        self.params['orire-fc']       = ''
        self.params['orire-tau']      = ''
        self.params['orire-fitgrp']   = ''
        self.params['nstorireout']    = '100'
  
        # Free energy control stuff = 
        self.params['free-energy']    = 'no'
        self.params['init-lambda']    = '0'
        self.params['delta-lambda']   = '0'
        self.params['sc-alpha']       = '0'
        self.params['sc-sigma']       = '0.3'
  
        # Non-equilibrium MD stuff = 
        self.params['acc-grps']       = ''
        self.params['accelerate']     = ''
        self.params['freezegrps']     = '' 
        self.params['freezedim']      = '' 
  
        self.params['cos-acceleration']         = '0'
  
        # Electric fields       = 
        # ; Format is number of terms (int) and for all terms an amplitude (real) = 
        # ; and a phase angle (real) = 
        self.params['E-x']            = ''
        self.params['E-xt']           = ''
        self.params['E-y']            = ''
        self.params['E-yt']           = ''
        self.params['E-z']            = ''
        self.params['E-zt']           = ''
  
        # User defined thingies = 
        self.params['user1-grps']     = ''
        self.params['user2-grps']     = ''
        self.params['userint1']       = ''
        self.params['userint2']       = ''
        self.params['userint3']       = ''
        self.params['userint4']       = ''
        self.params['userreal1']      = ''
        self.params['userreal2']      = ''
        self.params['userreal3']      = ''
        self.params['userreal4']      = ''      
        
        # Finally, compile a list of the default keys, we can later cross-check to see if we missed any non-default params
        self.defaultKeys = copy.copy(self.params.keys())

        
    def template(self):
        """Returns a template string for formatting the output"""
        
        return """;
  ;	mdp parameters 
  ;
  
  ; VARIOUS PREPROCESSING OPTIONS = 
  title                    = %(title)s
  cpp                      = %(cpp)s
  include                  = %(include)s
  define                   = %(define)s
  
  ; RUN CONTROL PARAMETERS = 
  integrator               = %(integrator)s
  ; start time and timestep in ps = 
  tinit                    = %(tinit)s
  dt                       = %(dt)s
  nsteps                   = %(nsteps)s
  ; mode for center of mass motion removal = 
  comm_mode                = %(comm_mode)s
  ; number of steps for center of mass motion removal = 
  nstcomm                  = %(nstcomm)s
  ; group(s) for center of mass motion removal = 
  comm_grps                = %(comm_grps)s
  
  ; LANGEVIN DYNAMICS OPTIONS = 
  ; Temperature, friction coefficient (amu/ps) and random seed = 
  bd-temp                  = %(bd-temp)s
  bd-fric                  = %(bd-fric)s
  ld-seed                  = %(ld-seed)s
  
  ; ENERGY MINIMIZATION OPTIONS = 
  ; Force tolerance and initial step-size = 
  emtol                    = %(emtol)s
  emstep                   = %(emstep)s
  ; Max number of iterations in relax_shells = 
  niter                    = %(niter)s
  ; Step size (1/ps^2) for minimization of flexible constraints = 
  fcstep                   = %(fcstep)s
  ; Frequency of steepest descents steps when doing CG = 
  nstcgsteep               = %(nstcgsteep)s
  
  ; OUTPUT CONTROL OPTIONS = 
  ; Output frequency for coords (x), velocities (v) and forces (f) = 
  nstxout                  = %(nstxout)s
  nstvout                  = %(nstvout)s
  nstfout                  = %(nstfout)s
  ; Output frequency for energies to log file and energy file = 
  nstlog                   = %(nstlog)s
  nstenergy                = %(nstenergy)s
  ; Output frequency and precision for xtc file = 
  nstxtcout                = %(nstxtcout)s
  xtc-precision            = %(xtc-precision)s
  ; This selects the subset of atoms for the xtc file. You can = 
  ; select multiple groups. By default all atoms will be written. = 
  xtc_grps                 = %(xtc_grps)s
  ; Selection of energy groups = 
  energygrps               = %(energygrps)s
  
  ; NEIGHBORSEARCHING PARAMETERS = 
  ; nblist update frequency = 
  nstlist                  = %(nstlist)s
  ; ns algorithm (simple or grid) = 
  ns_type                  = %(ns_type)s
  ; Periodic boundary conditions: xyz or no = 
  pbc                      = %(pbc)s
  ; nblist cut-off         = 
  rlist                    = %(rlist)s
  domain-decomposition     = %(domain-decomposition)s
  
  ; OPTIONS FOR ELECTROSTATICS AND VDW = 
  ; Method for doing electrostatics = 
  coulombtype              = %(coulombtype)s
  rcoulomb-switch          = %(rcoulomb-switch)s
  rcoulomb                 = %(rcoulomb)s
  
  ; Dielectric constant (DC) for cut-off or DC of reaction field = 
  epsilon_r                = %(epsilon_r)s
  ; Method for doing Van der Waals = 
  vdwtype                  = %(vdwtype)s
  ; cut-off lengths        = 
  rvdw_switch              = %(rvdw_switch)s
  rvdw                     = %(rvdw)s
  ; Apply long range dispersion corrections for Energy and Pressure = 
  DispCorr                 = %(DispCorr)s
  ; Spacing for the PME/PPPM FFT grid = 
  fourierspacing           = %(fourierspacing)s
  ; FFT grid size, when a value is 0 fourierspacing will be used = 
  fourier_nx               = %(fourier_nx)s
  fourier_ny               = %(fourier_ny)s
  fourier_nz               = %(fourier_nz)s
  ; EWALD/PME/PPPM parameters = 
  pme_order                = %(pme_order)s
  ewald_rtol               = %(ewald_rtol)s
  ewald_geometry           = %(ewald_geometry)s
  epsilon_surface          = %(epsilon_surface)s
  optimize_fft             = %(optimize_fft)s
  
  ; OPTIONS FOR WEAK COUPLING ALGORITHMS = 
  ; Temperature coupling   = 
  tcoupl                   = %(tcoupl)s
  ; Groups to couple separately = 
  tc_grps                  = %(tc_grps)s
  ; Time constant (ps) and reference temperature (K) = 
  tau_t                    = %(tau_t)s
  ref_t                    = %(ref_t)s
  ; Pressure coupling      = 
  Pcoupl                   = %(Pcoupl)s
  Pcoupltype               = %(Pcoupltype)s
  ; Time constant (ps), compressibility (1/bar) and reference P (bar) = 
  tau-p                    = %(tau-p)s
  compressibility          = %(compressibility)s
  ref-p                    = %(ref-p)s
  
  ; SIMULATED ANNEALING CONTROL = 
  annealing                = %(annealing)s
  ; Time at which temperature should be zero (ps) = 
  zero-temp_time           = %(zero-temp_time)s
  
  ; GENERATE VELOCITIES FOR STARTUP RUN = 
  gen_vel                  = %(gen_vel)s
  gen_temp                 = %(gen_temp)s
  gen_seed                 = %(gen_seed)s
  
  ; OPTIONS FOR BONDS     = 
  constraints              = %(constraints)s
  ; Type of constraint algorithm = 
  constraint-algorithm     = %(constraint-algorithm)s
  ; Do not constrain the start configuration = 
  unconstrained-start      = %(unconstrained-start)s
  ; Use successive overrelaxation to reduce the number of shake iterations = 
  Shake-SOR                = %(Shake-SOR)s
  ; Relative tolerance of shake = 
  shake-tol                = %(shake-tol)s
  ; Highest order in the expansion of the constraint coupling matrix = 
  lincs-order              = %(lincs-order)s
  ; Lincs will write a warning to the stderr if in one step a bond = 
  ; rotates over more degrees than = 
  lincs-warnangle          = %(lincs-warnangle)s
  ; Convert harmonic bonds to morse potentials = 
  morse                    = %(morse)s  
  ; ENERGY GROUP EXCLUSIONS = 
  ; Pairs of energy groups for which all non-bonded interactions are excluded = 
  energygrp_excl           = %(energygrp_excl)s 
  
  ; NMR refinement stuff  = 
  ; Distance restraints type: No, Simple or Ensemble = 
  disre                    = %(disre)s
  ; Force weighting of pairs in one distance restraint: Conservative or Equal = 
  disre-weighting          = %(disre-weighting)s
  ; Use sqrt of the time averaged times the instantaneous violation = 
  disre-mixed              = %(disre-mixed)s
  disre-fc                 = %(disre-fc)s
  disre-tau                = %(disre-tau)s
  ; Output frequency for pair distances to energy file = 
  nstdisreout              = %(nstdisreout)s
  ; Orientation restraints: No or Yes = 
  orire                    = %(orire)s
  ; Orientation restraints force constant and tau for time averaging = 
  orire-fc                 = %(orire-fc)s
  orire-tau                = %(orire-tau)s
  orire-fitgrp             = %(orire-fitgrp)s
  ; Output frequency for trace(SD) to energy file = 
  nstorireout              = %(nstorireout)s
  
  ; Free energy control stuff = 
  free-energy              = %(free-energy)s
  init-lambda              = %(init-lambda)s
  delta-lambda             = %(delta-lambda)s
  sc-alpha                 = %(sc-alpha)s
  sc-sigma                 = %(sc-sigma)s
  
  ; Non-equilibrium MD stuff = 
  acc-grps                 = %(acc-grps)s
  accelerate               = %(accelerate)s
  freezegrps               = %(freezegrps)s
  freezedim                = %(freezedim)s
  
  cos-acceleration         = %(cos-acceleration)s
  
  ; Electric fields       = 
  ; Format is number of terms (int) and for all terms an amplitude (real) = 
  ; and a phase angle (real) = 
  E-x                      = %(E-x)s
  E-xt                     = %(E-xt)s
  E-y                      = %(E-y)s
  E-yt                     = %(E-yt)s
  E-z                      = %(E-z)s
  E-zt                     = %(E-zt)s
  
  ; User defined thingies = 
  user1-grps               = %(user1-grps)s
  user2-grps               = %(user2-grps)s
  userint1                 = %(userint1)s
  userint2                 = %(userint2)s
  userint3                 = %(userint3)s
  userint4                 = %(userint4)s
  userreal1                = %(userreal1)s
  userreal2                = %(userreal2)s
  userreal3                = %(userreal3)s
  userreal4                = %(userreal4)s
  """
        
    def write(self, mdpfile):
        """Write an mdpfile to file, using the template.  If there are non-default parameters, these will get tacked on at the end."""
        
        fout = open(mdpfile,'w')
        
        # write params corresponding to the default template
        fout.write( self.template()%self.params )
        
        # check to see if there's any non-standard params
        nonstandard_keys = [key for key in self.params.keys() if self.defaultKeys.count(key) == 0 ]
        
        # write these non-standard parms at the end of the mdpfile
        fout.write('\n  ; Non-standard parameters\n')
        for key in nonstandard_keys:
            fout.write('  %s   = %r\n'%(key,self.params[key]) )

        # Close the file
        fout.close()

  
    def read(self, mdpfile):
        """Read the contents of a gromacs .mdp file.
    
        ARGUMENTS
          mdpfile (string or list of strings) - if a string is provided, the contents of the file are read;
                                                if a list is provided, this is interpreted as the contents of an mdp file
        """
    
        # handle multiple argument types
        if (type(mdpfile) == list):
            lines = mdpfile
        elif (type(mdpfile) == str):
            # if a filename is provided, get its contents
            fin = open(mdpfile, 'r')
            lines = fin.readlines()
            # remove extraneous newlines 
            lines = map(lambda x:x.rstrip("\n"),lines)
            fin.close()
        else:
            raise "Unrecognized argument: " + repr(mdpfile)
          
        # build an initial dictionary of key:values from the lines
        self.params = self.parseParametersFromMdpLines(lines)    
    
    
    def parseParametersFromMdpLines(self, lines):
        """NEED TO FIX so that argumentsare read in with the correct default type!!!! VAV!
        
        """
        
        params = dict()
        for line in lines:
            # strip comments
            index = line.find(';')
            if index > -1:
                line = line[0:index]
            # strip whitespace from both ends of line
            line = line.strip()
            # split off keyword
            fields = line.split()
            # if key = values ... is found, store
            if (len(fields) >= 3) and (fields[1] == '='):
                # extract keywords
                keyword = fields[0]
                # extract value(s)
                index = line.find('=')
                values = line[index+1:len(line)].strip()
                # store
                params[keyword] = values
            # increment line index
            line_index += 1
    
        return [params, keywords, keyword_line_index]
                
    def showParams(self):
        """Display dictionary of recognized key:value pairs from mdp file.
        """
        outstr = ''
        for key in self.keywords:
            outstr = outstr + self.keyword2line(key)
        print outstr
    
        return
    
  
    def keyword2line(self, keyword):
      """Takes in a parameter keyword and returns a formatted mdp file line.
  
      ARGUMENTS
        keyword (string) - the keyword specifying the formatted mdp line to retrieve
  
      RETURNS
        the formatted mdp line corresponding to the keyword, if found
  
      TODO
        What if keyword isn't found?
        We should store the corresponding original line from the file, or store the complete contents of the line (including comments?)      
  
      """ 
      return '%-20s = %-20s\n'%(keyword, self.params[keyword])
      
    def setParameter(self, keyword, valuestr):
      """Sets a particular mdp keyword to a value (string) in the parameter dictionary, updates mdp lines.
  
      ARGUMENTS
        keyword (string) - the keyword to which the parameter is to be assigned
        valuestr (string) - the corresponding parameter value to assign to the keyword
  
      TODO
        Searching in the parameter file can be facilitated if we store a dictionary of which line a keyword is specified on.
      """
  
      # determine whether the keyword already exists
      keyword_exists = self.keyword_line_index.has_key(keyword)
  
      # set keyword in dictionary (creating or overriding)
      self.params[keyword] = valuestr
  
      # create a new line from this keyword-value pair
      line = self.keyword2line(keyword)
      
      if keyword_exists:
        # if keyword exists, replace the line      
        line_index = self.keyword_line_index[keyword]
        self.lines[line_index] = line
      else:
        # if keyword does not exists, add a line
        self.lines.append(line)      
        self.keywords.append(keyword) # add the keyword to the list            
        self.keyword_line_index[keyword] = len(self.lines) - 1 # store line number of the appended line
  
      return
  
    def getParameter(self, keyword):
        """Returns the value associated with a particular parameter; or None if the parameter is not set"""
        if keyword not in self.params:
            return None
        return self.params[keyword]
          
    def randomizeSeed(self):
      """Randomize the random seed, adding one if none exists.
  
      TODO
        Uses range from 1 to 2**30 for now, but should figure out acceptable range for gromacs.
      """
  
      # choose a new random seed
      import random
      self.setParameter('gen_seed', str(random.randint(1, self.MAXSEED)))
  
      return
    

      
      
#####################
  
class ImplicitOptions314(object):
    """The options for implicit-solvent simulations in gromacs3.1.4 (and gromacs3.3) are encoded in the userint and userreal parameters. 
    There is a tricky set of rules for determining what these parameters should be, so this object was created
    as a data structure with methods to calculate the correct settings.
    """
  
    def __init__(self):
      """Initializes the ImplicitOptions() object.  This provides the parent class for a set of 
      objects all with differnt defaults.  The main default here is:
      
      VVSD
      Still GB
      
  
      Here is the README_ED, as of August 2007:
      
  ; Complete description of parameters (from 17 Jan 2007 README_ED from gromacs-3.1.4-ed)
  ; NOTE not all are enabled on PS3.  Those that are enabled on the PS3 are marked with [PS3]
  ;
  ; userint4 - BIT 0 - AGBNP on/off
  ; userint4 - BIT 1 - GB on/off [PS3]
  ; userint4 - BIT 2 - SASD for GB on/off
  ; userint4 - BIT 3 - Andersen on/off
  ; userint4 - BIT 4 - SASD for AGBNP on/off [PS3]
  ; userint4 - BIT 5 - Sheffield on/off
  ; userint4 - BIT 6 - use the velocity Verlet version of the SD algorithm (vvSD) [PS3]
  - userint4 - BIT 7 - scale charges via sdelta  => gradually tuning off electrostatics (AGBNP & EdGB)
  - userint4 - BIT 8 - turn on/off the GCMC part
  ; userint4 - BIT 9 - Still SA part on/off [PS3]
  ; userint4 - BIT 10 - turns on Still GB [PS3]
  ;
  ; userint1 - GB frequency (used for GB and AGBNP) [PS3]
  ; userint2 - BIT 0 - unused
  ; userint2 - BIT 1 - GB logging on/off
  ; userint2 - BIT 2 - AGBNP logging on/off
  ; usetint2 - BIT 3 - do NOT apply the solvation forces, this is helpful for comparisons
  ; userint3 - Andersen frequency [steps]
  ; userreal1 - inner/vacuum dielectric [PS3]
  ; userreal2 - outer/solvent dielectric [PS3]
  ; userreal3 - Sheffield: A parameter 
  ;           - Andersen:  target temperature
  ;           - vvSD: target temperature [PS3]
  ; userreal4 - Sheffield B parameter 
  ;           - friction [1/ps] thus 91 would represent water at room temperature [PS3]
  
  userint1    = 1    ; update frequency for GB forces (should be set equal to nstlist)
  userreal1   = 1.   ; interior dielectric constant
  userreal2   = 80.  ; exterior dielectric constant
  userint4    = 1600  ; Still GB on (2^10 = 1024), ACE SA term in the GB/SA (2^9 = 512), VVSD (2^6 = 64)
  userint2    = 0    ; no logging, apply solvation forces
  userreal3   = 300. ; temperature for VVSD integrator (in K)
  userreal4   = 91.0 ; viscosity for VVSD integrator (in 1/ps)
  """
  
      # implicit solvation models     
      self.use_StillGB = True           # ; userint4 - BIT 10 - turns on Still GB [PS3]
      self.use_GB = False               # ; userint4 - BIT 1 - GB on/off [PS3]
      self.use_AGBNP = False            # ; userint4 - BIT 0 - AGBNP on/off
      self.use_Sheffield = False        # ; userint4 - BIT 5 - Sheffield on/off
      
      # solvent model options
      self.use_SASD = False    # SASD is velocity scaling based on solventg accessibility. Setting this option will modify GB, AGBNP
                                        # ; userint4 - BIT 2 - SASD for GB on/off [PS3]
                                        # ; userint4 - BIT 4 - SASD for AGBNP on/off 
      self.use_ACE_SA = True            # ; userint4 - BIT 9 - Still SA part f GB/SA on/off [PS3]
          
      # integration
      self.use_VVSD = True              # ; userint4 - BIT 6 - use the velocity Verlet version of the SD algorithm (vvSD) [PS3]
      self.use_Andersen = False         # ; userint4 - BIT 3 - Andersen on/off
      
      # Additional Parameters 
      self.GB_frequency = 1             # ; userint1 - GB frequency (used for GB and AGBNP) [PS3]  (should be set equal to nstlist)
      self.scaleCharges = False         # ; userint4 - BIT 7 - scale charges via sdelta  => gradually tuning off electrostatics (AGBNP & EdGB)
      self.use_GCMC = False             # ; userint4 - BIT 8 - turn on/off the GCMC part (Grand Canonical)  [not ready yet!]
      
      ### userint2
      self.GB_logging = False           # ; userint2 - BIT 1 - GB logging on/off
      self.AGBNP_logging = False        # ; userint2 - BIT 2 - AGBNP logging on/off
      self.DontApplySolvationForces = False  # ; usetint2 - BIT 3 - do NOT apply the solvation forces, this is helpful for comparisons
      
      ### userint3
      self.Andersen_frequency = 1       # ; userint2 - BIT 1 - GB logging on/off
  
      
      # Continuous parameters
      self.innerDielectric = 1.0
      self.outerDielectric = 80.0
      self.Sheffield_A = 1.0            # ; userreal3 - Sheffield A parameter 
      self.Sheffield_B = 1.0            # ; userreal4 - Sheffield B parameter 
      self.targetTemperature = 300.     # ; userreal3 - Andersen:  target temperature
                                        # ;           - vvSD: target temperature [PS3]
      self.viscosity = 91.              # ; userreal4 - friction [1/ps] thus 91 would represent water at room temperature [PS3]
  
  
      # user1_groups    
      self.q = {}          # ; A dictionary of key:value (strings) that can be given to MdpFile.setParameter()
     
  
      ### groups for visualization
      self.user1_grps = None            # show-lines, show-spacefill, show-sticks, ca-trace, or show-balls
      self.user1_residues = None        # This will be used by make_PS3_ndx()
      self.user2_grps = None
      self.user2_residues = None
      ### NOTE: asdfghjk
      # JDC sez:  To include more groups is difficult, as this requires overriding the default behavior of grompp
      # to prevent unnecessary groups from being included in the .tpr file. A modified grompp for the 3.1.4 EdGB
      # version of gromacs is available that allows specification of up to five groups, though this overrides the
      # ACC, FREEZE, and ORFIT capabilities (which we don't use anyway). The modified grompp requires recompilation
      # with this modified kernel/readir.c file.
      #
      # self.user1_grps = show-lines
      # self.user2_grps = show-spacefill
      # self.user3_grps = show-sticks ; requires modified grompp
      # self.user4_grps = ca-trace ; requires modified grompp
      # self.user5_grps = show-balls ; requires modified grompp
  
      ### Do the setup
      self.setup()                                      # derived classes can change the parms, if desired
      self.buildMdpParams()       # compiles a list of params to set in the mdpfile
      
      
    def setup(self):
      """This is the customizable parrt of derived classes, used to set parameters for a particular type of implicit solvation """
      pass
    
  
    def buildMdpParams(self, PS3_Visualization=False):
      """compiles a list of params to set in the mdpfile.
      Returns a dictionary of key:value (strings) that can be given to MdpFile.setParameter()"""
      
      params = {}
  
      # userint 1, 2 and 3
      params['userint1'] = '%d'%self.GB_frequency
      params['userint2'] = '%d'%( self.GB_logging*2 + self.AGBNP_logging*4 + self.DontApplySolvationForces*8  )
      if self.use_Andersen:
        params['userint3'] = '%d'%self.Andersen_frequency
      
      # userint4
      userint4 = 0
      if self.use_StillGB:
          userint4 += (2**10)
          if self.use_ACE_SA:
              userint4 += (2**9)
      elif self.use_GB:
          userint4 += (2**1) 
          if self.use_SASD:
              userint4 += (2**2)
      elif self.use_AGBNP: 
          userint4 += 1
          if self.use_SASD:
              userint4 += (2**4) 
      elif self.use_Sheffield:   
          userint4 += (2**5) 
  
      if self.use_VVSD:
          userint4 += (2**6)
      if self.use_Andersen:
          userint4 += (2**3)
          
      params['userint4'] = '%d'%userint4
          
      # userreals
      params['userreal1'] = '%3.2f'%self.innerDielectric
      params['userreal2'] = '%3.2f'%self.outerDielectric
      if self.use_Sheffield:
          params['userreal3'] = '%3.2f'%self.Sheffield_A
          params['userreal4'] = '%3.2f'%self.Sheffield_B
      elif self.use_Andersen or self.use_VVSD:
          params['userreal3'] = '%3.2f'%self.targetTemperature
      params['userreal4'] = '%3.2f'%self.viscosity    
    
      if PS3_Visualization: 
          # user1_grps and user2_grps for PS3 Advanced Visualization  
          if self.user1_grps != None:
              params['user1_grps'] = self.user1_grps
          if self.user2_grps != None:
              params['user2_grps'] = self.user2_grps
  
      self.mdpParams =  params
     
      return
      
      
  
  # Derived Classes
  
class StillGBOptions314(ImplicitOptions314):
    def setup(self):
        """Customize the options according to StillGB + VVSD, eprot = 1.0, esolv = 80.0"""
        pass   # This is the same as our base class.
    
class GBOptions314(ImplicitOptions314):
    def setup(self):
        """Customize the options according to GB + SASD + VVSD, eprot = 1.0, esolv = 80.0"""
        self.use_GB = True
        self.use_StillGB = False         
    
class AGBNPOptions314(ImplicitOptions314):
    def setup(self):
        """Customize the options according to AGBNP + SASD + VVSD, eprot = 1.0, esolv = 80.0"""
        self.use_GB = False
        self.use_StillGB = False         
        self.use_AGBNP = True      
    
class SheffieldOptions314(ImplicitOptions314):
    def setup(self):
        """Customize the options according to Sheffield + VVSD, eprot = 1.0, esolv = 80.0"""
        self.use_GB = False
        self.use_StillGB = False         
        self.use_AGBNP = False
        self.use_Sheffield = True
        
class SigmoidalDielectricOptions314(ImplicitOptions314):
    def setup(self):
        """Customize the options according to  VVSD, eprot = 1.0, esolv = 80.0"""
        self.use_GB = False
        self.use_StillGB = False         
        self.use_AGBNP = False
        self.use_Sheffield = False
 
# MdpFile.py exceptions
class MdpFileFormatError(Exception): pass


        
