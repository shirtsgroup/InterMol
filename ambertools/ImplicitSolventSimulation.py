# hacks.py 
#
# An protein topology is constructed in the extended conformation in AMBER LEaP, and an implicit solvent simulation
# is conducted at high temperature with a Langevin thermostat to ensure rapid sampling of thermally denatured states.
#
# We start with a PDB file that comes from MCCE with the appropriate AMBER residues names already substituted.
#
# John D. Chodera
# 6 Jun 2007
#
# TODO:
# - Determine how to use substitution with both class methods and local variables.

import mmtools.utilities.Units as Units

def writeFileContents(filename, contents):
    """Write the given contents to a file.
    
    ARGUMENTS
        filename - the name of the file to write to
        contents - the contents of the file to be written.
    """
    
    outfile = open(filename, 'w')
    outfile.write(contents)
    outfile.close()
    
    return

class ImplicitSolventSimulation(object):
    
    def getDefaultSanderOptions(self):
        """
        Create a dictionary with default options for AMBER 'sander' for implicit solvent simulation.
        """
        # Create a new dictionary to store sander options.        
        options = dict()
        # dynamics options
        options['imin'] = 0       # run dynamics
        options['nstlim'] = 500   # number of steps to run dynamics
        options['ntwr'] = 500     # number of steps betwee restart files
        options['dt'] = 0.002     # 2 fs timestep
        options['nrespa'] = 2     # period for GB force update (AMBER manual max recommended for GB)
        # GB options
        options['igb'] = 1        # GB model selection
        options['gbsa'] = 1       # flag to add surface area penalty
        options['surften'] = 0.005 # surface tension in cal/mol (AMBER manual default)
        options['cut'] = 16.0     # nonbonded cutoff in A (AMBER manual minimum recommended)
        options['rgbmax'] = 10.0  # GB radius calc cutoff in A (AMBER manual minimum recommended)
        # constraints
        options['ntc'] = 2        # constrain bonds to hydrogen
        options['ntf'] = 2        # do not compute bond energy terms for constrained bonds to hydrogen
        options['tol'] = 1.0e-8   # recommended SHAKE tolerance
        # periodic boundary conditions
        options['ntb'] = 0        # non-periodic
        options['iwrap'] = 0      # no wrap
        options['nscm'] = 500     # remove COM translation every 500 steps
        # output options
        options['ioutfm'] = 0     # output format text
        options['ntwx'] = 500     # coordinate output interval in steps
        options['ntwe'] = 500     # energy output interval in stepes
        options['ntpr'] = 500     # interval for printing output to terminal
        options['ntave'] = 5000   # interval for printing average energies
        # random number seed
        import random
        random.seed()
        options['ig'] = random.randint(0, 0x0000ffffL)

        return options
    
    def createSanderInputFile(self, options, title = ''):
        """Create sander mdin.sander.in file contents from key-value pairs.
        
        ARGUMENTS
            options (dictionary) - key-value pairs to be placed in &cntrl...&end block.
            
        OPTIONAL ARGUMENTS
            title (string) - one-line title string (default: "")
        """
        
        mdin = "" # string to hold contents of mdin file
        
        # Title
        mdin += "%s\n" % title
        
        # &cntrl block begin
        mdin += " &cntrl\n  "        
        
        # key-value pairs, limit number per line
        keys_this_line = 0
        keys_per_line = 4
        for (key, value) in options.items():
            if (keys_this_line > 0): mdin += ", "
            mdin += "%s = %s" % (key, options[key].__str__())
            keys_this_line += 1
            if (keys_this_line == keys_per_line):
                mdin += "\n  "
                keys_this_line = 0
        
        # &cntrl block end
        mdin += "\n &end\n"

        # Return contents of mdin file.
        return mdin
    
    def __init__(self, initial_pdb_filename = None, sequence = None, leaprc_filename = "leaprc.ff03", gbradii = "amber6", workdir = ".", amberhome = None, debug = False):
        """\
        Setup an implicit solvent simulation using AMBER with LEaP.
        
        REQUIRED ARGUMENTS
          One of the following is required:
            sequence - Protein sequence (three-letter code, separated by spaces, all caps) used to generate protein in extended conformation.
            initial_pdb_filename - Name of PDB file with AMBER-compatible three-letter residues names used to generate protein in given conformation.
          If both are specified, the given sequence is used instead of the sequence present initial_pdb_filename.
          In both cases, hydrogen atoms are stripped out.
        
        OPTIONAL ARGUMENTS
            leaprc_filename - leaprc file used to select forcefield (default: leaprc.ff03)
            gbradii - GB radii set (default: "amber6")        
            sander_options - dictionary of override options for sander (default: None)
            workdir - name of directory to set up AMBER simulations in (default: .)
            amberhome - location of AMBER root installation (default: taken from AMBERHOME environment variable)
            debug - print extra debug information (default: False)
        """

        ### Check arguments.

        if not (sequence or initial_pdb_filename):
            raise "Either 'sequence' or 'initial_pdb_filename' must be specified."

        if initial_pdb_filename:
            # get absolute path for specified initial pdb filename
            import os.path
            initial_pdb_filename = os.path.abspath(initial_pdb_filename)            

        ### Check arguments and store them

        import os.path
        if not os.path.exists(workdir):
            raise "Specified work directory %s does not exist." % workdir        
        self.workdir = workdir
        self.debug = debug        

        ### Determine location of AMBER
        
        import os
        if amberhome != None:
            # Use provided argument to set the location of AMBER.
            self.amberhome = amberhome
        elif os.environ.has_key("AMBERHOME"):
            # See if we can find AMBERHOME environment variable
            self.amberhome = os.environ["AMBERHOME"]
        else:
            # No AMBERHOME defined
            raise "No amberhome specified."
        
        self.leap = os.path.join(self.amberhome, "exe/tleap")
        self.sander = os.path.join(self.amberhome, "exe/sander")
        self.ambpdb = os.path.join(self.amberhome, "exe/ambpdb")
        # TODO: Ensure tleap and sander both exist.
        if not os.path.exists(self.leap):
            raise "Could not find 'tleap' executable at %s" % self.leap
        if not os.path.exists(self.sander):
            raise "Could not find 'sander' executable at %s" % self.sander
        if not os.path.exists(self.ambpdb):
            raise "Could not find 'ambpdb' executable at %s" % self.ambpdb

        ### Copy PDB file excluding hydrogen atoms.
        if initial_pdb_filename:
            # Read the PDB file into memory.
            initial_pdb_file = open(initial_pdb_filename, 'r')            
            lines = initial_pdb_file.readlines()            
            initial_pdb_file.close()
            # Write non-hydrogen atoms.
            import os.path
            outfile = open(os.path.join(self.workdir, 'stripped.pdb'), 'w')
            for line in lines:
                if (line[0:6] == 'ATOM  ') and (line[13:14] != 'H') and (line[13:15] != 'OC') and (not (line[17:21] == 'CILE' and line[13:15] == 'CD')):
                    outfile.write(line)
            outfile.close()

        ### Set up the system in LEaP
        
        # Create a LEaP input file.
        aux_gb_leap_commands = ""
        if gbradii == "amber6":
            # Append AMBER6 readii directive.
            aux_gb_leap_commands += "set default PBradii amber6\n\n"
        elif gbradii == "amber7":
            # AMBER7 radii are used by default.
            pass
        else:
            raise "Unrecognized gbradii argument '%s'." % gbradii

        if sequence and not initial_pdb_filename:
            # Build the protein in the extended conformation from the specified sequence.
            create_system_leap_commands = """\
# create a sequence in the extended conformation
system = sequence { %(sequence)s }
""" % vars()
        elif sequence and initial_pdb_filename:
            # Build the protein from the given initial PDB file using given sequence.
            create_system_leap_commands = """\
# load the sequence and initial conformation from the given PDB filename
system = loadPdbUsingSeq stripped.pdb { %(sequence)s }
""" % vars()
        elif initial_pdb_filename and not sequence:
            # Instruct LEaP to build the protein from the given initial PDB file.            
            create_system_leap_commands = """\
system = loadPdb stripped.pdb
""" % vars()
            
        date = ""        
        
        leap_in = """\
# Generate sequence in extended conformation for simulation in implicit solvent.
# Automatically generated by ambertools.hacks at %(date)s.

# load in amber force field
source %(leaprc_filename)s

%(aux_gb_leap_commands)s

%(create_system_leap_commands)s

# check unit before writing
check system

# write files
saveAmberParm system prmtop initial.crd

# exit
quit
""" % vars()
        writeFileContents(os.path.join(workdir, "leap.in"), leap_in)
            
        import commands
        command = "cd %(workdir)s; %(leap)s -f leap.in" % { 'workdir' : self.workdir, 'leap' : self.leap }
        output = commands.getoutput(command)
        print "> %s" % command
        print output    

        # Convert to PDB
        command = "cd %(workdir)s ; cat initial.crd | %(ambpdb)s -p prmtop > initial.pdb" % { 'workdir' : self.workdir, 'ambpdb' : self.ambpdb }
        print "Executing > %s" % command
        output = commands.getoutput(command)
        
        return

    def minimize(self, useropts = None):
        """\
        Minimize the system in preparation for dynamics.    
    
        OPTIONAL ARGUMENTS
            useropts (dictionary) - additional user options to override defaults in sander mdin file
            
        RETURNS
            pdb_filename - the absolute pathname to the PDB file generated by minimization
        """
        
        # Create default options and override with user-supplied options.
        options = self.getDefaultSanderOptions()
        options['imin'] = 1 # minimize
        options['maxcyc'] = 1000 # 1000 steps of minimization
        options['ncyc'] = 500 # 500 initial steps of gradient descent    
        options['nrespa'] = 1
        options['dx0'] = 1.0e-6
        if useropts != None:
            options.update(useropts) # override options with user options
        
        # Create mdin file
        import os.path
        mdin = self.createSanderInputFile(options)
        writeFileContents(os.path.join(self.workdir,"min.sander.in"), mdin)    
        
        # Minimize the system.    
        import commands
        command = "cd %(workdir)s ; %(sander)s -O -i min.sander.in -o min.sander.out -p prmtop -c initial.crd -r minimized.crd -inf mdinfo" % { 'workdir' : self.workdir, 'sander' : self.sander }
        print "Executing > %s" % command
        output = commands.getoutput(command)
        print output    

        # Convert to PDB
        command = "cd %(workdir)s ; cat minimized.crd | %(ambpdb)s -p prmtop > minimized.pdb" % { 'workdir' : self.workdir, 'ambpdb' : self.ambpdb }
        print "Executing > %s" % command
        output = commands.getoutput(command)
                
        # Return absolute pathname to minimized PDB file.
        pdb_filename = os.path.join(self.workdir, 'minimized.pdb')
        return pdb_filename
    
    def dynamics(self, temperature = 300.0 * Units.K, duration = 1.0 * Units.ns, useropts = None):
        """\
        Run dynamics at the desired temperature using Langevin integrator.
        
        If the dynamics have already been run and ther e
        
        OPTIONAL ARGUMENTS
            temperature [units:temperature]- the temperature of the heat bath, in K (default: 300 * Units.K)
            duration [units:time] - the simulation duration (default: 1.0 * Units.ns)
            useropts (dictionary) - additional user options to override defaults in sander mdin file

        RETURNS
            pdb_filename - the absolute pathname to the PDB file generated by minimization
        """

        # Create default options and override with user-supplied options.
        options = self.getDefaultSanderOptions()
        options['temp0'] = temperature / Units.K # temperature of heat bath
        options['ntt'] = 3 # Langevin thermostat
        options['gamma_ln'] = 50.0 # viscosity (in 1/ps)
        if useropts != None:
            options.update(useropts) # override options with user options

        # Set simulation duration.
        import math
        options['nstlim'] = int(math.ceil( duration / (options['dt'] * Units.ps) )) # number of steps
              
        # Resume from previous run if restart file is present.
        import os.path
        if os.path.exists(os.path.join(self.workdir, 'md.crd')):
            initial_crd_filename = 'md.crd'            
            options['ntx'] = 5 # resume from previous velocities and coordinates
            options['irest'] = 1 # restart from given conformation
            initial_crd = 'md.crd'        
        else:
            initial_crd_filename = 'minimized.crd' # restart from minimized conformation
            options['tempi'] = temperature / Units.K # initial temperature for velocity reassignment    
        
        # Create mdin file
        import os.path        
        mdin = self.createSanderInputFile(options)
        writeFileContents(os.path.join(self.workdir,"md.sander.in"), mdin)
        
        # Run dynamics
        import commands
        command = "cd %(workdir)s; %(sander)s -O -i md.sander.in -o md.sander.out -p prmtop -c %(initial_crd)s -r md.crd -x md.trj -e md.ene -inf mdinfo" % { 'workdir' : self.workdir, 'sander' : self.sander, 'initial_crd' : initial_crd_filename }
        print "Executing > %s" % command
        output = commands.getoutput(command)
        print output    

        # Return absolute pathname to PDB version of current restart file.
        pdb_filename = os.path.join(self.workdir, 'md.pdb')
        return pdb_filename
    
    
