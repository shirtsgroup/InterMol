import sys, os, string

# JDC: It makes things much simpler for other developers if there is only one class per file, and if the filename and class match.  That way, we don't have to hunt for particular classes in files.  The files are less cluttered that way.
# JDC: Classes should be named starting with a capital letter, such as "IonCalculator".  Maybe 'IonicEnvironmentCalculation' would be more appropriate.
class ionCalculator(object):
    # JDC: This class lacked a docstring, so I created one.
    """A utility object to calculate numbers, concentrations of ions.
    
    AUTHOR
      Vincent Voelz - primary author
      John Chodera - secondary cleanup, commenting

    REQUIREMENTS
      gromacs must be installed, and GMXLIB, GMXPATH, and MMTOOLSPATH set appropriately.
      
    EXAMPLES
      # Create a gromacs system from a given PDB file, using the 'ffamber99p' forcefield.
      import mmtools.gromacstools.system as system
      gromacsSystem = system.GromacsSystem('test.pdb' useff='ffamber99p')
      
      # Add ions to bring up to 150 mM NaCl salt conditions.
      gromacsSystem.setup.setSaltConditions('NaCl', 0.150)

      # Prepare a system, writing TPR and GRO files, and *mdrun* script in the current directory
      import os, sys
      gromacsSystem.prepare(outname='prod', os.path.abspath(os.curdir))      
    """        

    # JDC: 'useff' might be named something more clear, like 'forcefieldName', or perhaps 'ffname'.
    def __init__(self, useff):
        # JDC: Added documentation on what required argument 'useff' signified.
        """Initialize the ion calculator given the gromacs forcefield of choice.
        
        REQUIRED ARGUMENTS
          useff - the name of the forcefield to use in calling grompp (e.g. 'ffamber99p')
        """
        
        # Store the chosen forcefield.
        # JDC: Shouldn't we check here to see if this forcefield is supported, so we can throw an exception of not?
        self.useff = useff
        
        # All possible ion names that we know about and could occur in gromacs forcefield files.
        # JDC: Should this be moved to static class data?
        self.possibleIons = ['Ca','Cl','Na','Mg', 'K', 'Rb', 'CS', 'Li', 'Zn', 'Sr', 'Ba' ]

        # Create empty dictionary to contain all the ions and ion information provided by the forcefield.
        self.ions = {}    # a dictionary of Ion objects {'Cl', <Ion object> }
        
        # Read list of ions, their charges, and masses from gromacs forcefield files.
        # JDC: We should through a Exception of some sort if this fails.
        self.getIons()  
        
        return

    # JDC: Is this a private method, or part of the public API?
    # JDC: Should this use arguments/return data, rather than operating on just internal object data?
    # JDC: Does this need to be instantiated separately for each object, or should it be a Singleton, only performed once for each available forcefield upon first class startup?
    def getIons(self):
        """Get atomic mass constants for all ion types listed in self.possibleIons from the gromacs forcefield files, storing them in self.ions.        
        
        """

        # Read the contents of the .rtp (residue topology) file for this forcefield.
        rtpfile = os.path.join(os.environ['GMXLIB'],(self.useff+'.rtp'))
        fin = open(rtpfile,'r')
        rtplines  = fin.readlines()
        fin.close()
        
        # Build up dictionary of ions present in .rtp file.
        for line in rtplines:
            fields = line.split()
            # JDC: Instead of using fields[0], fields[1], fields[2], the line should be parsed into meaningfully-named tokens.
            # This makes the comprehension of what this codes much more straightforward.
            if len(fields) > 0:
                if (fields[0] in self.possibleIons):                
                    # JDC: This idiom of parsing into named variables and then processing is much more readable, if more verbose, to other coders.
                    # Parse fields.
                    # Example format:
                    # Ca     amber99_15   2.00000     1
                    atomName = fields[0]
                    atomType = fields[1]
                    charge = float(fields[2])
                    # Store this ion in object ion dictionary.
                    self.ions[atomName] = Ion(atomName, atomType, charge)  
                        
        # Read contents of .atp (atom type masses) file for this forcefield.
        atpfile = os.path.join(os.environ['GMXLIB'],(self.useff+'.atp'))
        fin = open(atpfile,'r')
        atplines  = fin.readlines()
        fin.close()
        
        # Parse .atp file to determine atomic masses.
        for line in atplines:
            fields = line.split()
            if len(fields) > 0:
                # find the right ion
                for key in self.ions.keys():
                    # JDC: This could be more clear.
                    if self.ions[key].atomType == fields[0]:  
                        self.ions[key].atomicMass = float(fields[1])
                        self.ions[key].atomDescription = string.joinfields(fields[3:])

    # JDC: Should use Units class for concentration.  There is a 'Units.M' or 'Units.Molar' unit.
    def molarityToFraction(self, concentration, solvent='water'):
        """Compute the number fraction of ion molecules given a Molar concentration.
        
        REQUIRED ARGUMENTS
          concentration - the concentration of ions (Molar)
          
        OPTIONAL ARGUMENTS
          solvent - the name of the solvent ('water' is only solvent currently supported, defaults to 'water')
        
        RETURNS
          ???????????
        """

        # JDC: Should raise an exception if unrecognized solvent is passed.
        if solvent != 'water':
            print 'Only solvent=\'water\' is supported.  Exiting.'
            sys.exit(1)

        # JDC: You should use Units class for computations like this.
        # The calculation would be something like the following (if I understand what you're trying to do):
        #
        # import mmtools.utilities.Units as Units
        # solventMass = 18.016 * Units.g / Units.mol
        # solventDensity = 1.0 * Units.g / Units.cm**3
        # solventConcentration = solventDensity / solventMass
        # ionFraction = ionConcentration / solventConcentration        
        #
        # See how easy that is!!!
        waterGramPerMol = 18.016
        waterMolPerGram = 1./waterGramPerMol
        waterGramPerLiter = 1000.0
        waterMolPerLiter = waterMolPerGram * waterGramPerLiter
        
        # JDC: I scratched my head for a while and couldn't figure out what the quantity in parenthesis represented.  Maybe I'm just tired, but it would be much less confusing if you assigned it to a variable like 'numberOfIonsInBox' and then returned *this* variable instead.        
        return (concentration/waterMolPerLiter)

        
# JDC: This class should be in a separate file.
class Ion(object):
    """A monoatomic ion.    
    
    Stores atom name, forcefield atom type, charge, a textual description, and the atomic mass of a single-atom ion.

    """
    
    def __init__(self, atomName, atomType, charge):
        """Initialize the data structures.

        REQUIRED ARGUMENTS
          atomName - the ion name
          atomType - the forcefield atom type
          charge - the partial atomic charge of the ion (in units of electron charge)
          
        """
        
        # Initialize object data.              # JDC: Every variable should be described the first time it is defined.
        self.atomName = atomName               # the atom name (e.g. 'Cl')
        self.atomType = atomType               # the forcefield atom type (e.g. 'amber99_30')
        self.charge = charge                   # the partial atomic charge for the ion (in units of electron charge)
        self.atomDescription = None            # a textual description of the ion (e.g. "chloride ion")
        self.atomicMass = None                 # the mass of the ion (in amu)

        return
    
    def __repr__(self):
        """Self-documented string representation of the internal data structures.
        
        """
        
        # Construct string representation
        # JDC: Descriptions are much clearer to construct line-by-line.
        tmpstr = '%s:\n' % self.atomName
        tmpstr += '\tatomtype: %s\n' % self.atomType
        tmpstr += '\tcharge %2.1f\n' % self.charge # JDC: Is this a sufficient number of decimal places?
        tmpstr += '\tdescription = %s\n' % self.atomDescription
        tmpstr += '\tatomic mass %2.1f\n' % self.atomicMass # JDC: Is this a sufficint number of decimal places?

        # Return string representation.
        return tmpstr
    
        
        
