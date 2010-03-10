import sys,os,tempfile, string

from Ion import *
import mmtools.utilities.Units as Units

class IonEnvironmentCalculation(object):
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

    def __init__(self, system):
        """Initialize the ion calculator given the gromacs forcefield of choice.
        
        REQUIRED ARGUMENTS
          system - a SimulationPreparation.SystemSetup() object containing salt concentration, forcefield information, etc. 
        """
        
        # Store all the system info
        self.system = system  # a SystemSetup() object containing salt concentration, forcefields, etc,
        
        # All possible ion names that we know about and could occur in gromacs forcefield files.
        # JDC: Should this be moved to static class data?
        self.possibleIons = ['Ca','Cl','Na','Mg', 'K', 'Rb', 'CS', 'Li', 'Zn', 'Sr', 'Ba' ]

        # Read list of ions, their charges, and masses from gromacs forcefield files, and
        # create a dictionary to contain all the ions and ion information provided by the forcefield.
        
        #try:
        if (1):
            self.ions = self.getIons()    # a dictionary of Ion objects {'Cl', <Ion object> }
        #except:
        #   raise IonParsingError
        
        return


    def getIons(self):
        """Get atomic mass constants for all ion types listed in self.possibleIons from the gromacs forcefield files
        RETURNS
          a dictionary of Ion() objects, indexed by the ion name.        
        """

        # Initialize a dictionary to return
        iondict = {}
        
        # Read the contents of the .rtp (residue topology) file for this forcefield.
        rtpfile = os.path.join(os.environ['GMXLIB'],(self.system.forcefield+'.rtp'))
        fin = open(rtpfile,'r')
        rtplines  = fin.readlines()
        fin.close()
        
        # Build up dictionary of ions present in .rtp file.
        for line in rtplines:
            fields = line.split()
            if len(fields) > 0:
                if (fields[0] in self.possibleIons):                
                    # Parse fields.
                    # Example format:
                    # Ca     amber99_15   2.00000     1
                    atomName, atomType, charge = fields[0], fields[1], float(fields[2])
                    # Store this ion in object ion dictionary.
                    iondict[atomName] = Ion(atomName, atomType, charge)  
                        
        # Read contents of .atp (atom type masses) file for this forcefield.
        atpfile = os.path.join(os.environ['GMXLIB'],(self.system.forcefield+'.atp'))
        fin = open(atpfile,'r')
        atplines  = fin.readlines()
        fin.close()
        
        # Right now the Ion() objects in the iondict do not have an atomicMass or atomDescription
        # We will parse .atp file to determine these
        for line in atplines:
            # Parse fields
            # Example format:
            # amber94_30        35.45000	; Cl-
            fields = line.split()
            if len(fields) > 0:
                # For each line of the *.atp file, try to find an Ion() object that matches it.
                for key in iondict.keys():                    
                    if iondict[key].atomType == fields[0]:  
                        iondict[key].atomicMass = float(fields[1])
                        iondict[key].atomDescription = string.joinfields(fields[3:])

        # Return a dictionary of Ion() objects
        return iondict


    def molarityToNumberFraction(self, ionConcentration, solvent='water'):
        """Compute the number fraction of ion molecules given a Molar concentration.
        
        REQUIRED ARGUMENTS
          ionConcentration - the concentration of ions (Molar)
          
        OPTIONAL ARGUMENTS
          solvent - the name of the solvent ('water' is only solvent currently supported, defaults to 'water')
        
        RETURNS
          numberFraction - the number fraction (0. < numberFraction < 1.) representing the ratio of total numbers of 
                           salt molecules to water molecules 
        """

        if solvent != 'water':
            raise NonWaterSolventNotSupported
        
        solventMass = 18.016 * Units.g / Units.mol
        solventDensity = 1.0 * Units.g / (0.001 * Units.l)   # water is 1 g/mL
        solventConcentration = solventDensity / solventMass
        ionFraction = (ionConcentration * Units.mol / Units.l) / solventConcentration        
        
        print 'solventMass = 18.016 * Units.g / Units.mol = ', solventMass
        print 'solventDensity = 1.0 * Units.g / Units.cm**3 = ', solventDensity
        print 'solventConcentration = solventDensity / solventMass =', solventConcentration
        print 'ionFraction = ionConcentration / solventConcentration = ', ionFraction
        
        return ionFraction


# Exceptions
class IonParsingError(Exception):
    print "There was an error parsing the information about each ion from the forcefield files."
        
class NonWaterSolventNotSupported(Exception):
    print "Currently, water is the only supported solvent."
        
