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
    
        
        
