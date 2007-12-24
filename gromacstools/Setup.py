# Change Log:
# 06/26/07 GRB - added parameters/functions for using absolute box size
# 12/21/07 JDC - added salt information for CaCl2 for kinase project

import sys,os,tempfile, string

class Setup(object):

  def __init__(self, infile=None): 
    """A class to keep track of parameters important to the system setup."""

    # editconf parms
    self.boxType = 'octahedron'            #  -bt   enum   tric  Box type for -box and -d: 'tric', 'cubic', 'dodecahedron' or 'octahedron'
    self.boxSoluteDistance = '0.9'    #  -d   real      0  Distance between the solute and the box
    self.useAbsBoxSize = False      # whether use absolute or relative box size
    self.absBoxSize = '3.0'
    
    # salt conditions (defaults)
    self.salt = 'NaCl'                #  Supported:  'NaCl' or 'SodiumChloride'
    self.saltconc = 0.050             #  Molar salt concentration
    self.positiveIonName = 'Na'
    self.negativeIonName = 'Cl'
    self.positiveIonCharge = 1
    self.negativeIonCharge = -1 
    self.positiveStoichiometry = 1
    self.negativeStoichiometry = 1

    # cosolvent conditions
    self.cosolvent = None
    self.cosolventconc = 0.0
   
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

    elif saltname == 'CaCl2' or saltname=='calcium chloride':
      self.salt = 'CaCl2'                #  Supported:  'CaCl2' or 'CalciumChloride'
      self.saltconc = saltconcentration            #  Molar salt concentration
      self.positiveIonName = 'Ca'
      self.negativeIonName = 'Cl'
      self.positiveIonCharge = 2
      self.negativeIonCharge = -1
      self.positiveIonStoichiometry = 2
      self.negativIoneStoichiometry = 1

    else:
        print 'salt type', salt, 'not supported!  Try \'NaCl\' or \'CaCl2\'.  Exiting.'
	sys.exit(1)
    return None
 
  def setCosolventConditions(self, cosolventname, concentration):
    """Set ion names and stoichiometry for the chosen salt."""
    if cosolventname == 'gnd':
      self.cosolvent = 'gnd'                #  Supported:  'NaCl' or 'SodiumChloride'
      self.cosolventconc = concentration            #  Molar salt concentration

    else:
        print 'salt type', salt, 'not supported!  Try \'NaCl\'.  Exiting.'
        sys.exit(1)
    return None
 
  def set_boxSoluteDistance(self, boxSoluteDistance):
    """Set the distance between the solute (protein) and the periodic box.
    INPUT
    boxSoluteDistance - distance in nm"""
    self.boxSoluteDistance = str(boxSoluteDistance)
    return None
 
  def set_boxType(self, boxType):
    """Set the distance between the solute (protein) and the periodic box.
    INPUT
    boxType - can be either ''tric', 'cubic', 'dodecahedron' or 'octahedron'."""
    self.boxType = str(boxType)
    return None
 
  def setUseAbsBoxSize(self, useAbsBoxSize):
    """Set whether or not to use an absolute box size.
    INPUT
    useAbsBoxSize - bool"""
    self.useAbsBoxSize = useAbsBoxSize

  def setAbsBoxSize(self, size):
    """Set the absolute size of the box.  Size is applied to all dimensions.
    INPUT
    size - box size in nm, string""" 
    self.absBoxSize = size

  def read(self, filename): 
    """Read parameters in from a GromacsSystemSetup text file."""
    return None
  
  def write(self, filename):
    """Write parameters to a GromacsSystemSetup text file."""
    return None

  def __repr__(self): 
    """A class to keep track of parameters important to the system setup."""

    # editconf parms
    outstr = ''
    outstr = outstr + '%-30s %s\n'%( 'self.boxType',repr(self.boxType) )
    outstr = outstr + '%-30s %s\n'%( 'self.boxSoluteDistance', repr(self.boxSoluteDistance ) )
    outstr = outstr + '%-30s %s\n'%( 'self.useAbsBoxSize', repr(self.useAbsBoxSize) )
    outstr = outstr + '%-30s %s\n'%( 'self.absBoxSize', repr(self.absBoxSize) )
    
    # salt conditions
    outstr = outstr + '%-30s %s\n'%(  'self.salt', repr(self.salt  ) )
    outstr = outstr + '%-30s %s\n'%(  'self.saltconc', repr(self.saltconc  ) )
    outstr = outstr + '%-30s %s\n'%(  'self.positiveIonName', repr(self.positiveIonName  ) )
    outstr = outstr + '%-30s %s\n'%(  'self.negativeIonName', repr(self.negativeIonName  ) )
    outstr = outstr + '%-30s %s\n'%(  'self.positiveStoichiometry', repr(self.positiveStoichiometry  ) )
    outstr = outstr + '%-30s %s\n'%(  'self.negativeStoichiometry', repr(self.negativeStoichiometry  ) )
    
    return outstr

