# Change Log:
# 06/26/07 GRB - added parameters/functions for using absolute box size

import sys,os,tempfile, string

class Setup(object):

  def __init__(self, infile=None): 
    """A class to keep track of parameters important to the system setup."""

    # editconf parms
    self.boxType = 'octahedron'            #  -bt   enum   tric  Box type for -box and -d: tric, cubic, dodecahedron or octahedron
    self.boxSoluteDistance = '0.9'    #  -d   real      0  Distance between the solute and the box
    self.useAbsBoxSize = False      # whether use absolute or relative box size
    self.absBoxSize = '3.0'
    
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
  
  def set_boxSoluteDistance(self, boxSoluteDistance):
    """Set the distance between the solute (protein) and the periodic box.
    INPUT
    boxSoluteDistance - distance in nm"""
    self.boxSoluteDistance = str(boxSoluteDistance)
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
  