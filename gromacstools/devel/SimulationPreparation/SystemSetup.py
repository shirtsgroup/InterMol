
# HISTORY
#
# 2/11/2010 VAV -- rewrote as gromacstools2.0 object
#
# 12/21/07 JDC - added salt information for CaCl2 for kinase project
# 06/26/07 GRB - added parameters/functions for using absolute box size

import sys, os, tempfile, string

class SystemSetup(object):

    def __init__(self, infile=None): 
	"""A class to keep track of parameters important to the system setup."""
	
	
	# forcefield parameters
	self.forcefieldCodes = self.getForcefieldCodes()
	    # ff dictionary {'ffamber03':'3'} e.g. is stored in self.ForceFieldCodes
	    # The index comes from the FF.dat file in $GMXLIB
	    
	self.forcefield = 'ffamber03'      # string to specify forcefield (default is 'ffamber99p').  'ffamber03' also supported
	self.gmx_version = '3.1.4'         # gmx version.  Supported are '3.1.4' and '3.3'
	
	
	# A list of paths to look in for mdpfiles  
	self.mdpPaths = [ os.path.join( os.environ['MMTOOLSPATH'], 'gromacstools/devel/SimulationPreparation/mdp') ]
	
	
	# editconf pararameters for periodic box conditions and solvent boxes 
	self.boxType = 'octahedron'        #  -bt   enum   tric  Box type for -box and -d: 'tric', 'cubic', 'dodecahedron' or 'octahedron'
	self.boxSoluteDistance = '0.9'     #  -d   real      0  Distance between the solute and the box (in nm)
	self.useAbsBoxSize = False         # whether use absolute or relative box size
	self.absBoxSize = '3.0'            # absolute box size  (in nm)
	
	# salt conditions (defaults)
	self.salt = 'NaCl'                #  Supported:  'NaCl' or 'SodiumChloride'
	self.saltconc = 0.050             #  Molar salt concentration
	self.positiveIonName = 'Na'
	self.negativeIonName = 'Cl'
	self.positiveIonCharge = 1
	self.negativeIonCharge = -1 
	self.positiveStoichiometry = 1
	self.negativeStoichiometry = 1
    
	# water box file used for explicit solvation
	self.WaterBoxFile = os.path.join(os.environ['GMXLIB'],'ffamber_tip3p.gro')
	
	# cosolvent conditions
	# NOTE:  Any cosolvent simulation must also have pre-constructed water box that will be used for solvation
	self.cosolventName = None
	self.cosolventConcentration = 0.0
	self.cosolventWaterBoxFile = None
	
	
      
	# setup from file?  THIS is still a DUMMY function....
	if infile != None:
	    self.read(infile)
    
  
    def getForcefieldCodes(self):
	"""Find the numeric code for specifying the forcefield listed in $GMXLIB/FF.dat."""
	
	codes = {}  # a dictionary {'ffamber03':'3'}  of the forcefield name and the (string) code to specifiy it 
	
	try:
	    fin = open(os.path.join(os.environ['GMXLIB'],'FF.dat'),'r')
	    
	except GMXForcefieldCodeError:
	    print "File FF.dat cannot be found in GMXLIB=%s   ..... \nExiting....\n"%self.files.GMXLIB
	    sys.exit(1)
	    
	lines = fin.readlines()
	numff = int(lines.pop(0)) 
	for i in range(0,numff):
	    fields = lines[i].split()
	    key = fields.pop(0) #DLM: Looks like there is a bug in indentation here.
	    codes[ key ] = str(i)
	    
	return codes

    
    def setForcefield(self, forcefield):
	"""Sets the forcefield to the given name.  
	
	NOTES
	The forcefield must be one of the names matching those listed in $GMXLIB/FF.dat.
	This method will cross-check the FF.dat list and raise an exception if it can't
	find it.  ***If you have recently installed a forcefield parameter file, and get
	this error, you must add your new forcefield to the $GMXLIB/FF.dat file.
	"""
	
	if not (forcefield in self.forcefieldCodes.keys()):
	    raise ForcefieldNotFoundError
	
	self.forcefield = forcefield
	
	  
    def setSaltConditions(self, saltname, saltconcentration): 
	"""Set ion names and stoichiometry for the chosen salt."""
	if saltname == 'NaCl' or saltname=='SodiumChloride':
	  self.salt = 'NaCl'                #  Supported:  'NaCl' or 'SodiumChloride'
	  self.saltconc = saltconcentration            #  Molar salt concentration
	  self.positiveIonName = 'Na'
	  self.negativeIonName = 'Cl'
	  self.positiveIonCharge = 1
	  self.negativeIonCharge = -1
	  self.positiveIonStoichiometry = 1
	  self.negativIoneStoichiometry = 1
    
	elif saltname == 'CaCl2' or saltname=='CalciumChloride':
	  self.salt = 'CaCl2'                #  Supported:  'CaCl2' or 'CalciumChloride'
	  self.saltconc = saltconcentration            #  Molar salt concentration
	  self.positiveIonName = 'Ca'
	  self.negativeIonName = 'Cl'
	  self.positiveIonCharge = 2
	  self.negativeIonCharge = -1
	  self.positiveIonStoichiometry = 2
	  self.negativeIonStoichiometry = 1
    
	elif saltname == 'MgCl2' or saltname=='MagnesiumChloride':
	  self.salt = 'MgCl2'                #  Supported:  'MgCl2' or 'MagnesiumChloride'
	  self.saltconc = saltconcentration            #  Molar salt concentration
	  self.positiveIonName = 'Mg'
	  self.negativeIonName = 'Cl'
	  self.positiveIonCharge = 2
	  self.negativeIonCharge = -1
	  self.positiveIonStoichiometry = 2
	  self.negativeIonStoichiometry = 1
    
    
	else:
	    print 'salt type', saltname, 'not supported!  Try \'NaCl\' or \'CaCl2\' or \'MgCl2\'.  Exiting.'
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
	pass
    
    def write(self, filename):
	"""Write parameters to a GromacsSystemSetup text file."""
	pass
  
    def __repr__(self): 
	"""Return a string listing parameters important to the system setup."""
    
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
  

# SystemSetup.py exceptions
class ForcefieldNotFoundError(Exception):
    
    print """
  NOTE: The forcefield must be one of the names matching those listed in $GMXLIB/FF.dat.
	This method will cross-check the FF.dat list and raise an exception if it can't
	find it.  ***If you have recently installed a forcefield parameter file, and get
	this error, you must add your new forcefield to the $GMXLIB/FF.dat file.
	"""
    
    
    
class SystemSetupError(Exception): pass
