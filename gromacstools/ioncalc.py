import sys, os, string

class ionCalculator(object):
	
    def __init__(self, useff):
        """A utility object to calculate numbers, concentrations of ions"""
	
	self.useff = useff
	self.possibleIons = ['Ca','Cl','Na','Mg', 'K', 'Rb', 'CS', 'Li', 'Zn', 'Sr', 'Ba' ]

	self.ions = {}    # a dictionary of Ion objects {'Cl', <Ion object> }
        self.getIons()	
	
	return
		
	
    def getIons(self):
        """get atomic mass constants from the forcefields files"""

        rtpfile = os.path.join(os.environ['GMXLIB'],(self.useff+'.rtp'))
        fin = open(rtpfile,'r')
        rtplines  = fin.readlines()
        fin.close()
	
	for line in rtplines:
	    fields = line.split()
	    if len(fields) > 0:
	      if (fields[0] in self.possibleIons):
		self.ions[fields[0]] = Ion(fields[0], fields[1], float(fields[2]))  #Ca     amber99_15   2.00000     1
			
        atpfile = os.path.join(os.environ['GMXLIB'],(self.useff+'.atp'))
        fin = open(atpfile,'r')
        atplines  = fin.readlines()
        fin.close()
	
	for line in atplines:
	    fields = line.split()
	    if len(fields) > 0:
		# find the right ion
		for key in self.ions.keys():
		    if self.ions[key].atomType == fields[0]:  
			self.ions[key].atomicMass = float(fields[1])
			self.ions[key].atomDescription = string.joinfields(fields[3:])


    def molarityToFraction(self, concentration, solvent='water'):
        """Computes the number fraction of ion molecules given a Molar concentration"""

        if solvent != 'water':
	    print 'Only solvent=\'water\' is supported.  Exiting.'
	    sys.exit(1)

            
	waterGramPerMol = 18.016
	waterMolPerGram = 1./waterGramPerMol
	waterGramPerLiter = 1000.0
	waterMolPerLiter = waterMolPerGram * waterGramPerLiter
		
        return (concentration/waterMolPerLiter)

	

	
class Ion(object):
    """A data structure for an ion"""
    
    def __init__(self, atomname, atomtype, charge):
	"""Initialize the data structures"""
	
	self.atomName = atomname               #'Cl', e.g.
	self.atomType = atomtype               #'amber99_30', e.g.
	self.charge = charge
	self.atomDescription = None
	self.atomicMass = None
	
    def __repr__(self):
	"""string representation of the data structures"""
	
	tmpstr = '%s:\n\tatomtype: %s\n\tcharge: %2.1f\n\tdescription: %s\n\tatomic mass: %2.1f\n'%(self.atomName, self.atomType, self.charge, self.atomDescription, self.atomicMass)
	
	return tmpstr
    
	
	
