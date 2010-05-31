#=============================================================================================
# Container Class for Gromacs Topology Parameter Info
#=============================================================================================

"""
GromacsParameter.py

    DESCRIPTION
    
    These are container objects to store forcefield parameters.  There should be one derived class for each of the 
    Parameter directives that are possible to list in a Gromacs topology file:
    
	defaults
	atomtypes
	bondtypes
	constrainttypes
	pairtypes
	angletypes
	dihedraltypes1
	dihedraltypes2
	dihedraltypes3
	nonbond params1
	nonbond params2
	
    REQUIREMENTS 	
    
	* Needs an installation of the simtk.units package.  See: https://simtk.org/home/python_units
    
	* The environment variable $GMXLIB needs to be defined, listing the directory where Gromacs forcefield and *.itp files
	  are stored.  This can be added to .bash_profile, and sourced, for example:
	  
	  export GMXLIB="/home/my-account/local/gromacs-3.1.4/share/top"
	
    NOTES 
    
    VAV (May 30, 2010):  I've written templates for these, and one full example (GromacsBondTypeParameterInfo),
    but these obviously need to expounded on.  Hopefully this is a job for MRS's students!  Each GromacsParameterInfo
    needs to have an __init__ method that is wrapped by a unit-enforcing decorator.
     
    The directiveString() method of these classes can output their contents into an appropriately-formatted Directive line
    for output to a Gromacs *.top file.
    
    There also need to be helper functions to do the opposite: i.w. use the information contained in the ParameterInfo
    derived classes to create a  series of Force objects with the correct parameters for given atomtypes.
    These helper functions are part of the GromacsTopology class. 

"""

class ParameterInfo(object):
    """Parent class for ParameterInfo objects."""


class GromacsParameterInfo(ParameterInfo):
    """Parent class for GromacsParameterInfo objects."""

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
	to a line as in a GromacsTopologyFile.Directive object """
        pass


class GromacsBondParameterInfo(GromacsParameterInfo):
    """
    [ bondtypes ]
    ; i    j  func       b0          kb
    CT H0         1    0.10900   284512.0 ; 03GLY changed from 331 bsd on NMA nmodes; AA, SUGARS

    """
    
    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole / units.nanometers**2, comment='')
    def __init__(self, atomtype1, atomtype2, func, distance, kspring, comment=''):
	"""
	"""
	self.atomtype1 = atomtype1
	self.atomtype2 = atomtype2
	self.func = func
	self.distance = distance
	self.kspring = kspring
        self.comment = comment
	return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
	to a line as in a GromacsTopologyFile directive  """
	return '%s %s  %d    %f %f ; %s\n'%(atomtype1, atomtype2, func, distance, kspring, comment)


# Bond parameters
class GromacsG96BondParameterInfo(GromacsParameterInfo):
    pass
class GromacsMorseBondParameterInfo(GromacsParameterInfo):
    pass
class GromacsCubicBondParameterInfo(GromacsParameterInfo):
    pass
class GromacsConnectionBondParameterInfo(GromacsParameterInfo):
    pass
class GromacsHarmonicBondParameterInfo(GromacsParameterInfo):
    pass

# LJ/Coul nonbonded parameters
class GromacsLJNonbondedParameterInfo(GromacsParameterInfo):
    pass
class GromacsBuckinghamNonbondedParameterInfo(GromacsParameterInfo):
    pass

# Dihedral parameters
class GromacsProperDihedralParameterInfo(GromacsParameterInfo):
    pass
class GromacsImproperDihedralParameterInfo(GromacsParameterInfo):
    pass
class GromacsRBDihedralParameterInfo(GromacsParameterInfo):
    pass

# constraint parameters
class GromacsConstraintParameterInfo(GromacsParameterInfo):
    pass
class GromacsConstraintNCParameterInfo(GromacsParameterInfo):
    pass

# settle parameters
class GromacsSettleParameterInfo(GromacsParameterInfo):
    pass

# dummy parameters
class GromacsDummy2ParameterInfo(GromacsParameterInfo):
    pass
class GromacsDummy3ParameterInfo(GromacsParameterInfo):
    pass
class GromacsDummy3fdParameterInfo(GromacsParameterInfo):
    pass
class GromacsDummy3fadParameterInfo(GromacsParameterInfo):
    pass
class GromacsDummy3outParameterInfo(GromacsParameterInfo):
    pass

# restraint parameters
class GromacsPositionRestraintParameterInfo(GromacsParameterInfo):
    pass
class GromacsDistanceRestraintParameterInfo(GromacsParameterInfo):
    pass
class GromacsOrientationRestraintParameterInfo(GromacsParameterInfo):
    pass
class GromacsAngleRestraintParameterInfo(GromacsParameterInfo):
    pass
class GromacsAngleZRestraintParameterInfo(GromacsParameterInfo):
    pass

# Exclusion parameters
class GromacsExclusionParameterInfo(GromacsParameterInfo):
    pass
