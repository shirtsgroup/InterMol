from System import *


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
        """Output a correctly-formatted string (with an ending newline) corresponding to a line as in a GromacsTopologyFile.Directive object

        """

        pass


# Default Parameters
class GromacsDefaultParameterInfo(GromacsParameterInfo):


    def __init__(self, func, cr, genpairs, fudgeLJ, fudgeQQ, comment=''):
        """
        """
        self.func = func
        self.cr = cr
        self.genpairs = genpairs
        self.fudgeLJ = fudgeLJ
        self.fudgeQQ = fudgeQQ
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding to a line as in a GromacsTopologyFile directive

        """

        return '%d%16d%18s%16.1f%11.4f ; %s\n'%(self.func, self.cr, self.genpairs, self.fudgeLJ, self.fudgeQQ, self.comment)


# Atom Parameters
class GromacsAtomParameterInfo(GromacsParameterInfo):
    """
    [ atomtypes ]
    ; name     bond_type  mass   charge ptype  sigma      epsilon
    amber99_0     H0     0.0000  0.0000  A   2.47135e-01  6.56888e-02

    """

    # If combination rule = 1
    #@accepts_compatible_units(None, None, units.mass(amu?), units.elementary_charge, None, units.kilojoules_per_mole*units.nanometers**6, units.kilojoules_per_mole*units.nanometers**12, None)

    # If combination rule = 2 or 3
    #@accepts_compatible_units(None, None, units.mass(amu?), units.elementary_charge, None, units.nanometers, units.kilojoules_per_mole, None)
    def __init__(self, name, bondtype, mass, charge, ptype, sigma, epsilon, comment=''):
        """
        """
        self.name = name
        self.bontype = bondtype
        self.mass = mass
        self.charge = charge
        self.ptype = ptype
        self.sigma = sigma
        self.epsilon = epsilon
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s%7s%11.4f%8.4f%3s%14e%13e; %s\n'%(self.name, self.bondtype, self.mass._value, self.charge._value, self.ptype, self.sigma._value, self.epsilon._value, self.comment)


# Bond parameters
class GromacsBondParameterInfo(GromacsParameterInfo):
    """
    [ bondtypes ]
    ; i    j  func       b0          kb
    CT H0         1    0.10900   284512.0 ; 03GLY changed from 331 bsd on NMA nmodes; AA, SUGARS

    """

    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole / units.nanometers**2, None)
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
        """Output a correctly-formatted string (with an ending newline) corresponding to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d    %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.distance._value, self.kspring._value, self.comment)


class GromacsG96BondParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole / units.nanometers**4, None)
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
        """Output a correctly-formatted string (with an ending newline) corresponding to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d    %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.distance._value, self.kspring._value, self.comment)


class GromacsMorseBondParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole, units.nanometers**(-1), None)
    def __init__(self, atomtype1, stomtype2, func, distance, D, beta, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.func = func
        self.distance = distance
        self.D = D
        self.beta = beta
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d    %f %f %f; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.distance._value, self.D._value, self.beta._value, self.comment)



class GromacsCubicBondParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole / units.nanometers**2, units.kilojoules_per_mole / units.nanometers**3, None)
    def __init__(self, atomtype1, atomtype2, func, distance, C2, C3, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.func = func
        self.distance = distance
        self.C2 = C2
        self.C3 = C3
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d    %f %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.distance._value, self.C2._value, self.C3._value, self.comment)


class GromacsConnectionBondParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None)
    def __init__(self, atomtype1, atomtype2, func, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.func = func
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.comment)



class GromacsHarmonicBondParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole / units.nanometers**2, None)
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
        """Output a correctly-formatted string (with an ending newline) corresponding to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d    %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.distance._value, self.kspring._value, self.comment)

# Pair parameters
class GromacsPairTypes1ParameterInfo(GromacsParameterInfo):

    # If combination rule = 1
    #@accepts_compatible_units(None, None, None, units.kilojoules_per_mole*units.nanometers**6, units.kilojoules_per_mole*units.nanometers**12, None)
    # If combination rule = 2 and 3
    #@accepts_compatible_units(None, None, None, units.nanometers**6, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, func, V, W, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.func = func
        self.V = V
        self.W = W
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d    %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.V, self.W, self.comment)

class GromacsPairTypes2ParameterInfo(GromacsParameterInfo):

    # If combination rule = 1
    #@accepts_compatible_units(None, None, None, None, units.elementary_charge, units.elementary_charge, units.nanometers**6, units.kilojoules_per_mole, None)

    # If combination rule = 2 or 3
    #@accepts_compatible_units(None, None, None, None, units.elementary_charge, units.elementary_charge, units.kilojoules_per_mole*units.nanometers**6, units.kilojoules_per_mole*units.nanometers**12, None)

    def __init__(self, atomtype1, atomtype2, func, fudgeQQ, qi, qj, V, W, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.func = func
        self.fudgeQQ = fudgeQQ
        self.qi = qi
        self.qj = qj
        self.V = V
        self.W = W
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d    %s %f %f %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.fudgeQQ, self.qi._value, self.qj._value, self.V._value, self.W._value, self.comment)

class GromacsPairTypesNBParameterInfo(GromacsParameterInfo):

    # If combination rule = 1 or 3
    #@accepts_compatible_units(None, None, None, units.elementary_charge, units.elementary_charge, units.kilojoules_per_mole*units.nanometers**6, units.kilojoules_per_mole*units.nanometers**12, None)

    # If combination rule = 2
    #@accepts_compatible_units(None, None, None, units.elementary_charge, units.elementary_charge, units.nanometers, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, func, qi, qj, V, W, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.func = func
        self.qi = qi
        self.qj = qj
        self.V = V
        self.W = W
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d    %f %f %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.qi._value, self.qj._value, self.V._value, self.W._value, self.comment)

# Angle parameters
class GromacsAngleParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, units.degrees, units.kilojoules_per_mole * units.radians**(-2), None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, theta, k, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.func = func
        self.theta = theta
        self.k = k
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %d    %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.theta._value, self.k._value, self.comment)

class GromacsG96AngleParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, units.degrees, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, theta, k, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.func = func
        self.theta = theta
        self.k = k
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %d    %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.theta._value, self.k._value, self.comment)


class GromacsCrossBondBondAngleParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, units.nanometers, units.nanometers, units.kilojoules_per_mole*units.nanometers**(-2), None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, r1e, r2e, krr, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.func = func
        self.r1e = r1e
        self.r2e = r2e
        self.krr = krr
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %d    %f %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.r1e._value, self.r2e._value, self.krr._value, self.comment)


class GromacsCrossBondAngleAngleParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, units.nanometers, units.nanometers, units.nanometers, units.kilojoules_per_mole*units.nanometers**(-2), None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, r1e, r2e, r3e, kr, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.func = func
        self.r1e = r1e
        self.r2e = r2e
        self.r3e = r3e
        self.kr = kr
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %d    %f %f %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.r1e._value, self.r2e._value, self.r3e._value, self.kr._value, self.comment)

class GromacsUreyBradleyAngleParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, units.degrees, units.kilojoules_per_mole, units.nanometers, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, theta, k, r13, kUB, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.func = func
        self.theta = theta
        self.k = k
        self.r13 = r13
        self.kUB = kUB
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %d    %f %f %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.theta._value, self.k._value, self.r13._value, self.kUB._value, self.comment)

class GromacsQuarticAngleParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, units.degrees, units.kilojoules_per_mole, units.kilojoules_per_mole*units.radians**(-1), units.kilojoules_per_mole*units.radians**(-2), units.kilojoules_per_mole*units.radians**(-3), units.kilojoules_per_mole*units.radians**(-4), None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, theta, C0, C1, C2, C3, C4, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.func = func
        self.theta = theta
        self.C0 = C0
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %d    %f %f %f %f %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.theta._value, self.C0._value, self.C1._value, self.C2._value, self.C3._value, self.C4._value, self.comment)


# LJ/Coul nonbonded parameters
class GromacsLJNonbondedParameterInfo(GromacsParameterInfo):

    # If combination rule = 1 or 3
    #@accepts_compatible_units(None, None, None, units.kilojoules_per_mole*units.nanometers**6, units.kilojoules_per_mole*units.nanometers**12, None)

    # If combination rule = 2
    #@accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, func, V, W, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.func = func
        self.V = V
        self.W = W
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d    %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.V, self.W, self.comment)


class GromacsBuckinghamNonbondedParameterInfo(GromacsParameterInfo):

    """
    [ nonbond_params ]
    ; ai    aj             funct   a        b        c
    C       C       2 2.07861e+04   2.94117e+01   1.17244e-03
    """

    @accepts_compatible_units(None, None, None, units.kilojoules_per_mole, units.nanometers**(-1), units.kilojoules_per_mole * units.nanometers**6, None)
    def __init__(self, atomtype1, atomtype2, func, a, b, c6, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.func = func
        self.a = a
        self.b = b
        self.c6 = c6
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d    %f %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.a._value, self.b._value, self.c6._value, self.comment)



# Dihedral parameters
class GromacsProperDihedralParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, None, units.degrees, units.kilojoules_per_mole, None, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, phi, kphi, multiplicity, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.atomtype4 = atomtype4
        self.func = func
        self.phi = phi
        self.kphi = kphi
        self.multiplicity = multiplicity
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %s  %d    %f %f %d ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.phi._value, self.kphi._value, self.multiplicity, self.comment)


class GromacsImproperDihedralParameterInfo(GromacsParameterInfo):
    """
    [ dihedrals ]
    ;  ai  aj      ak      al   funct     xi         kxi
    1        2        4        3        2 35.264     334.720        ; default params
    """

    @accepts_compatible_units(None, None, None, None, None, units.degrees, units.kilojoules_per_mole / units.radians**2, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, xi, kxi, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.atomtype4 = atomtype4
        self.func = func
        self.xi = xi
        self.kxi = kxi
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %s  %d    %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.xi._value, self.kxi._value, self.comment)



class GromacsRBDihedralParameterInfo(GromacsParameterInfo):
    """
    [ dihedrals ]
    ;   ai    aj    ak    al funct       c0          c1          c2          c3          c4          c5
    1     5     8     9   3   -4.524996    2.558516    3.932960   -1.966480    0.000000    0.000000
    """

    @accepts_compatible_units(None, None, None, None, None, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, C0, C1, C2, C3, C4, C5, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.atomtype4 = atomtype4
        self.func = func
        self.C0 = C0
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4
        self.C5 = C5
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %s  %d    %f %f %f %f %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.C0._value, self.C1._value, self.C2._value, self.C3._value, self.C4._value, self.C5._value, self.comment)


class GromacsFourierDihedralParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, None, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, C1, C2, C3, C4, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.atomtype4 = atomtype4
        self.func = func
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %s  %d    %f %f %f %f %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.C1._value, self.C2._value, self.C3._value, self.C4._value, self.comment)


# Constraint parameters
class GromacsConstraintParameterInfo(GromacsParameterInfo):
    """
    [ constraints ]
    ;  ai    aj funct         dist
    2    1    1            0.16432
    """

    @accepts_compatible_units(None, None, None, units.nanometers, None)
    def __init__(self, atomtype1, atomtype2, func, distance, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.func = func
        self.distance = distance
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d    %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.distance._value, self.comment)


class GromacsConstraintNCParameterInfo(GromacsParameterInfo):
    """
    [ constraints ]
    ;  ai    aj funct         dist
    2    1    2           0.16432
    """

    @accepts_compatible_units(None, None, None, units.nanometers, None)
    def __init__(self, atomtype1, atomtype2, func, distance, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.func = func
        self.distance = distance
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d    %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.distance._value, self.comment)



# Settle parameters
class GromacsSettleParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, units.nanometers, units.nanometers, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, distanceOH, distanceHH, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.func = func
        self.distanceOH = distanceOH
        self.distanceHH = distanceHH
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s  %d    %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.distanceOH._value, self.distanceHH._value, self.comment)



# Dummy parameters
class GromacsDummy2ParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, None, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, a, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.func = func
        self.a = a
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s  %d    %s ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.a, self.comment)



class GromacsDummy3ParameterInfo(GromacsParameterInfo):
    """
    [ dummies3 ]
    ; Dummy from                        funct        a                b
    4        1        2        3        1        0.128012065        0.128012065
    """

    @accepts_compatible_units(None, None, None, None, None, None, None, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, a, b, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.atomtype4 = atomtype4
        self.func = func
        self.a = a
        self.b = b
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %s  %d    %s %s ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.a, self.b, self.comment)


class GromacsDummy3fdParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, None, None, units.nanometers, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, a, d, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.atomtype4 = atomtype4
        self.func = func
        self.a = a
        self.d = d
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %s  %d    %s %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.a, self.d._value, self.comment)


class GromacsDummy3fadParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, None, units.degrees, units.nanometers, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, theta, d, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.atomtype4 = atomtype4
        self.func = func
        self.theta = theta
        self.d = d
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %s  %d    %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.theta._value, self.d._value, self.comment)


class GromacsDummy3outParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, None, None, None, units.nanometers**(-1), None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, a, b, c, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.atomtype4 = atomtype4
        self.func = func
        self.a = a
        self.b = b
        self.c = c
        self.comment = comment
        return


    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %s  %d    %s %s %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.a, self.b, self.c._value, self.comment)



# restraint parameters
class GromacsPositionRestraintParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, units.kilojoules_per_mole / units.nanometers**2, units.kilojoules_per_mole / units.nanometers**2, units.kilojoules_per_mole / units.nanometers**2, None)
    def __init__(self, atomtype1, func, kx, ky, kz, comment=''):

        """
        """
        self.atomtype1 = atomtype1
        self.func = func
        self.kx = kx
        self.ky = ky
        self.kz = kz
        self.comment = comment
        return


    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s  %d    %f %f %f ; %s\n'%(self.atomtype1, self.func, self.kx._value, self.ky._value, self.kz._value, self.comment)

class GromacsDistanceRestraintParameterInfo(GromacsParameterInfo):

    @accepts_compatible_units(None, None, None, None, None, units.nanometers, units.nanometers, units.nanometers, None, None)
    def __init__(self, atomtype1, atomtype2, func, Type, label, low, up1, up2, weight, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.func = func
        self.Type = Type
        self.label = label
        self.low = low
        self.up1 = up1
        self.up2 = up2
        self.weight = weight
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d    %s %s %f %f %f %s ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.Type, self.label, self.low._value, self.up1._value, self.up2._value, self.weight, self.comment)

class GromacsOrientationRestraintParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, None, None, None, units.amu, units.amu**(-1), None)
    #c should be in units "U nm^alpha"
    #need some way to implement this if you want a unit check on this value
    def __init__(self, atomtype1, atomtype2, func, exp, label, alpha, c, obs, weight, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.func = func
        self.exp = exp
        self.label = label
        self.alpha = alpha
        self.c = c
        self.obs = obs
        self.weight = weight
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s  %d    %s %s %s %f %f %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.exp, self.label, self.alpha, self.c, self.obs._value, self.weight._value, self.comment)


class GromacsAngleRestraintParameterInfo(GromacsParameterInfo):


    @accepts_compatible_units(None, None, None, None, None, units.degrees, units.kilojoules_per_mole, None, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, theta, kc, multiplicity, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.atomtype3 = atomtype3
        self.atomtype4 = atomtype4
        self.func = func
        self.theta = theta
        self.kc = kc
        self.multiplicity
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %s %s %d    %f %f %d ; %s\n'%(self.atomtype1, atomtype2, self.atomtype3, self.atomtype4, self.func, self.theta._value, self.kc._value, self.multiplicity, self.comment)


class GromacsAngleZRestraintParameterInfo(GromacsParameterInfo):

    @accepts_compatible_units(None, None, None, units.degrees, units.kilojoules_per_mole, None, None)
    def __init__(self, atomtype1, atomtype2,  func, theta, kc, multiplicity, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.func = func
        self.theta = theta
        self.kc = kc
        self.multiplicity
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s %s %d    %f %f %d ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.theta._value, self.kc._value, self.multiplicity, self.comment)


# Exclusion parameters
class GromacsExclusionParameterInfo(GromacsParameterInfo):
    """
    [ exclusions ]
    1        2        3  ;sadfsdfad
    """


    def __init__(self, atomtype1, indices, comment=''):
        """
        """
        self.atomtype1 = atomtype1
        self.indices = indices
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive
        """

        return '%s %d ; %s\n'%(self.atomtype1, self.indices, self.comment)
