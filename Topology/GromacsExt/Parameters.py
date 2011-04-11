from Decorators import *


# Default Parameters
class DefaultParameter(object):
    """
    [ defaults ]
    ; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
    1               2               yes             0.5     0.8333

    """

    @accepts_compatible_units(None, None, None, None, None, None)
    def __init__(self, func, cr, genpairs, fudgeLJ, fudgeQQ, comment=''):

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
class AtomType1Parameter(object):
    """
    [ atomtypes ]
    ; name     bond_type  mass   charge ptype  sigma      epsilon
    amber99_0     H0     0.0000  0.0000  A   2.47135e-01  6.56888e-02

    """

    @accepts_compatible_units(None, units.amu, units.elementary_charge, None, units.kilojoules_per_mole*units.nanometers**6, units.kilojoules_per_mole*units.nanometers**12, None,  None, None)
    def __init__(self, name, mass, charge, ptype, V, W, bondtype='', Z='', comment=''):

        self.name = name
        self.bondtype = bondtype
        self.Z = Z
        self.mass = mass
        self.charge = charge
        self.ptype = ptype
        self.V = V
        self.W = W
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s%7s%4s%11.4f%8.4f%3s%14e%13e; %s\n'%(self.name, self.bondtype, self.Z, self.mass._value, self.charge._value, self.ptype, self.V._value, self.W._value, self.comment)

class AtomType23Parameter(object):
    """
    [ atomtypes ]
    ; name     bond_type  mass   charge ptype  sigma      epsilon
    amber99_0     H0     0.0000  0.0000  A   2.47135e-01  6.56888e-02

    """

    @accepts_compatible_units(None, units.amu, units.elementary_charge, None, units.nanometers, units.kilojoules_per_mole, None,  None, None)
    def __init__(self, name, mass, charge, ptype, V, W, bondtype='', Z='',comment=''):

        self.name = name
        self.bondtype = bondtype
        self.Z = Z
        self.mass = mass
        self.charge = charge
        self.ptype = ptype
        self.V = V
        self.W = W
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive

        """

        return '%s%7s%4s%11.4f%8.4f%3s%14e%13e; %s\n'%(self.name, self.bondtype, self.Z, self.mass._value, self.charge._value, self.ptype, self.V._value, self.W._value, self.comment)


# Bond parameters
class Bond(object):
    """
    [ bondtypes ]
    ; i    j  func       b0          kb
    CT   H0       1    0.10900   284512.0 ; 03GLY changed from 331 bsd on NMA nmodes; AA, SUGARS

    """

    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole*units.nanometers**(-2), None)
    def __init__(self, atomtype1, atomtype2, func, distance, kspring, comment=''):

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

        return '%s   %s       %d    %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.distance._value, self.kspring._value, self.comment)


class G96Bond(object):
    """
    [ bondtypes ]
    ; i    j  func       b0          kb
    CT   H0       2    0.10900   284512.0 ;

    """

    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole*units.nanometers**(-4), None)
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

        return '%s   %s       %d    %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.distance._value, self.kspring._value, self.comment)


class MorseBond(object):
    """
    [ bondtypes ]
    ; i    j  func       b0          D         beta
    CT   H0       3    0.10900   284512.0   0.10900 ;

    """

    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole, units.nanometers**(-1), None)
    def __init__(self, atomtype1, stomtype2, func, distance, D, beta, comment=''):

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

        return '%s   %s       %d    %f   %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.distance._value, self.D._value, self.beta._value, self.comment)


class CubicBond(object):
    """
    [ bondtypes ]
    ; i    j  func       b0         C2         C3
    CT   H0       4    0.10900   284512.0   284512.0 ;

    """

    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole*units.nanometers**(-2), units.kilojoules_per_mole*units.nanometers**(-3), None)
    def __init__(self, atomtype1, atomtype2, func, distance, C2, C3, comment=''):

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

        return '%s   %s       %d    %f   %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.distance._value, self.C2._value, self.C3._value, self.comment)


class ConnectionBond(object):
    """
    [ bondtypes ]
    ; i    j  func
    CT   H0       5 ;

    """

    @accepts_compatible_units(None, None, None, None)
    def __init__(self, atomtype1, atomtype2, func, comment=''):

        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype2
        self.func = func
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding to a line as in a GromacsTopologyFile directive

        """

        return '%s   %s       %d ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.comment)


class HarmonicBond(object):
    """
    [ bondtypes ]
    ; i    j  func       b0          kb
    CT   H0       6    0.10900   284512.0 ;

    """

    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole*units.nanometers**(-2), None)
    def __init__(self, atomtype1, atomtype2, func, distance, kspring, comment=''):

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

        return '%s   %s       %d    %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.distance._value, self.kspring._value, self.comment)


# Pair parameters
class PairType11(object):
    """
    [ pairtypes ]
    ; i    j  func       sigma         epsilon
    CT   H0       1    3.56359e-01  1.04600e+00 ;

    """

    @accepts_compatible_units(None, None, None, units.kilojoules_per_mole*units.nanometers**6, units.kilojoules_per_mole*units.nanometers**12, None)
    def __init__(self, atomtype1, atomtype2, func, V, W, comment=''):

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

        return '%s   %s       %d    %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.V, self.W, self.comment)


class PairType123(object):
    """
    [ pairtypes ]
    ; i    j  func       sigma         epsilon
    CT   H0       1    3.56359e-01  1.04600e+00 ;

    """

    @accepts_compatible_units(None, None, None, units.nanometers**6, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, func, V, W, comment=''):

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

        return '%s   %s       %d    %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.V, self.W, self.comment)


class PairType21(object):
    """
    [ pairtypes ]
    ; i    j  func       fudgeQQ     qi       qj       sigma         epsilon
    CT   H0       2      0.8333    0.0000   0.0000   3.56359e-01   1.04600e+00 ;

    """

    @accepts_compatible_units(None, None, None, None, units.elementary_charge, units.elementary_charge, units.nanometers**6, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, func, fudgeQQ, qi, qj, V, W, comment=''):

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

        return '%s   %s       %d    %f   %f   %f   %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.fudgeQQ, self.qi._value, self.qj._value, self.V._value, self.W._value, self.comment)


class PairType223(object):
    """
    [ pairtypes ]
    ; i    j  func       fudgeQQ     qi       qj       sigma         epsilon
    CT   H0       2      0.8333    0.0000   0.0000   3.56359e-01   1.04600e+00 ;

    """

    @accepts_compatible_units(None, None, None, None, units.elementary_charge, units.elementary_charge, units.kilojoules_per_mole*units.nanometers**6, units.kilojoules_per_mole*units.nanometers**12, None)
    def __init__(self, atomtype1, atomtype2, func, fudgeQQ, qi, qj, V, W, comment=''):

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

        return '%s   %s       %d    %f    %f   %f   %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.fudgeQQ, self.qi._value, self.qj._value, self.V._value, self.W._value, self.comment)


class PairTypeNB13(object):
    """
    [ pairtypes_nb ]
    ; i    j  func         qi       qj       sigma         epsilon
    CT   H0       1      0.0000   0.0000   3.56359e-01   1.04600e+00 ;

    """

    @accepts_compatible_units(None, None, None, units.elementary_charge, units.elementary_charge, units.kilojoules_per_mole*units.nanometers**6, units.kilojoules_per_mole*units.nanometers**12, None)
    def __init__(self, atomtype1, atomtype2, func, qi, qj, V, W, comment=''):

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

        return '%s   %s       %d    %f   %f   %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.qi._value, self.qj._value, self.V._value, self.W._value, self.comment)


class PairTypesNB2(object):
    """
    [ pairtypes_nb ]
    ; i    j  func         qi       qj       sigma         epsilon
    CT   H0       1      0.0000   0.0000   3.56359e-01   1.04600e+00 ;

    """

    @accepts_compatible_units(None, None, None, units.elementary_charge, units.elementary_charge, units.nanometers, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, func, qi, qj, V, W, comment=''):

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

        return '%s   %s       %d    %f   %f   %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.qi._value, self.qj._value, self.V._value, self.W._value, self.comment)


# Angle parameters
class Angle(object):
    """
    [ angletypes ]
    ;  i    j    k  func       th0       cth
    H0  CT  H0           1   109.500    292.880 ; 03GLY

    """

    @accepts_compatible_units(None, None, None, None, units.degrees, units.kilojoules_per_mole * units.radians**(-2), None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, theta, k, comment=''):

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

        return '%s   %s   %s       %d    %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.theta._value, self.k._value, self.comment)


class G96Angle(object):
    """
    [ angletypes ]
    ;  i    j    k  func       th0       cth
    H0  CT  H0           2   109.500    292.880 ;

    """

    @accepts_compatible_units(None, None, None, None, units.degrees, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, theta, k, comment=''):

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

        return '%s   %s   %s       %d    %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.theta._value, self.k._value, self.comment)


class CrossBondBondAngle(object):
    """
    [ angletypes ]
    ;  i    j    k  func       r1e       r2e         cth
    H0  CT  H0           3   0.10900    0.10900    292.880 ;

    """

    @accepts_compatible_units(None, None, None, None, units.nanometers, units.nanometers, units.kilojoules_per_mole*units.nanometers**(-2), None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, r1e, r2e, krr, comment=''):

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

        return '%s   %s   %s       %d    %f   %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.r1e._value, self.r2e._value, self.krr._value, self.comment)


class CrossBondAngleAngle(object):
    """
    [ angletypes ]
    ;  i    j    k  func       r1e       r2e         r3e        cth
    H0  CT  H0           4   0.10900    0.10900    0.10900    292.880 ;

    """

    @accepts_compatible_units(None, None, None, None, units.nanometers, units.nanometers, units.nanometers, units.kilojoules_per_mole*units.nanometers**(-2), None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, r1e, r2e, r3e, kr, comment=''):

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

        return '%s   %s   %s       %d    %f   %f   %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.r1e._value, self.r2e._value, self.r3e._value, self.kr._value, self.comment)


class UreyBradleyAngle(object):
    """
    [ angletypes ]
    ;  i    j    k  func       th0       cth         r13        kUB
    H0  CT  H0           5   109.500    292.880    0.10900    292.880 ;

    """

    @accepts_compatible_units(None, None, None, None, units.degrees, units.kilojoules_per_mole, units.nanometers, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, theta, k, r13, kUB, comment=''):

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

        return '%s   %s   %s       %d    %f   %f   %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.theta._value, self.k._value, self.r13._value, self.kUB._value, self.comment)


class QuarticAngle(object):
    """
    [ angletypes ]
    ;  i    j    k  func       th0        c0         c1         c2         c3         c4
    H0  CT  H0           6   109.500    292.880    292.880    292.880    292.880    292.880 ;

    """

    @accepts_compatible_units(None, None, None, None, units.degrees, units.kilojoules_per_mole, units.kilojoules_per_mole*units.radians**(-1), units.kilojoules_per_mole*units.radians**(-2), units.kilojoules_per_mole*units.radians**(-3), units.kilojoules_per_mole*units.radians**(-4), None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, theta, C0, C1, C2, C3, C4, comment=''):

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

        return '%s   %s   %s       %d    %f   %f   %f   %f   %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.theta._value, self.C0._value, self.C1._value, self.C2._value, self.C3._value, self.C4._value, self.comment)


# LJ/Coul nonbonded parameters
class LJNonbonded13(object):
    """
    [ nonbond_params ]
    ; i    j  func       sigma         epsilon
    CT   H0       1    3.56359e-01  1.04600e+00 ;

    """

    @accepts_compatible_units(None, None, None, units.kilojoules_per_mole*units.nanometers**6, units.kilojoules_per_mole*units.nanometers**12, None)
    def __init__(self, atomtype1, atomtype2, func, V, W, comment=''):

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

        return '%s   %s       %d    %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.V, self.W, self.comment)


class LJNonbonded2(object):
    """
    [ nonbond_params ]
    ; i    j  func       sigma         epsilon
    CT   H0       1    3.56359e-01  1.04600e+00 ;

    """

    @accepts_compatible_units(None, None, None, units.nanometers, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, func, V, W, comment=''):

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

        return '%s   %s       %d    %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.V, self.W, self.comment)


class BuckinghamNonbonded(object):
    """
    [ nonbond_params ]
    ; ai    aj             funct   a        b        c
    C       C       2 2.07861e+04   2.94117e+01   1.17244e-03

    """

    @accepts_compatible_units(None, None, None, units.kilojoules_per_mole, units.nanometers**(-1), units.kilojoules_per_mole * units.nanometers**6, None)
    def __init__(self, atomtype1, atomtype2, func, a, b, c6, comment=''):

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

        return '%s   %s       %d    %f   %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.a._value, self.b._value, self.c6._value, self.comment)


# Dihedral parameters
class ProperDihedral(object):
    """
    [ dihedraltypes ]
    ;  ai  aj      ak      al   funct     th0         kth0       mult
    C        C        C        C        1 35.264     334.720       2 ;

    """

    @accepts_compatible_units(None, None, None, None, None, units.degrees, units.kilojoules_per_mole, None, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, phi, kphi, multiplicity, comment=''):

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

        return '%s   %s   %s   %s       %d    %f   %f   %d ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.phi._value, self.kphi._value, self.multiplicity, self.comment)


class ImproperDihedral(object):
    """
    [ dihedraltypes ]
    ;  ai  aj      ak      al   funct     xi         kxi
    C        C        C        C        2 35.264     334.720        ; default params

    """

    @accepts_compatible_units(None, None, None, None, None, units.degrees, units.kilojoules_per_mole*units.radians**(-2), None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, xi, kxi, comment=''):

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

        return '%s   %s   %s   %s       %d    %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.xi._value, self.kxi._value, self.comment)


class RBDihedral(object):
    """
    [ dihedraltypes ]
    ;   ai    aj    ak    al funct       c0          c1          c2          c3          c4          c5
    C     C     C     C   3   -4.524996    2.558516    3.932960   -1.966480    0.000000    0.000000

    """

    @accepts_compatible_units(None, None, None, None, None, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, C0, C1, C2, C3, C4, C5, comment=''):

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

        return '%s   %s   %s   %s       %d    %f   %f   %f   %f   %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.C0._value, self.C1._value, self.C2._value, self.C3._value, self.C4._value, self.C5._value, self.comment)


class FourierDihedral(object):
    """
    [ dihedraltypes ]
    ;   ai    aj    ak    al funct       c0          c1          c2          c3          c4
    C     C     C     C   4   -4.524996    2.558516    3.932960   -1.966480    0.000000

    """

    @accepts_compatible_units(None, None, None, None, None, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, units.kilojoules_per_mole, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, C1, C2, C3, C4, comment=''):

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

        return '%s   %s   %s   %s       %d    %f   %f   %f   %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.C1._value, self.C2._value, self.C3._value, self.C4._value, self.comment)


# Constraint parameters
class Constraint(object):
    """
    [ constrainttypes ]
    ;  ai    aj funct         dist
    C    C    1            0.16432

    """

    @accepts_compatible_units(None, None, None, units.nanometers, None)
    def __init__(self, atomtype1, atomtype2, func, distance, comment=''):

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

        return '%s   %s       %d    %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.distance._value, self.comment)


class ConstraintNC(object):
    """
    [ constrainttypes ]
    ;  ai    aj funct         dist
    C    C    2           0.16432

    """

    @accepts_compatible_units(None, None, None, units.nanometers, None)
    def __init__(self, atomtype1, atomtype2, func, distance, comment=''):

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

        return '%s   %s       %d    %f ; %s\n'%(self.atomtype1, self.atomtype2, self.func, self.distance._value, self.comment)


# Settle parameters
class Settle(object):
    """
    [ settles ]
    ; OW     funct     doh     dhh
    1       1       0.09572 0.15139

    """

    @accepts_compatible_units(None, None, None, None, units.nanometers, units.nanometers, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, distanceOH, distanceHH, comment=''):

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

        return '%s   %s   %s       %d    %f    %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.distanceOH._value, self.distanceHH._value, self.comment)


# Dummy parameters
class Dummy2(object):
    """
    [ dummies2 ]
    ; Dummy from                        funct        a
    4        1        2        3        1        0.128012065

    """

    @accepts_compatible_units(None, None, None, None, None, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, func, a, comment=''):

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

        return '%s   %s   %s       %d    %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.func, self.a, self.comment)


class Dummy3(object):
    """
    [ dummies3 ]
    ; Dummy from                        funct        a                b
    4        1        2        3        1        0.128012065        0.128012065

    """

    @accepts_compatible_units(None, None, None, None, None, None, None, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, a, b, comment=''):

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

        return '%s   %s   %s   %s       %d    %f    %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.a, self.b, self.comment)


class Dummy3fd(object):
    """
    [ dummies3fd ]
    ; Dummy from                        funct        a                d
    4        1        2        3        2        0.128012065        0.128012065

    """

    @accepts_compatible_units(None, None, None, None, None, None, units.nanometers, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, a, d, comment=''):

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

        return '%s   %s   %s   %s       %d    %f    %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.a, self.d._value, self.comment)


class Dummy3fad(object):
    """
    [ dummies3fad ]
    ; Dummy from                        funct        th0              d
    4        1        2        3        3        0.128012065        0.128012065

    """

    @accepts_compatible_units(None, None, None, None, None, units.degrees, units.nanometers, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, theta, d, comment=''):

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

        return '%s   %s   %s   %s       %d    %f    %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.theta._value, self.d._value, self.comment)


class Dummy3out(object):
    """
    [ dummies3out ]
    ; Dummy from                        funct        a                b                c
    4        1        2        3        4        0.128012065        0.128012065    0.128012065

    """

    @accepts_compatible_units(None, None, None, None, None, None, None, units.nanometers**(-1), None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, a, b, c, comment=''):

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

        return '%s   %s   %s   %s       %d    %f    %f    %f ; %s\n'%(self.atomtype1, self.atomtype2, self.atomtype3, self.atomtype4, self.func, self.a, self.b, self.c._value, self.comment)


# Restraint parameters
class PositionRestraint(object):
    """
    [ position_restraints ]


    """

    @accepts_compatible_units(None, None, units.kilojoules_per_mole*units.nanometers**(-2), units.kilojoules_per_mole*units.nanometers**(-2), units.kilojoules_per_mole*units.nanometers**(-2), None)
    def __init__(self, atomtype1, func, kx, ky, kz, comment=''):

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

        return '%s       %d    %f    %f    %f ; %s\n'%(self.atomtype1, self.func, self.kx._value, self.ky._value, self.kz._value, self.comment)


class DistanceRestraint(object):
    """
    [ distance_restraints ]


    """

    @accepts_compatible_units(None, None, None, None, None, units.nanometers, units.nanometers, units.nanometers, None, None)
    def __init__(self, atomtype1, atomtype2, func, Type, label, low, up1, up2, weight, comment=''):

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


class OrientationRestraint(object):
    """
    [ orientation_restraints ]


    """

    @accepts_compatible_units(None, None, None, None, None, None, None, units.amu, units.amu**(-1), None)
    #c should be in units "U nm^alpha"
    #need some way to implement this if you want a unit check on this value
    def __init__(self, atomtype1, atomtype2, func, exp, label, alpha, c, obs, weight, comment=''):

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


class AngleRestraint(object):
    """
    [ angle_restraints ]


    """

    @accepts_compatible_units(None, None, None, None, None, units.degrees, units.kilojoules_per_mole, None, None)
    def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, func, theta, kc, multiplicity, comment=''):

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


class AngleZRestraint(object):
    """
    [ angle_restraints_z ]


    """

    @accepts_compatible_units(None, None, None, units.degrees, units.kilojoules_per_mole, None, None)
    def __init__(self, atomtype1, atomtype2,  func, theta, kc, multiplicity, comment=''):

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
class Exclusion(object):
    """
    [ exclusions ]
    1        2        3  ;

    """

    def __init__(self, atomtype1, indices, comment=''):

        self.atomtype1 = atomtype1
        self.indices = indices
        self.comment = comment
        return

    def directiveString(self):
        """Output a correctly-formatted string (with an ending newline) corresponding
        to a line as in a GromacsTopologyFile directive
        """

        return '%s %d ; %s\n'%(self.atomtype1, self.indices, self.comment)
