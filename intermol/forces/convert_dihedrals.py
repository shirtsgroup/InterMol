from __future__ import print_function

import math

import parmed.unit as units

from intermol.exceptions import (UnimplementedFunctional, UnsupportedFunctional,
                                 UnimplementedSetting, UnsupportedSetting,
                                 GromacsError, InterMolError)

def convert_nothing(x):
    """ useful utility for not converting anything"""
    return x


def convert_dihedral_from_proper_to_trig(p):

    k = p['k']
    multiplicity = p['multiplicity']
    zu = 0*k.unit
    fcs = {
        'phi': p['phi'],
        'fc0': k,
        'fc1': zu,
        'fc2': zu,
        'fc3': zu,
        'fc4': zu,
        'fc5': zu,
        'fc6': zu
        }

    k # which force constant is nonzero because of the multiplicity?
    fk = "fc%d" % (multiplicity._value)
    fcs[fk] = k
    return fcs


def convert_dihedral_from_fourier_to_trig(f):

    fcs = dict()
    F1 = f['c1']
    F2 = f['c2']
    F3 = f['c3']
    F4 = f['c4']

    zu = 0*F1.unit

    fcs['phi'] = 0 * units.degrees
    fcs['fc0'] = 0.5*(F1+F2+F3+F4)
    fcs['fc1'] = 0.5*F1
    fcs['fc2'] = -0.5*F2
    fcs['fc3'] = 0.5*F3
    fcs['fc4'] = -0.5*F4
    fcs['fc5'] = zu
    fcs['fc6'] = zu

    return fcs


def convert_dihedral_from_trig_to_fourier(fcs):

    F = dict()
    F['F1'] = 2*fcs['fc1']
    F['F2'] = -2*fcs['fc2']
    F['F3'] = 2*fcs['fc3']
    F['F4'] = -2*fcs['fc4']
    F['F5'] = 0  # probably not correct?  No test cases.

    if fcs['fc0'] != 0.5*(F['F1']+F['F2']+F['F3']+F['F4']):
        print("This dihedral is inconsistent with OPLS format")

    return F


def convert_dihedral_from_trig_to_proper(fcs, convention='0'):

    # this has to be smarter, because there are two options; there's one, or there are multiple.
    # we handle this by returning multiple keywords

    # create a copy of the dictionary to strip out several parameters
    # TODO: Confusing and probably unneccesary. We only want the non-zero coeffs.
    # ftmp = fcs.copy()
    # ftmp.pop('phi')
    # ftmp.pop('fc0')
    # coefficients = ftmp.values()
    # ncount = sum(coeff._value != 0.0 for coeff in coefficients)

    plist = []
    for parameter_key, parameter_value in fcs.items():  # only one of these should be nonzero
        if parameter_value._value != 0.0 and parameter_key not in ['fc0', 'phi']:
            p = dict()
            p['phi'] = fcs['phi']
            p['multiplicity'] = int(parameter_key[2]) * units.dimensionless
            p['weight'] = 0 * units.dimensionless
            p['k'] = parameter_value
            plist.append(p)
    # all parameters are zero, so the force constant is zero if fc0 is zero
    # The other possibility is to just return no parameters here, since such
    # a dihedral has no energy

    if len(plist) == 0:
        # all zeros; check if fc0 is 0. if so, the force_constant is zero.
        if fcs['fc0']._value == 0.0:
            p = dict()
            p['phi'] = fcs['phi']
            p['multiplicity'] = int(parameter_key[2]) * units.dimensionless
            p['weight'] = 0 * units.dimensionless
            p['k'] = fcs['fc0']
            plist.append(p)
        else:
            raise InterMolError("Unable to convert dihedral from trig to proper")
    return plist


def convert_dihedral_from_RB_to_OPLS(c):

    c0 = c['C0']
    c1 = c['C1']
    c2 = c['C2']
    c3 = c['C3']
    c4 = c['C4']
    c5 = c['C5']

    f = dict()
    if (c5 !=0.0 * c0.unit and c1+c2+c3+c4 != 0.0 * c0.unit):
        print("This Rb dihedral is inconsistent with OPLS style")
        print("because C5 = ", c5)
        print(" (should be 0) and c1+c2+c3+c4 = ", c1+c2+c3+c4)
        print(" (should be 0)")
        # REALLY SHOULD ADD SOME SORT OF EXCEPTION HERE.
    # note - f1 and f3 are opposite sign as expected in GROMACS, probably because of angle conventions.
    f['f1'] = 2.0 * c1 + 3.0 * c3 / 2.0
    f['f2'] = -c2 - c4
    f['f3'] = c3 / 2.0
    f['f4'] = -c4 / 4.0
    return f


def convert_dihedral_from_OPLS_to_RB(f):

    f1 = f['f1']
    f2 = f['f2']
    f3 = f['f3']
    f4 = f['f4']

    c = dict()
    # Note: c1 and c3 are the negative of what is defined on equation 4.64 of Gromacs Manual 4.6.1
    c['C0'] = f2 + 0.5*(f1+f3)
    c['C1'] = 0.5 * (f1 - 3.0 * f3)
    c['C2'] = -f2 + 4.0 * f4
    c['C3'] = 2.0 * f3
    c['C4'] = -4.0 * f4
    c['C5'] = 0.0 * c['CO'].unit  # need to keep everything in units
    c['C6'] = 0.0 * c['CO'].unit
    return c


def convert_dihedral_from_trig_to_RB(fcs):

    # sign is -1 or 1
    # RB is \sum_n=0^6 cos(x)^n
    # Trig is f_0 + \sum_n=1^6 cos(nx-phi)
    # we restrict to phi = 0 or 180 (might need to generalize later), so we can write as
    #               f_0 + \sum_n=1^6 cos(nx)cos(phi) + sin(nx)sin(phi)
    #               f_0 + \sum_n=1^6 cos(nx)  (phi = 0)
    #               f_0 + \sum_n=1^6 -cos(nx)  (phi = 180)
    # phi = 180 corresponds to a sign of -1 on everything but f_0.  Easier to multiply just f_0, 
    # then multiply the rest. 
    # cos(2x) = 2cos^2(x)-1
    # cos(3x) = 4cos^3(x)-3cos(x)
    # cos(4x) = 8cos^4(X)-8cos^2(x)+1
    # cos(5x) = 16cos^5(x)-20cos^3(x)+5cos(x)
    # cos(6x) = 32cos^6(x)-48cos^4(x)+18cos^2(x)-1
    # Thus:
    #   f0 + f1*cos(x) + f2*cos(2x) + f3*cos(3x) + f4*cos(4x) + f5*cos(5x) + f6*cos(6x) 
    #   z = cos(x)
    # = f0 + f1*z + f2*(2z^2-1) + f3*(4*z^3-3*z) + f4*(8z^4-8z^2+1) + f5*(16z^5-20z^3+5*z) + f6*(32z^6-48z^4+18*z^2-1) 
    # c0 =  f0-f2+f4-f6
    # c1 =  f1-3f3+5f5
    # c2 = 2f2-8f4+18f6
    # c3 = 4f3-20f5
    # c4 = 8f4-48f6
    # c5 = 16f5
    # c6 = 32f6

    # phi = 180 corresponds to a sign of -1 on everything but f_0
    fc0 = fcs['fc0']
    sign = math.cos(fcs['phi'].value_in_unit(units.radians))
    fc1 = sign * fcs['fc1']
    fc2 = sign * fcs['fc2']
    fc3 = sign * fcs['fc3']
    fc4 = sign * fcs['fc4']
    fc5 = sign * fcs['fc5']
    fc6 = sign * fcs['fc6']

    c = dict()
    c['C0'] = fc0 - fc2 + fc4 - fc6
    c['C1'] = fc1 - 3.0*fc3 + 5.0*fc5
    c['C2'] = 2.0*fc2 - 8.0*fc4 + 18.0*fc6
    c['C3'] = 4.0*fc3 - 20.0*fc5
    c['C4'] = 8.0*fc4 - 48.0*fc6
    c['C5'] = 16.0*fc5
    c['C6'] = 32.0*fc6


    return c


def convert_dihedral_from_RB_to_trig(c):
    c0 = c['C0']
    c1 = c['C1']
    c2 = c['C2']
    c3 = c['C3']
    c4 = c['C4']
    c5 = c['C5']
    if 'C6' in c:  # program might not define this one, need to check it exists.
        c6 = c['C6']
    else:
        c6 = 0*c0.unit

    # See above for conversion; simply inverting the matrix.
    # Need to handle sign for 180.
    fcs = dict()
    # sign?  units?  Is there a way to get out of this?
    fcs['phi'] = 0 * units.degrees
    fcs['fc0'] = 1.0*c0 + 0.5*c2 + 0.3750*c4 + 0.3125*c6
    fcs['fc1'] = 1.0*c1 + 0.75*c3 + 0.6250*c5
    fcs['fc2'] = 0.5*c2 + 0.5*c4 + 0.46875*c6
    fcs['fc3'] = 0.25*c3 + 0.3125*c5
    fcs['fc4'] = 0.125*c4 + 0.1875*c6
    fcs['fc5'] = 0.0625*c5
    fcs['fc6'] = 0.03125*c6

    return fcs

