def ConvertDihedralFromProperDihedralToDihedralTrig(k, multiplicity):

    zu = 0*k.unit
    fc0 = k
    fc1 = zu
    fc2 = zu
    fc3 = zu
    fc4 = zu
    fc5 = zu
    fc6 = zu

    if (multiplicity == 1):
        fc1 = k
    elif (multiplicity == 2):
        fc2 = k
    elif (multiplicity == 3):
        fc3 = k
    elif (multiplicity == 4):
        fc4 = k
    elif (multiplicity == 5):
        fc5 = k
    elif (multiplicity == 6):
        fc6 = k
    return fc0, fc1, fc2, fc3, fc4, fc5, fc6

def ConvertDihedralFromFourierToDihedralTrig(F1, F2, F3, F4):

    zu = 0*F1.unit
    fc0 = 0.5*(F1+F2+F3+F4)
    fc1 = 0.5*F1
    fc2 = -0.5*F2
    fc3 = 0.5*F3
    fc4 = -0.5*F4
    fc5 = zu
    fc6 = zu

    return fc0, fc1, fc2, fc3, fc4, fc5, fc6

def ConvertDihedralFromDihedralTrigToFourier(fc0, fc1, fc2, fc3, fc4, fc5, fc6):

    F1 = 2*fc1
    F2 = -2*fc2
    F3 = 2*fc3
    F4 = -2*fc4

    if fc0 != 0.5*(F1+F2+F3+F4):
        print "This dihedral is inconsistent with OPLS format",

    return F1, F2, F3, F4
    
def ConvertDihedralFromRBToOPLS(c0,c1,c2,c3,c4,c5,c6):

    if (c5 !=0.0 * c0.unit and c1+c2+c3+c4 != 0.0 * c0.unit):
        print "This RB dihedral is inconsistent with OPLS style",
        print "because c5 = ",c5,
        print " (should be 0) and c1+c2+c3+c4 = ",c1+c2+c3+c4,
        print " (should be 0)"
        # REALLY SHOULD ADD SOME SORT OF EXCEPTION HERE.
    # note - f1 and f3 are opposite sign as expected in GROMACS, probably because of angle conventions.
    f1 = 2.0 * c1 + 3.0 * c3 / 2.0
    f2 = -c2 - c4
    f3 = c3 / 2.0
    f4 = -c4 / 4.0
    return f1, f2, f3, f4

def ConvertDihedralFromOPLSToRB(f1,f2,f3,f4):

    # Note: c1 and c3 are the negative of what is defined on equation 4.64 of Gromacs Manual 4.6.1
    c0 = f2 + 0.5*(f1+f3)
    c1 = 0.5 * (f1 - 3.0 * f3)
    c2 = -f2 + 4.0 * f4
    c3 = 2.0 * f3
    c4 = -4.0 * f4
    c5 = 0.0 * c0.unit  # need to keep everything in units
    c6 = 0.0 * c0.unit
    return c0, c1, c2, c3, c4, c5, c6

def ConvertDihedralFromDihedralTrigToRB(sign, phi, fc0, fc1, fc2, fc3, fc4, fc5, fc6):

    # sign is -1 or 1
    # RB is \sum_n=0^6 cos(x)^n
    # DihedralTrig is f_0 + \sum_n=1^6 cos(nx-phi)
    # we restrict to phi = 0 or 180 (might need to generalize later), so we can write as
    #               f_0 + \sum_n=1^6 cos(nx)cos(phi) + sin(nx)sin(phi)
    #               f_0 + \sum_n=1^6 cos(nx)  (phi = 0)
    #               f_0 + \sum_n=1^6 -cos(nx)  (phi = 180)
    # phi corresponds to a sign of -1 on everything but f_0.  Easier to multiply f_0 by sign, and
    # then multiply by everything by sign at the end. 
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

    c0 = sign*fc0 - fc2 + fc4 - fc6
    c1 = fc1 - 3.0*fc3 + 5.0*fc5
    c2 = 2.0*fc2 - 8.0*fc4 + 18.0*fc6
    c3 = 4.0*fc3 - 20.0*fc5
    c4 = 8.0*fc4 - 48.0*fc6
    c5 = 16.0*fc5
    c6 = 32.0*fc6

    # Multiplying by -1 on odd powers to switch between sign conventions
    return sign*c0,-sign*c1,sign*c2,-sign*c3,sign*c4,-sign*c5,sign*c6,

def ConvertDihedralFromRBToDihedralTrig(c0, c1, c2, c3, c4, c5, c6):

    # see above for conversion; simply inverting the matrix.  Need to handle sign for 180?

    fc0 =  1.0*c0 + 0.5*c2 + 0.3750*c4 + 0.3125*c6
    fc1 =  1.0*c1 + 0.75*c3 + 0.6250*c5
    fc2 =  0.5*c2 + 0.5*c4 + 0.46875*c6
    fc3 =  0.25*c3 + 0.3125*c5
    fc4 =  0.125*c4 + 0.1875*c6 
    fc5 =  0.0625*c5
    fc6 =  0.03125*c6

    # Multiplying by -1 on odd powers to switch between sign conventions
    return fc0, -fc1, fc2, -fc3, fc4, -fc5, fc6

