def ConvertDihedralFromRBToOPLS(c0,c1,c2,c3,c4,c5,c6):

    if (c5 !=0.0 and c1+c2+c3+c4 != 0.0):
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
    c5 = 0.0*c0.unit  # need to keep everything in units
    c6 = 0.0*c0.unit
    return c0, c1, c2, c3, c4, c5, c6

def ConvertDihedralFromProperTrigToRB(sign,f0,f1,f2,f3,f4,f5,f6):

    # sign is -1 or 1
    # RB is \sum_n=0^6 cos(x)^n
    # ProperTrig is f_0 + \sum_n=1^6 cos(nx-phi)
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

    c0 = sign*f0 - f2 + f4 - f6
    c1 = f1 - 3.0*f3 + 5.0*f5
    c2 = 2.0*f2 - 8.0*f4 + 18.0*f6
    c3 = 4.0*f3 - 20.0*f5
    c4 = 8.0*f4 - 48.0*f6
    c5 = 16.0*f5
    c6 = 32.0*f6
    return sign*c0, sign*c1, sign*c2, sign*c3, sign*c4, sign*c5, sign*c6

def ConvertDihedralFromRBToProperTrig(c0,c1,c2,c3,c4,c5,c6):

    # see above for conversion; simply inverting the matrix.
    f0 =  1.0*c0 + 0.5*c2 + 0.3750*c4 + 0.3125*c6
    f1 =  1.0*c1 + 0.75*c3 + 0.6250*c5
    f2 =  0.5*c2 + 0.5*c4 + 0.46875*c6
    f3 =  0.25*c3 + 0.3125*c5
    f4 =  0.125*c4 + 0.1875*c6 
    f5 =  0.0625*c5
    f6 =  0.03125*c6
    return f0, f1, f2, f3, f4, f5, f6

