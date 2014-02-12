def ConvertFromRBToOPLSDihedral(c0,c1,c2,c3,c4,c5,c6):

    if (c5 !=0.0 and c1+c2+c3+c4 != 0.0):
        print "This RB dihedral is inconsistent with OPLS style",
        print "because c5 = ",c5,
        print " (should be 0) and c1+c2+c3+c4 = ",c1+c2+c3+c4,
        print " (should be 0)"
        # REALLY SHOULD ADD SOME SORT OF EXCEPTION HERE.
    f1 = -2.0 * c1 - 3.0 * c3 / 2.0
    f2 = -c2 - c4
    f3 = -c3 / 2.0
    f4 = -c4 / 4.0
    return f1,f2,f3,f4

def ConvertFromOPLSToRBDihedral(f1,f2,f3,f4):

    c0 = f2 + 0.5*(f1+f3)
    c1 = 0.5*(-f1+3*f3)
    c2 = -f2 + 4*f4
    c3 = -2*f3
    c4 = -4*f4
    c5 = 0.0
    c6 = 0.0
    return c0,c1,c2,c3,c4,c5,c6

