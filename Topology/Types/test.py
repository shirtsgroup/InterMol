from cctools.Types import *



newBondType1 = BondType(1, 2, 1, 5.0, 1.0)
newBondType2 = BondType(2, 3, 1, 4.0, 2.0)
newBondType3 = BondType(3, 4, 1, 3.0, 3.0)

BondTypes = set([newBondType1, newBondType2, newBondType3])

newBond = AbstractBondType(1, 2, 1)

if newBond in BondTypes:
    print 'Found it!'

