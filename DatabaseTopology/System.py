from Decorators import *
from MoleculeType import MoleculeType
from OrderedSet import OrderedSet
from OrderedDict import OrderedDict
from HashMap import HashMap
from DBClass import DBClass

class System(object):
    _sys = None
    
    def __init__(self, name = None):
        """Initialize a new System object. This must be run before the system can be used.

        Args:
            name (str): The name of the system

        >>> __init__(name='sysname')
        """
        cursor = DBClass._db.cursor()
        cursor.execute("""INSERT INTO Systems (sys_name, 
                                                genPairs,
                                                coulombCorrection,
                                                ljCorrection,
                                                nbFunc,
                                                combinationRule,
                                                box_v1x,
                                                box_v2x,
                                                box_v3x,
                                                box_v1y,
                                                box_v2y,
                                                box_v3y,
                                                box_v1z,
                                                box_v2z,
                                                box_v3z) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);""", (name, True, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
        cursor.close()
        DBClass._db.commit()
        
        cursor = DBClass._db.cursor()
        cursor.execute("""SELECT sysID FROM Systems WHERE sys_name=%s;""", (name)) 
        
        self.sysID = cursor.fetchone()[0]
        DBClass._db.commit()
        self.combinationRule = None


    def addMoleculeType(self, moleculeType, nrexcl):
        """Append a molecule into the System.

        Args:
            molecule (:py:class:`Topology.Molecule`): The molecule object to be appended
        """
        cursor = DBClass._db.cursor()
        cursor.execute("""INSERT INTO MoleculeTypes (moleculeType, sysID, nrexcl) VALUES (%s, %s, %s);""", (moleculeType, self.sysID, nrexcl))
        cursor.close()
        DBClass._db.commit() 
   
    def addMolecule(self, molecule):
        cursor = DBClass._db.cursor()
        cursor.execute("""INSERT INTO Molecules (sysID, molIndex, molType) VALUES (%s, %s, %s);""", (self.sysID, molecule.molIndex, molecule.molType))
        cursor.close()
        DBClass._db.commit()
        cursor = DBClass._db.cursor()
        cursor.execute("""SELECT molID FROM Molecules ORDER BY molID DESC LIMIT 0,1;""")
        result = cursor.fetchone()
        cursor.close()
        DBClass._db.commit()
        return result[0]
        
    def addAtom(self, atom):
        print atom.atomNum
        print atom.mass
        print atom.charge
        cursor = DBClass._db.cursor()
        cursor.execute("""SELECT * FROM Residues WHERE resIndex=%s AND resName=%s LIMIT 0,1;""",(atom.resNum, atom.resName))    
        result = cursor.fetchone()
        cursor.close()
        DBClass._db.commit()
        if not result:
            cursor = DBClass._db.cursor()
            cursor.execute("""INSERT INTO Residues (resIndex, resName, molID) VALUES (%s, %s, %s);""",(atom.resNum, atom.resName, atom.molID))
            cursor.close()
            DBClass._db.commit()
            cursor = DBClass._db.cursor()
            cursor.execute("""SELECT resID FROM Residues ORDER BY resID DESC LIMIT 0,1;""")
            result = cursor.fetchone()
            cursor.close()
            DBClass._db.commit()
        cursor = DBClass._db.cursor()
        cursor.execute("""INSERT INTO Atoms (resID, atom_index, atom_type, pType, atom_sigma, atom_epsilon, atom_z, atom_mass, atom_charge) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s);""",(result[0], atom.atomNum, atom.atomName, atom.ptype, atom.sigma[0], atom.epsilon[0], atom.Z, atom.mass[0], atom.charge[0]))
        cursor.close()
        DBClass._db.commit() 
 
    def addAtomTypes(self,atomtype):
        cursor = DBClass._db.cursor()        
        cursor.execute("""INSERT INTO AtomTypes (atomType, 
                                                bondType, 
                                                z, 
                                                mass,
                                                charge,
                                                ptype,
                                                sigma,
                                                epsilon) VALUES (%s, %s, %s, %s, %s, %s, %s, %s);""",(atomtype.atomtype, 
                    atomtype.bondtype, 
                    atomtype.Z, 
                    atomtype.mass._value, 
                    atomtype.charge._value, 
                    atomtype.ptype, 
                    atomtype.sigma._value, 
                    atomtype.epsilon._value))
        cursor.close()
        DBClass._db.commit()

    def getAtomType(self, atomtype):
        cursor = DBClass._db.cursor()
        cursor.execute("""SELECT bondType, z, mass, charge, ptype, sigma, epsilon FROM AtomTypes WHERE atomType=%s;""", (atomtype))
        result = cursor.fetchone()
        cursor.close()
        DBClass._db.commit()
        return result
    
    @accepts_compatible_units(units.nanometers, units.nanometers, units.nanometers, units.nanometers, units.nanometers, units.nanometers,units.nanometers, units.nanometers, units.nanometers)
    def setBoxVector(self, v1x, v2x, v3x, v1y, v2y, v3y, v1z, v2z, v3z):
        """Sets the boxvector for the system. Assumes the box vector is in the correct form. [[v1x,v2x,v3x],[v1y,v2y,v3y],[v1z,v2z,v3z]]
        """
        cursor = DBClass._db.cursor()
        cursor.execute("""UPDATE Systems SET v1x=%s, 
                                            v2x=%s, 
                                            v3x=%s, 
                                            v1y=%s, 
                                            v2y=%s, 
                                            v3y=%s, 
                                            v1z=%s, 
                                            v2z=%s, 
                                            v3z=%s 
                                            WHERE sysID=%s;""",(v1x._value, v2x._value, v3x._value, v1y._value, v2y._value, v3y._value, v1z._value, v2z._value, v3z._value, self.sysID))
        cursor.close()
        DBClass._db.commit()

    def setNBFunc(self, nbfunc):
        """Sets the nbfunc
        """
        cursor = DBClass._db.cursor()
        cursor.execute("""UPDATE Systems SET nbFunc=%s WHERE sysID=%s;""", (nbfunc, self.sysID))        
        cursor.close()
        DBClass._db.commit()        
    
    def setCombinationRule(self,combinationRule):
        self.combinationRule = combinationRule
        cursor = DBClass._db.cursor()
        cursor.execute("""UPDATE Systems SET combinationRule=%s WHERE sysID=%s;""",(combinationRule, self.sysID))
        cursor.close()
        DBClass._db.commit()

    def setGenPairs(self, genPairs):
        cursor = DBClass._db.cursor()
        cursor.execute("""UPDATE Systems SET genPairs=%s WHERE sysID=%s;""",(genPairs, self.sysID))
        cursor.close()
        DBClass._db.commit()

    def setLJCorrection(self, ljCorrection):
        cursor = DBClass._db.cursor()
        cursor.execute("""UPDATE Systems SET ljCorrection=%s WHERE sysID=%s;""",
(ljCorrection, self.sysID))
        cursor.close()
        DBClass._db.commit()

    def setCoulombCorrection(self, coulombCorrection):
        cursor = DBClass._db.cursor()
        cursor.execute("""UPDATE Systems SET coulombCorrection=%s WHERE sysID=%s;""",(coulombCorrection, self.sysID))

    def __str__(self):
        """String representation of a System object
        """
        return "System: " + self.name
        
    def __repr__(self):
        """String representation of a System object
        """
        return "System: " + self.name
