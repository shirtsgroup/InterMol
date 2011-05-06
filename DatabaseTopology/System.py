from Decorators import *
from MoleculeType import MoleculeType
from OrderedSet import OrderedSet
from OrderedDict import OrderedDict
from HashMap import HashMap
from DBClass import DBClass
import time
import pdb

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
        cursor.close()
        DBClass._db.commit()
        self.combinationRule = None
        self.atomsQueue = list()
        self.bondsQueue = list()
        self.anglesQueue = list()
        self.dihedralsQueue = list()
        self.prepChargeGroupQueue = list()
        self.chargeGroupQueue = list()
        self.pairsQueue = list()

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
        cursor.execute("""INSERT INTO Molecules (molIndex, molType) VALUES (%s, %s);""", (molecule.molIndex, molecule.molType))
        cursor.close()
        DBClass._db.commit()
        cursor = DBClass._db.cursor()
        cursor.execute("""SELECT molID FROM Molecules ORDER BY molID DESC LIMIT 0,1;""")
        result = cursor.fetchone()
        cursor.close()
        DBClass._db.commit()
        return result[0]

    def updateAtomStructure(self, atom):
        pass            

    def commitAtoms(self):
        cursor = DBClass._db.cursor()
        cursor.executemany("""INSERT INTO Atoms (resID, atom_index, atom_type, pType, atom_sigma, atom_epsilon, atom_z, atom_mass, atom_charge) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s);""",self.atomsQueue)
        cursor.close()
        DBClass._db.commit() 
        self.atomsQueue = list()
        self.commitChargeGroups()
    
    def duplicateAtoms(self, molID, atom):
        cursor = DBClass._db.cursor()
        cursor.execute("""SELECT * FROM Residues WHERE resIndex=%s AND resName=%s AND molID = %s LIMIT 0,1;""",(atom[20], atom[21], molID))    
        result = cursor.fetchone()
        cursor.close()
        DBClass._db.commit()
        if not result:
            cursor = DBClass._db.cursor()
            cursor.execute("""INSERT INTO Residues (resIndex, resName, molID) VALUES (%s, %s, %s);""",(atom[20], atom[21], molID))
            cursor.close()
            DBClass._db.commit()
            cursor = DBClass._db.cursor()
            cursor.execute("""SELECT resID FROM Residues ORDER BY resID DESC LIMIT 0,1;""")
            result = cursor.fetchone()
            cursor.close()
            DBClass._db.commit()
        #self.prepChargeGroupQueue.append((atom.cgnr, atom.atomNum, molID))
        self.atomsQueue.append((result[0], atom[3], atom[4], atom[5], atom[6], atom[7], atom[8], atom[9], atom[10]))
 
    def addAtom(self, molID, atom):
        cursor = DBClass._db.cursor()
        cursor.execute("""SELECT * FROM Residues WHERE resIndex=%s AND resName=%s and molID = %s LIMIT 0,1;""",(atom.resNum, atom.resName, molID))    
        result = cursor.fetchone()
        cursor.close()
        DBClass._db.commit()
        if not result:
            cursor = DBClass._db.cursor()
            cursor.execute("""INSERT INTO Residues (resIndex, resName, molID) VALUES (%s, %s, %s);""",(atom.resNum, atom.resName, molID))
            cursor.close()
            DBClass._db.commit()
            cursor = DBClass._db.cursor()
            cursor.execute("""SELECT resID FROM Residues ORDER BY resID DESC LIMIT 0,1;""")
            result = cursor.fetchone()
            cursor.close()
            DBClass._db.commit()
        self.prepChargeGroupQueue.append((atom.cgnr, atom.atomNum, molID))
        self.atomsQueue.append((result[0], atom.atomNum, atom.atomName, atom.ptype, atom.sigma[0], atom.epsilon[0], atom.Z, atom.mass[0], atom.charge[0]))
    
    def commitChargeGroups(self):
        for cgnr, index, molID in self.prepChargeGroupQueue:
            cursor = DBClass._db.cursor()
            cursor.execute("""SELECT atomID FROM (Atoms NATURAL JOIN Residues NATURAL JOIN Molecules) WHERE molID = %s AND atom_index = %s;""",(molID, index))  
            atomID = cursor.fetchone()[0]
            cursor.close()
            self.chargeGroupQueue.append((cgnr, atomID))
        self.prepChargeGroupQueue = list() 
        cursor = DBClass._db.cursor()
        cursor.executemany("""INSERT INTO ChargeGroup (chargeGroupID, charge_group_atomID ) VALUES (%s, %s);""",self.chargeGroupQueue)
        cursor.close()
        DBClass._db.commit() 
        self.chargeGroupQueue = list()

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
    
    def commitBonds(self):
        cursor = DBClass._db.cursor()
        cursor.executemany("""INSERT INTO Bonds (bond_atomID1, bond_atomID2, bond_func, bond_length, bond_k, bond_alpha, bond_beta) VALUES (%s, %s, %s, %s, %s, %s, %s);""", self.bondsQueue)
        cursor.close()
        DBClass._db.commit()
        self.bondsQueue = list()
     
    def addBond(self, molID, bond):
        cursor = DBClass._db.cursor()
        cursor.execute("""SELECT atomID FROM (Atoms NATURAL JOIN Residues NATURAL JOIN Molecules) WHERE molID = %s AND (atom_index = %s OR atom_index = %s);""", (molID, bond.atom1, bond.atom2))    
        atomID1 = cursor.fetchone()[0]
        atomID2 = cursor.fetchone()[0]
        cursor.close()
        DBClass._db.commit()
        
        self.bondsQueue.append((atomID1, atomID2, bond.func, bond._length._value, bond._k._value, bond._alpha._value, bond._beta._value))
       

    def commitAngles(self):
        cursor = DBClass._db.cursor()
        cursor.executemany("""INSERT INTO Angles (angle_atomID1, angle_atomID2, angle_atomID3, angle_func, angle_theta, angle_k, angle_c1, angle_c2, angle_c3, angle_c4, angle_c5, angle_k2) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);""", self.anglesQueue) 
        cursor.close()
        DBClass._db.commit()
        self.anglesQueue = list()

    def addAngle(self, molID, angle):
        cursor = DBClass._db.cursor()
        cursor.execute("""SELECT atomID FROM (Atoms NATURAL JOIN Residues NATURAL JOIN Molecules) WHERE molID = %s AND (atom_index = %s OR atom_index = %s or atom_index = %s);""", (molID, angle.atom1, angle.atom2, angle.atom3))
        atomID1 = cursor.fetchone()[0]
        atomID2 = cursor.fetchone()[0]
        atomID3 = cursor.fetchone()[0]
        cursor.close()
        DBClass._db.commit()
        self.anglesQueue.append((atomID1, atomID2, atomID3, angle.func, angle._theta._value, angle._k._value, angle._c1._value, angle._c2._value, angle._c3._value, angle._c4._value, angle._c5._value, angle._k2._value))

    def commitPairs(self):
        cursor = DBClass._db.cursor()
        cursor.executemany("""INSERT INTO Pairs (pairs_atomID1, pairs_atomID2, pairs_func) VALUES (%s, %s, %s);""", self.pairsQueue) 
        cursor.close()
        DBClass._db.commit()
        self.pairsQueue = list()
     
    def addPair(self, molID, pair):
        cursor = DBClass._db.cursor()
        cursor.execute("""SELECT atomID FROM (Atoms NATURAL JOIN Residues NATURAL JOIN Molecules) WHERE molID = %s AND (atom_index = %s OR atom_index = %s);""", (molID, pair.atom1, pair.atom2))
        atomID1 = cursor.fetchone()[0]
        atomID2 = cursor.fetchone()[0]
        cursor.close()
        DBClass._db.commit()
        self.pairsQueue.append((atomID1, atomID2, pair.func)) 

    def commitDihedrals(self):   
        cursor = DBClass._db.cursor()
        cursor.executemany("""INSERT INTO Dihedrals (dihedral_atomID1, dihedral_atomID2, dihedral_atomID3, dihedral_atomID4, dihedral_func, dihedral_phi, dihedral_k, dihedral_c1, dihedral_c2, dihedral_c3, dihedral_c4, dihedral_c5, dihedral_c6, dihedral_multiplicity) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);""", self.dihedralsQueue) 
        cursor.close()
        DBClass._db.commit()
        self.dihedralsQueue = list()

    def addDihedral(self, molID, dihedral):
        cursor = DBClass._db.cursor()
        cursor.execute("""SELECT atomID FROM (Atoms NATURAL JOIN Residues NATURAL JOIN Molecules) WHERE molID = %s AND (atom_index = %s OR atom_index = %s or atom_index = %s or atom_index = %s);""", (molID, dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4))
        atomID1 = cursor.fetchone()[0]
        atomID2 = cursor.fetchone()[0]
        atomID3 = cursor.fetchone()[0]
        atomID4 = cursor.fetchone()[0]
        cursor.close()
        DBClass._db.commit()
        self.dihedralsQueue.append((atomID1, atomID2, atomID3, atomID4, dihedral.func, dihedral._phi._value, dihedral._k._value, dihedral._c1._value, dihedral._c2._value, dihedral._c3._value, dihedral._c4._value, dihedral._c5._value, dihedral._c6._value, dihedral.multiplicity))
 
    def setSettles(self, molID, settles):
        cursor = DBClass._db.cursor()
        cursor.execute("""SELECT molType FROM Molecules WHERE molID = %s;""", (molID))
        molType = cursor.fetchone()[0]
        cursor.close()
        cursor = DBClass._db.cursor()
        cursor.execute("""UPDATE MoleculeTypes SET dOH = %s, dHH = %s WHERE moleculeType = %s;""", (settles.dOH, settles.dHH, molType))
        cursor.close()
        DBClass._db.commit()

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
        cursor.close()
        DBClass._db.commit()

    def writeXML(self, filename):
        cursor = DBClass._db.cursor()
        cursor.execute("""SELECT * FROM MoleculeOverview;""")
        result = cursor.fetchall()
        cursor.close()
        fout = open(filename, 'w')
        for item in result:    
            fout.write(str(item) + '\n')
        fout.close()
    
    def __str__(self):
        """String representation of a System object
        """
        return "System: " + self.name
        
    def __repr__(self):
        """String representation of a System object
        """
        return "System: " + self.name
