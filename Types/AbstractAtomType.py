import sys
sys.path.append('..')
from Decorators import Decorators

class AbstractAtomType(object):

    def __init__(self, 
            atomtype, 
            bondtype = None,
            Z =  None, 
            mass = None, 
            charge = None, 
            ptype = None):
        """An abstract representation of a generic atom type.

        Args:
            atomtype (str): The type of the atom
            
            bondtype (str): The type of the bond the atom is involved in

            Z (int): The atomic number of the atom

            mass (float): The mass of the atom in 'amu'

            charge (float): The charge of the atom in 'elementary charge units'

            ptype (str): The free energy state of the type

        >>> __init__(atomtype='H0', bondtype='H0', Z=1, mass=1.0080, charge=0.2329, ptype='A')
        """

        self.atomtype = atomtype
        self.bondtype = bondtype
        self.Z = Z
        self.mass = mass
        self.charge = charge
        self.ptype = ptype

    def __eq__(self, object):
        return self.atomtype == object.atomtype

    def __hash__(self):
        return hash(self.atomtype[0])

