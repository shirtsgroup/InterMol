from intermol.forces.abstract_type import AbstractType


class AbstractAtomType(AbstractType):
    __slots__ = ['atomtype', 'bondtype', 'atomic_number', 'mass', 'charge', 'ptype']

    def __init__(self, atomtype, bondtype=None, atomic_number=None, mass=None,
                 charge=None, ptype=None):
        """An abstract representation of a generic atom type.

        Args:
            atomtype (str): The type of the atom
            bondingtype (str): The type of the bond the atom is involved in
            atomic_number (int): The atomic number of the atom
            mass (float): The mass of the atom in 'amu'
            charge (float): The charge of the atom in 'elementary charge units'
            ptype (str): The free energy state of the type
        """
        super(AbstractAtomType, self).__init__()
        self.atomtype = atomtype
        self.bondtype = bondtype
        self.atomic_number = atomic_number
        self.mass = mass
        self.charge = charge
        self.ptype = ptype
