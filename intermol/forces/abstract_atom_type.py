class AbstractAtomType(object):
    __slots__ = ['atomtype', 'bondtype', 'atomic_number', 'mass', 'charge', 'ptype']
    def __init__(self, atomtype, bondtype=None, atomic_number=None, mass=None,
            charge=None, ptype=None):
        """An abstract representation of a generic atom type.

        Args:
            atomtype (str): The type of the atom
            bondtype (str): The type of the bond the atom is involved in
            atomic_number (int): The atomic number of the atom
            mass (float): The mass of the atom in 'amu'
            charge (float): The charge of the atom in 'elementary charge units'
            ptype (str): The free energy state of the type

        >>> __init__(atomtype='H0', bondtype='H0', atomic_number=1, mass=1.0080, charge=0.2329, ptype='A')
        """

        self.atomtype = atomtype
        self.bondtype = bondtype
        self.atomic_number = atomic_number
        self.mass = mass
        self.charge = charge
        self.ptype = ptype


    def __repr__(self):
        return  "{0}({1})".format(self.__class__.__name__, ", ".join(["{0}={1}".format(x, getattr(self, x, "Undefined")) for x in dir(self) if not(x.endswith("__") or hasattr(getattr(self, x, "Undefined"), '__call__'))]))
