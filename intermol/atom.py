class Atom(object):
    """  """
    def __init__(self, index, name=None, residue_index=-1, residue_name=None):
        """Create an Atom object

        Args:
            index (int): index of atom in the molecule
            name (str): name of the atom (eg., N, CH)
            residue_index (int): index of residue in the molecule
            residue_name (str): name of the residue (eg., THR, CYS)
        """
        self.index = index
        self.name = name
        self.residue_index = residue_index
        self.residue_name = residue_name

        self._position = list()
        self._velocity = list()
        self._force = list()

        self._atomtype = dict()
        self.bondingtype = None
        self.atomic_number = None
        self.cgnr = None
        self._mass = dict()
        self._charge = dict()
        self.ptype = "A"
        self._sigma = dict()
        self._epsilon = dict()

    @property
    def atomtype(self):
        return self._atomtype

    @atomtype.setter
    def atomtype(self, index_atomtype):
        """Sets the atomtype

        Args:
            index_atomtype (tuple): A or B state and atomtype
        """
        try:
            idx, val = index_atomtype
        except ValueError:
            raise ValueError("Pass an iterable with two items.")
        else:
            self._atomtype[idx] = val

    @property
    def sigma(self):
        return self._sigma

    @sigma.setter
    def sigma(self, index_sigma):
        """Sets the sigma

        Args:
            index_sigma (tuple): A or B state and sigma
        """
        try:
            idx, val = index_sigma
        except ValueError:
            raise ValueError("Pass an iterable with two items.")
        else:
            self._sigma[idx] = val

    @property
    def epsilon(self):
        return self._epsilon

    @epsilon.setter
    def epsilon(self, index_epsilon):
        """Sets the epsilon

        Args:
            index_epsilon (tuple): A or B state and epsilon
        """
        try:
            idx, val = index_epsilon
        except ValueError:
            raise ValueError("Pass an iterable with two items.")
        else:
            self._epsilon[idx] = val

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, index_mass):
        """Sets the mass

        Args:
            index_mass (tuple): A or B state and mass
        """
        try:
            idx, val = index_mass
        except ValueError:
            raise ValueError("Pass an iterable with two items.")
        else:
            self._mass[idx] = val

    @property
    def charge(self):
        return self._charge

    @charge.setter
    def charge(self, index_charge):
        """Sets the charge

        Args:
            index_charge (tuple): A or B state and charge
        """
        try:
            idx, val = index_charge
        except ValueError:
            raise ValueError("Pass an iterable with two items.")
        else:
            self._charge[idx] = val

    @property
    def position(self):
        """Return the cartesian coordinates of the atom """
        return self._position

    @position.setter
    def position(self, xyz):
        """Sets the position of the atom

        Args:
            xyz (list, float): x, y and z coordinates
        """
        self._position = xyz

    @property
    def velocity(self):
        """Return the velocity of the atom"""
        return self._velocity

    @velocity.setter
    def velocity(self, vxyz):
        """Sets the velocity of the atom

        Args:
            vxyz (list, float): x-, y- and z-directed velocity
        """
        self._velocity = vxyz

    @property
    def force(self):
        """Return the force on the atom """
        return self._force

    @force.setter
    def force(self, fxyz):
        """Sets the force on the atom

        Args:
            fxyz (list, float): x-, y- and z-directed force
        """
        self._force = fxyz

    def __repr__(self):
        return 'Atom{0}({1}, {2})'.format(id(self), self.index, self.name)

    def __str__(self):
        return 'Atom({0}, {1})'.format(self.index, self.name)
