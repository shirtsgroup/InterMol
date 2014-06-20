from orderedset import OrderedSet
from hashmap import HashMap


class MoleculeType(object):
    """An abstract container for molecules of one type
    """
    def __init__(self, name):
        """Initialize the MoleculeType container

        Args:
            name (str): the name of the moleculetype to add
        """
        self.name = name
        self.moleculeSet = OrderedSet()
        self.bondForceSet = HashMap()
        self.pairForceSet = HashMap()
        self.angleForceSet = HashMap()
        self.dihedralForceSet = HashMap()
        self.torsiontorsionForceSet = HashMap()
        self.constraints = HashMap()
        self.exclusions = HashMap()
        self.settles = None
        self.nrexcl = None

    def addMolecule(self, molecule):
        """Add a molecule into the moleculetype container

        Args:
            molecule (Molecule): the molecule to append
        """
        self.moleculeSet.add(molecule)

    def removeMolecule(self, molecule):
        """Remove a molecule from the system.

        Args:
            molecule (Molecule): remove a molecule from the moleculeType
        """
        self.moleculeSet.remove(molecule)

    def getMolecule(self, molecule):
        """Get a molecule from the system

        Args:
            molecule (Molecule): retrieve an equivalent molecule from the moleculetype
        """
        return get_equivalent(self.moleculeSet, molecule, False)

    def addForce(self, force):
        """Add a forces to the moleculeType

        Args:
            forces (AbstractForce): Add a forces or contraint to the moleculeType
        """
        self.forceSet.add(force)

    def removeForce(self, force):
        """Remove a forces from the moleculeType

        Args:
            forces (AbstractForce): Remove a forces from the moleculeType
        """
        self.forceSet.remove(force)

    def getForce(self, force):
        """Get a forces from the moleculeType

        Args:
            forces (AbstractForce): Retrieve a forces from the moleculeType
        """
        return get_equivalent(self.forceSet, force, False)

    def setNrexcl(self, nrexcl):
        """Set the nrexcl

        Args:
            nrexcl (int): the value for nrexcl
        """
        self.nrexcl = nrexcl

    def getNrexcl(self):
        """Gets the nrexcl
        """
        return self.nrexcl
