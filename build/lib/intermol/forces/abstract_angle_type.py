from intermol.forces.abstract_type import AbstractType


class AbstractAngleType(AbstractType):
    __slots__ = ['bondingtype1', 'bondingtype2', 'bondingtype3', 'c']

    def __init__(self, bondingtype1, bondingtype2, bondingtype3, c=False):
        """An abstract representation of a generic angle type. """
        super(AbstractAngleType, self).__init__()
        self.bondingtype1 = bondingtype1
        self.bondingtype2 = bondingtype2
        self.bondingtype3 = bondingtype3
        self.c = c   # Is the bond constrained or not? Desmond only.
