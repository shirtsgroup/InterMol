from intermol.forces.abstract_type import AbstractType


class AbstractBondType(AbstractType):
    __slots__ = ['bondingtype1', 'bondingtype2', 'order', 'c']

    def __init__(self, bondingtype1, bondingtype2, order=1, c=False):
        super(AbstractBondType, self).__init__()
        self.bondingtype1 = bondingtype1
        self.bondingtype2 = bondingtype2
        self.order = order  # Bond order. Desmond only.
        self.c = c  # Is the bond constrained or not? Desmond only.
