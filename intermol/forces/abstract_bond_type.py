from intermol.forces.abstract_type import AbstractType


class AbstractBondType(AbstractType):
    __slots__ = ['bondingtype1', 'bondingtype2', 'order', 'c']

    bondlist = ['harmonic',
                'fene',
                'fene_expandable',
                'connection',
                'morse',
                'nonlinear',
                'quartic',
                'quartic_breakable',
                'cubic',
                'harmonic_potential',
                'g96'
                ]

    subtypes = []

    def __init__(self, bondingtype1, bondingtype2, order=1, c=False):
        super(AbstractBondType, self).__init__()
        self.bondingtype1 = bondingtype1
        self.bondingtype2 = bondingtype2
        self.order = order  # Bond order. Desmond only.
        self.c = c  # Is the bond constrained or not? Desmond only.
