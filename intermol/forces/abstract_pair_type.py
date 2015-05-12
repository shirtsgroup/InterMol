from intermol.forces.abstract_type import AbstractType


class AbstractPairType(AbstractType):
    __slots__ = ['bondingtype1', 'bondingtype2', 'scaleLJ', 'scaleQQ', 'long']

    pairlist = ['lj_c',
                'lj_sigeps',
                'ljq_c',
                'ljq_sigeps',
                'lj_default',
                'ljq_default'
                ]

    subtypes = []

    def __init__(self, bondingtype1, bondingtype2, scaleLJ=None, scaleQQ=None,
                 long=False):
        super(AbstractPairType, self).__init__()
        self.bondingtype1 = bondingtype1
        self.bondingtype2 = bondingtype2
        self.scaleLJ = scaleLJ
        self.scaleQQ = scaleQQ
        self.long = long
