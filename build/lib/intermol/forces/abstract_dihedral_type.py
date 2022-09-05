from intermol.forces.abstract_type import AbstractType


class AbstractDihedralType(AbstractType):
    __slots__ = ['bondingtype1', 'bondingtype2', 'bondingtype3', 'bondingtype4',
                 'improper']

    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4,
                 improper=False):
        super(AbstractDihedralType, self).__init__()
        self.bondingtype1 = bondingtype1
        self.bondingtype2 = bondingtype2
        self.bondingtype3 = bondingtype3
        self.bondingtype4 = bondingtype4
        self.improper = improper
