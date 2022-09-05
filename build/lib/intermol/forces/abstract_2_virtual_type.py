from intermol.forces.abstract_type import AbstractType


class Abstract2VirtualType(AbstractType):
    __slots__ = ['bondingtype1', 'bondingtype2', 'bondingtype2', 'placeholder']

    def __init__(self, bondingtype1, bondingtype2, bondingtype3, placeholder):
        super(Abstract2VirtualType, self).__init__()
        self.bondingtype1 = bondingtype1
        self.bondingtype2 = bondingtype2
        self.bondingtype3 = bondingtype3
        self.placeholder = placeholder
