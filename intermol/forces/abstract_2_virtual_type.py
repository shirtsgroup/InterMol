from intermol.forces.abstract_type import AbstractType
from intermol.forces.abstract_virtual_type import AbstractVirtualType


class Abstract2VirtualType(AbstractVirtualType):
    __slots__ = ['bondingtype1', 'bondingtype2', 'bondingtype3', 'placeholder']

    def __init__(self, bondingtype1, bondingtype2, bondingtype3, placeholder=False):
        super(AbstractVirtualType, self).__init__(placeholder=placeholder)
        self.bondingtype1 = bondingtype1
        self.bondingtype2 = bondingtype2
        self.bondingtype3 = bondingtype3
