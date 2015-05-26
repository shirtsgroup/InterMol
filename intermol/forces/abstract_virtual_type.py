from intermol.forces.abstract_type import AbstractType


class AbstractVirtualType(AbstractType):
    __slots__ = ['placeholder']

    def __init__(self, placeholder=False):
        super(AbstractVirtualType, self).__init__()
        self.placeholder = placeholder
