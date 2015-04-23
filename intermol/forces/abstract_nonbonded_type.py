from intermol.forces.abstract_type import AbstractType


class AbstractNonbondedType(AbstractType):
    __slots__ = ['atom1', 'atom2', 'type']

    def __init__(self, atom1, atom2, type):
        super(AbstractNonbondedType, self).__init__()
        self.atom1 = atom1
        self.atom2 = atom2
        self.type = type
        # TODO: refactor to avoid using builtin 'type'
