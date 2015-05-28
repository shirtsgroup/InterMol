from intermol.forces.abstract_type import AbstractType


class AbstractNonbondedType(AbstractType):
    __slots__ = ['atom1', 'atom2', 'form']

    def __init__(self, atom1, atom2, form):
        super(AbstractNonbondedType, self).__init__()
        self.atom1 = atom1
        self.atom2 = atom2
        self.form = form
