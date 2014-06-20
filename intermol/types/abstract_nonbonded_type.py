class AbstractNonbondedType(object):
    __slots__ = ['atom1', 'atom2', 'type']
    def __init__(self, atom1, atom2, type):
        self.atom1 = atom1
        self.atom2 = atom2
        self.type = type

    def __eq__(self, nb_type):
        return ((self.atom1 == nb_type.atom1 and
                 self.atom2 == nb_type.atom2)
                or
                (self.atom2 == nb_type.atom2 and
                 self.atom3 == nb_type.atom1)
                and
                self.type == nb_type.type)

    def __hash__(self):
        return hash(tuple([self.atom1, self.atom2, self.type]))

