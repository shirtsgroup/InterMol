class AbstractAngleType(object):
    __slots__ = ['atom1', 'atom2', 'atom3']
    def __init__(self, atom1, atom2, atom3):
        """An abstract representation of a generic angle type."""
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3

    def __eq__(self, angle_type):
        return ((self.atom1 == angle_type.atom1 and
                 self.atom2 == angle_type.atom2 and
                 self.atom3 == angle_type.atom3)
                or
                (self.atom1 == angle_type.atom3 and
                 self.atom2 == angle_type.atom2 and
                 self.atom3 == angle_type.atom1))

    def __hash__(self):
        return hash(tuple([self.atom1, self.atom2, self.atom3]))

