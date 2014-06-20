class AbstractDihedral(object):
    __slots__ = ['atom1', 'atom2', 'atom3', 'atom4', 'improper']

    def __init__(self, atom1, atom2, atom3, atom4, improper = False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.improper = improper
    def __eq__(self, object):
        if ((self.atom1 == object.atom1
                and self.atom2 == object.atom2
                and self.atom3 == object.atom3
                and self.atom4 == object.atom4)
                or
                (self.atom1 == object.atom4
                and self.atom2 == object.atom3
                and self.atom3 == object.atom2
                and self.atom4 == object.atom1)) and self.improper == object.improper:
            return True
        else:
            return False


    def __hash__(self):
        return hash(tuple([self.atom1, self.atom2, self.atom3, self.atom4, self.improper]))

    def __repr__(self):
        return str(self.atom1) + " " + str(self.atom2) + " " + str(self.atom3) + " " + str(self.atom4) + " " + str(self.improper)
