class AbstractDihedralType(object):
    def __init__(self, atom1, atom2, atom3, atom4, type):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.type = type


    def __eq__(self, object):
        if ((self.atom1 == object.atom1
                and self.atom2 == object.atom2
                and self.atom3 == object.atom3
                and self.atom4 == object.atom4)
                or
                (self.atom1 == object.atom4
                and self.atom2 == object.atom3
                and self.atom3 == object.atom2
                and self.atom4 == object.atom1)):
            return True
        else:
            return False


    def __hash__(self):
        if (self.atom1 == '*') and (self.atom4 == '*'):
            return hash(tuple([self.atom2, self.atom3, self.type]))
        else:
            return hash(tuple([self.atom1, self.atom2, self.atom3, self.atom4, self.type]))

