class AbstractDihedral(object):
    def __init__(self, atom1, atom2, atom3, atom4, multiplicity):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.multiplicity = multiplicity


    def __eq__(self, object):
        if ((self.atom1 == object.atom1
                and self.atom2 == object.atom2
                and self.atom3 == object.atom3
                and self.atom4 == object.atom4
                and self.multiplicity == object.multiplicity)
                or
                (self.atom1 == object.atom4
                and self.atom2 == object.atom3
                and self.atom3 == object.atom2
                and self.atom4 == object.atom1
                and self.multiplicity == object.multiplicity)):
            return True
        else:
            return False


    def __hash__(self):
        return hash(tuple([self.atom1, self.atom2, self.atom3, self.atom4, self.multiplicity]))

    def __repr__(self):
        return str(self.atom1) + " " + str(self.atom2) + " " + str(self.atom3) + " " + str(self.atom4)
