class AbstractDihedralType(object):
    __slots__ = ['atom1', 'atom2', 'atom3', 'atom4', 'improper']
    def __init__(self, atom1, atom2, atom3, atom4, improper=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.improper = improper

    def __eq__(self, dihedral_type):
        return ((self.atom1 == dihedral_type.atom1 and
                 self.atom2 == dihedral_type.atom2 and
                 self.atom3 == dihedral_type.atom3 and
                 self.atom4 == dihedral_type.atom4)
                or
                (self.atom1 == dihedral_type.atom4 and
                 self.atom2 == dihedral_type.atom3 and
                 self.atom3 == dihedral_type.atom2 and
                 self.atom4 == dihedral_type.atom1)
                and
                self.improper == dihedral_type.improper)

    def __hash__(self):
        if (self.atom1 == '*') and (self.atom4 == '*'):
            return hash(tuple([self.atom2, self.atom3,self.improper]))
        else:
            return hash(tuple([self.atom1, self.atom2, self.atom3, self.atom4, self.improper]))

