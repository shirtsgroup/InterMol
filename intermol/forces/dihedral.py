class Dihedral(object):

    def __init__(self, atom1, atom2, atom3, atom4, forcetype=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.forcetype = forcetype
