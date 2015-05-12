class Bond(object):

    def __init__(self, atom1, atom2, bondtype=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.bondtype = bondtype