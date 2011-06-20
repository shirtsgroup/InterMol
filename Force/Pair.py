class Pair(object):
    def __init__(self, atom1, atom2, func):
        self.atom1 = atom1
        self.atom2 = atom2
        self.func = func

    def __eq__(self, object):
        if (self.atom1 == (object.atom1 or object.atom2)) and (self.atom2 == (object.atom2 or object.atom1)):
            return True
        else:
            return False

    def __hash__(self):
        return hash(tuple([self.atom1, self.atom2]))

