from intermol.Decorators import *

class VSite2(object):
    
    def __init__(self, atom1, atom2, atom3, a):
        if atom1:
            self.atom1 = atom1
        if atom2:
            self.atom2 = atom2
        if atom3:
            self.atom3 = atom3
        self.a = a

    def getarameters(self):
        return(self.atom1, self.atom2, self.atom3, self.a)
