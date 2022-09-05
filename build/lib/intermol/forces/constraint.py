from intermol.decorators import *

class Constraint(object):

    @accepts_compatible_units(None, None, None, None, None, None, None, None, None,
            units.nanometers, units.nanometers, units.nanometers, units.nanometers,
            units.nanometers, units.nanometers, units.nanometers, units.nanometers,
            None)

    def __init__(self, atom1, atom2, length1, type, atom3=None, length2=None, atom4=None, length3=None, atom5=None, length4=None, atom6=None, length5=None, atom7=None, length6=None, atom8=None, length7=None, atom9=None, length8=None):
        """
        """
        self.type = type
        if type=='HOH':
            self.n = 2
            self.atom1 = atom1
            self.atom2 = atom2
            self.atom3 = atom3
            self.length1 = length1
            self.length2 = length2
            self.length3 = length3
        elif str(type)[0:2]=='AH':
            self.n = int(list(type)[-1])
            if self.n >= 1:
                self.atom1 = atom1
                self.atom2 = atom2
                self.length1 = length1
            if self.n >= 2:
                self.atom3 = atom3
                self.length2 = length2
                self.length3 = None
            if self.n >= 3:
                self.atom4 = atom4
                self.length3 = length3
            if self.n >= 4:
                self.atom5 = atom5
                self.length4 = length4
            if self.n >= 5:
                self.atom6 = atom6
                self.length5 = length5
            if self.n >= 6:
                self.atom7 = atom7
                self.length6 = length6
            if self.n >= 7:
                self.atom8 = atom8
                self.length7 = length7
            if self.n == 8:
                self.atom9 = atom9
                self.length8 = length8

    def __repr__(self):
        return  "{0}({1})".format(self.__class__.__name__, ", ".join(["{0}={1}".format(x, getattr(self, x, "Undefined")) for x in dir(self) if not(x.endswith("__") or hasattr(getattr(self, x, "Undefined"), '__call__'))]))

