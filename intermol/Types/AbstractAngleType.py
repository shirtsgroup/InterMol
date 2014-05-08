class AbstractAngleType(object):
    __slots__ = ['atom1', 'atom2', 'atom3', 'type']
    def __init__(self, atom1, atom2, atom3, type):
        """An abstract representation of a generic angle type.

        Args:
            atom1 (int): Index of the first atom involved in the angle

            atom2 (int): Index of the center atom involved in the angle

            atom3 (int): Index of the final atom involved in the angle

            type  (int): A number specifying the type of the angle (follows GROMACS convention)

        >>> __init__(atom1=1, atom2=2, atom3=3, type=1)
        """
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.type = type


    def __eq__(self, object):
        if ((self.atom1 == object.atom1
                and self.atom2 == object.atom2 
                and self.atom3 == object.atom3) 
                or
                (self.atom1 == object.atom3 
                and self.atom2 == object.atom2 
                and self.atom3 == object.atom1)): 
            return True
        else:
            return False

    def __hash__(self):
        return hash(tuple([self.atom1, self.atom2, self.atom3, self.type]))

