from cctools.Decorators import *

class Exclusions(object):

    def __init__(self, exclusions):
        """
        """
        if (type(exclusions).__name__ == 'list') and (len(exclusions) > 0):
            self.exclusions = exclusions

    def getForceParameters(self):
        return (self.exclusions)

    def __repr__(self):
        print self.exclusions

    def __str__(self):
        print self.exclusions
