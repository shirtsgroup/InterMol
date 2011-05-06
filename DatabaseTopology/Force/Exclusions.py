from Topology.Decorators import *

class Exclusions(object):

    def __init__(self, exclusions):
        """
        """
        if (type(exclusions).__name__ == 'list') and (len(exclusions) > 0):
            self.exclusions = exclusions
