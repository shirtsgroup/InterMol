class Exclusions(object):

    def __init__(self, exclusions):
        """
        """
        if isinstance(exclusions, list) and len(exclusions) > 0:
            self.exclusions = exclusions

    def get_parameters(self):
        return self.exclusions

    def __repr__(self):
        return self.exclusions

    def __str__(self):
        return self.exclusions
