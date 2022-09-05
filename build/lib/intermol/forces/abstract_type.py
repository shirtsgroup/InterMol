class AbstractType(object):

    def __repr__(self):
        """Print the object and all of its non-magic attributes. """
        attributes = ["{0}={1}".format(x, getattr(self, x)) for x in dir(self)
                      if not (x.startswith('__') or x.endswith('__'))]
        printable_attributes = ', '.join(attributes)
        return "{0}({1})".format(self.__class__.__name__, printable_attributes)
