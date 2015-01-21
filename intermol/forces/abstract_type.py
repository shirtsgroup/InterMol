class AbstractType(object):
    def __init__(self):
        pass

    def __repr__(self):
        """Print the object and all of its "private" attributes. """
        attributes = ["{0}={1}".format(x, getattr(self, x, "Undefined")) for x in dir(self)
                      if not (x.endswith("__") or hasattr(getattr(self, x, "Undefined"), '__call__'))]
        return "{0}({1})".format(self.__class__.__name__, ", ".join(attributes))
