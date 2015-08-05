import collections
from copy import deepcopy


class OrderedSet(collections.Set):

    def __init__(self, iterable=()):
        self.d = collections.OrderedDict.fromkeys(iterable)

    def add(self, key):
        self.d[key] = None

    def discard(self, key):
        del self.d[key]

    def difference_update(self, *args, **kwargs):
        intersection = set(self.d.keys()).intersection(args[0])
        self.intersection_update(intersection)

    def intersection_update(self, *args, **kwargs):
        for part in args[0]:
            del self.d[part]

    def __len__(self):
        return len(self.d)

    def __contains__(self, element):
        return element in self.d

    def __iter__(self):
        return self.d.__iter__()

    def __le__(self, other):
        if not isinstance(other, collections.Set):
            return NotImplemented
        if len(self) > len(other):
            return False

        for e1, e2 in zip(self, other):
            if e1 != e2:
                return False
        return True

    def __repr__(self):
        class_name = self.__class__.__name__
        if not self:
            return '{0!s}()'.format(class_name)
        return '{0!s}({1!r})'.format(class_name, list(self))

    def __deepcopy__(self, memo):
        result = OrderedSet()
        for elt in self:
            result.add(deepcopy(elt,memo))
        return result
