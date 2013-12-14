import copy

class OrderedSet(object):
    def __init__(self):
        """Initialize the orderedSet. Essentially a map coupled with a list.
        """
        self.list = list()
        self.map = {}
        
    def __len__(self):
        """Return the size of the OrderedSet
        """
        return len(self.list)

    def __contains__(self, key):
        """Check if the orderedSet contains some key
        
        Args:
            key: the key to check
        """
        return key in self.map

    def __iter__(self):
        """Return an iterator over the container
        """
        return iter(self.list)

    def __reverse__(self):
        """Return a reverse iterator over the container
        """
        return reverse(self.list)

    def add(self, key):
        """Add a key to the container
        
        Args:
            key: key to add
        """
        if key not in self.map:
            self.map[key] = key
            self.list.append(key)

    def remove(self, key):
        """Remove a value from the container

        Args:
            key: key to remove:w
        """
        if key in self.map:
            self.map.pop(key)
            self.list.remove(key)
    
    def get(self, key):
        """Get a key from the container

        Args:
            key: key to retrieve from the container
        """
        if key in self.map:
            return copy.deepcopy(self.map[key])
        else:
            return False
     
    def pop(self, last=True):
        """Pop a value from the container

        Args:
            last (boolean): if true it pops the last element otherwise the first
        """
        if not self:
            raise KeyError('Set is empty')
        if last:
            key = self.list[-1]
        else:
            key = self.list[0]
        self.remove(key)
        return key    

    def __repr__(self):
        return repr(self.list)

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)

    def __getitem__(self, key):
        if isinstance(key, int):
            return copy.deepcopy(self.list[key])
        else:
            return copy.deepcopy(self.map[key])

    def __setitem__(self, key, value):
        if isinstance(key, int):
            temp = self.list[key]
            self.map.pop(temp)
            self.map[value] = value
            list[key] = value   
        else:
            idx = self.list.index(key)
            self.list[idx] = value
            self.map.pop(key)
            self.map[value] = value
              
