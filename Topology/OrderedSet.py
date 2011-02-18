import copy

class OrderedSet(object):
    def __init__(self):
        self.list = list()
        self.map = {}
        
    def __len__(self):
        return len(self.list)

    def __contains__(self, key):
        return key in self.map

    def __iter__(self):
        return iter(self.list)

    def __reverse__(self):
        return reverse(self.list)

    def add(self, key):
        if key not in self.map:
            self.map[key] = key
            self.list.append(key)

    def remove(self, key):
        if key in self.map:
            self.map.pop(key)
            self.list.remove(key)
    
    def get(self, key):
        if key in self.map:
            return copy.deepcopy(self.map[key])
        else:
            return False
     
    """
    def removeForUpdate(self, key):
        if key in self.map:
            return self.map.pop(key)
        
    def update(self, key):
        self.map[key] = key   
    """

    def pop(self, last=True):
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
              
