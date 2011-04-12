class HashMap(object):
    def __init__(self):
        self.map = {}

    def __len__(self):
        return len(self.map)
    
    def __contains__(self, key):
        return key in self.map
   
    def add(self, key):
        if key not in self.map:
            self.map[key] = key
    
    def remove(self,key): 
        if key in self.map:
            self.map.pop(key)

    def get(self, key):
        if key in self.map:
            return self.map[key]
        else:
            return False

    def itervalues(self):
        return self.map.values()        

    def __repr__(self):
        return repr(self.map)

    def __eq__(self, other):
        return set(self) == set(other)
                
