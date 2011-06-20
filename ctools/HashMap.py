class HashMap(object):
    def __init__(self):
        """Initializes the HashMap class
        """
        self.map = {}

    def __len__(self):
        """Get the length of the container
        """
        return len(self.map)
    
    def __contains__(self, key):
        """Checks if key is in container

        Args:
            key: the key to check
        """
        return key in self.map
   
    def add(self, key):
        """Add a value to the container
        
        Args:
            key: the key to add
        """
        if key not in self.map:
            self.map[key] = key
    
    def remove(self,key):
        """Remove a key from the container
    
        Args:
            key: key to remove
        """ 
        if key in self.map:
            self.map.pop(key)

    def get(self, key):
        """Retrieve a key from the container

        Args:
            key: key to retrieve
        """
        if key in self.map:
            return self.map[key]
        else:
            return False

    def itervalues(self):
        """Return a list of values
        """
        return self.map.values()        

    def __repr__(self):
        return repr(self.map)

    def __eq__(self, other):
        return set(self) == set(other)
                
