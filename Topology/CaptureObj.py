import pdb
class CaptureObj:
    'Object wrapper that remembers "other" for successful equality tests.'
    def __init__(self, obj):
        print 'Init'
        self.obj = obj
        self.match = obj

    def __eq__(self, other):
        print 'Capture Obj Equals'
        result = (self.obj == other)
        if result:
            pdb.set_trace()
            self.match = other
        return result

    def __getattr__(self, name):  # support hash() or anything else needed by __contains__
        print 'getattr'
        print name
        return getattr(self.obj, name)

def get_equivalent(container, item, default=None):
    '''Gets the specific container element matched by: "item in container".

    Useful for retreiving a canonical value equivalent to "item".  For example, a
    caching or interning application may require fetching a single representative
    instance from many possible equivalent instances).

    >>> get_equivalent(set([1, 2, 3]), 2.0)             # 2.0 is equivalent to 2
    2
    >>> get_equivalent([1, 2, 3], 4, default=0)
    0
    '''
    t = CaptureObj(item)
    if t in container:
        return t.match
    return default


