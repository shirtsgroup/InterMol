class Error(Exception):
    '''Base class for exceptions in InterMol'''
    pass

class ReadError(Error):
    '''Exception raised for errors in reading in input

    Attributes:
        type -- input file type
        file -- input file in which the error occurred
        msg  -- explanation of the error
    '''
    pass
#    def __init__(self, type, file, msg):
#        self.expr = expr
#        self.msg = msg

class WriteError(Error):
    '''Exception raised for errors in writing output file

    Attributes:
        type -- output file type
        file -- output file in which the error occurred
        msg  -- explanation of the error
    '''

    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg
