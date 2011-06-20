from GromacsctoolsParser import *

Parser = GromacsctoolsParser()

for line in Parser.preprocess('2PPN.top'):
    print line

