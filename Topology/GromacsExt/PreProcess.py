from GromacsTopologyParser import *

Parser = GromacsTopologyParser()

for line in Parser.preprocess('2PPN.top'):
    print line

