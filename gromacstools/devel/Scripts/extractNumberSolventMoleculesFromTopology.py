from mmtools.gromacstools.devel.Helpers.TopologyFile import *

def extractNumberSolventMoleculesFromTopology(intopfile):
    """Extract the number of water molecules from a topology file."""
    
    top = TopologyFile(topfile=intopfile)
    # print top.nmols
    
    try: 
	numberWaterMolecules = top.nmols['SOL']
    except:
	raise TopologyFileFormatError

    return int(numberWaterMolecules)

