import sys
from ctools.System import System

def initSystem(name):
    """Initialize the System. This must be called prior to reading in anything
    
    Args:
        name (str): name of the system
    """
    System._sys = System(name)
    print "System initialized\n"

def loadStructure(*files):

    for file in files:
        try:
            filename = str(file)
        except:
            sys.stderr.write("ERROR(loadStructure): Invalid argument! Arguments should be string filenames\n")

        extension = file.split(".")[-1]
       
        if extension.lower() == 'gro':
            import ctools.GromacsExt.GromacsStructureParser as GromacsStructureParser
            print 'Reading in Gromacs structure "%s"...' %(filename)
            GromacsStructureParser.readStructure(filename)
        
        elif extension.lower() == '':
            pass

        else:
            sys.stderr.write(" Error: '%s' is not a supported structure file format\n" %(extension.lower()))
        
        print "Structure loaded\n"
        
def loadTopology(*files):
     
    for file in files:
        try:
            filename = str(file)
        except:
            sys.stderr.write("ERROR(loadTopology): Invalid argument! Arguments should be string filenames\n")

        extension = file.split(".")[-1]


        if extension.lower() == 'top':
            from ctools.GromacsExt.GromacsTopologyParser import GromacsTopologyParser
            print 'Reading in Gromacs topology "%s"...' %(filename)
            if not GromacsTopologyParser._GroTopParser:
                GromacsTopologyParser._GroTopParser = GromacsTopologyParser()
            GromacsTopologyParser._GroTopParser.parseTopology(filename)
        elif extension.lower() == '':
            pass

        else:
            sys.stderr.write(" Error: '%s' is not a supported topology file format\n" %(extension.lower()))
        
        print "Topology loaded\n"

def writeStructure(*files):
    
    for file in files:
        try:
            filename = str(file)
        except:
            sys.stderr.write("ERROR(writeStructure): Invalid argument! Arguments should be string filenames\n")

        extension = file.split(".")[-1]

        if extension.lower() == 'gro':
            import ctools.GromacsExt.GromacsStructureParser as GromacsStructureParser
            print "Writing Gromacs structure file..."
            GromacsStructureParser.writeStructure(filename)
        
        elif extension.lower() == 'pdb':
            print "Writing PDB structure file..."

        else:
            sys.stderr.write(" Error: '%s' is not a supported structure file format\n" %(extension.lower()))

        print 'Finished writing "%s"\n' %(filename)


def writeTopology(*files):

    for file in files:
        try:
            filename = str(file)
        except:
            sys.stderr.write("ERROR(writeTopology): Invalid argument! Arguments should be string filenames\n")

        extension = file.split(".")[-1]

        if extension.lower() == 'top':
            from ctools.GromacsExt.GromacsTopologyParser import GromacsTopologyParser
            print "Writing Gromacs topology file..."
            GromacsTopologyParser._GroTopParser.writeTopology(filename)

        elif extension.lower() == '':
            pass

        else:
            sys.stderr.write(" Error: '%s' is not a supported topology file format\n" %(extension.lower()))

        print 'Finished writing "%s"\n' %(filename)
