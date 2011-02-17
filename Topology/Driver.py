import pdb
from  cctools import *

sys = None
GroTopParser = None 

def initSystem(name):
    global sys, GroTopParser
    sys = System(name)
    GroTopParser = GromacsTools.GromacsTopologyParser(sys)

    print "System initialized\n"

def loadStructure(*files):

    for file in files:
        try:
            filename = str(file)
        except:
            sys.stderr.write("ERROR(loadStructure): Invalid argument! Arguments should be string filenames\n")

        extension = file.split(".")[-1]
       
        if extension.lower() == 'gro':
            print 'Reading in Gromacs structure "%s"...' %(filename)
            value = GromacsTools.readStructure(filename)
            sys.name = value[0]
            v1x = value[1]
            v2x = value[2]
            v3x = value[3]
            v1y = value[4]
            v2y = value[5]
            v3y = value[6]
            v1z = value[7]
            v2z = value[8]
            v3z = value[9]
            sys.setBoxVector(v1x, v2x, v3x, v1y, v2y, v3y, v1z, v2z, v3z)
            for molecule in value[10]:
                sys.addMolecule(molecule, molecule.name) 

            sys.atoms = value[11]

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
            print 'Reading in Gromacs topology "%s"...' %(filename)
            global GroTopParser
            GroTopParser.parseTopology(filename)
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
            print "Writing Gromacs structure file..."
            GromacsTools.writeStructure(sys, filename)
        
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
            print "Writing Gromacs topology file..."
            GroTopParser.writeTopology(filename)

        elif extension.lower() == '':
            pass

        else:
            sys.stderr.write(" Error: '%s' is not a supported topology file format\n" %(extension.lower()))

        print 'Finished writing "%s"\n' %(filename)
