import pdb
import MySQLdb
from Topology import *
from Topology.GromacsExt.GromacsTopologyParser import GromacsTopologyParser
from Topology.System import System
from Topology.DBClass import DBClass

def initSystem(name, host, username,password):
    DBClass._db = MySQLdb.connect(host=host, user=username, passwd=password, db='cs4750pmc8p')
    System._sys = System(name)
    GromacsTopologyParser._GroTopParser = GromacsExt.GromacsTopologyParser()
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
            GromacsExt.readStructure(filename, sys, db)
        
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
            print "Writing Gromacs structure file..."
            GromacsExt.writeStructure(sys, filename)
        
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

        elif extension.lower() == 'array':
            System._sys.writeXML(filename)

        else:
            sys.stderr.write(" Error: '%s' is not a supported topology file format\n" %(extension.lower()))

        print 'Finished writing "%s"\n' %(filename)
 
