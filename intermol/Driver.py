from os.path import splitext
import pdb

from intermol.System import System


def initSystem(name):
    """Initialize the System. This must be called prior to reading in anything

    Args:
        name (str): name of the system
    """
    System._sys = System(name)
    print "System initialized\n"

def load(*files):
    """
    """
    for f in files:
        try:
            filename = str(f)
        except:
            raise Exception("Arguments should be string filenames")

        extension = splitext(filename)[1].lower()

        if extension == '.gro':
            import intermol.GromacsExt.GromacsStructureParser as GromacsStructureParser
            print "Reading in Gromacs structure '{0}'...".format(filename)
            GromacsStructureParser.readStructure(filename)
            print "Structure loaded\n"

        elif extension == '.top':
            from intermol.GromacsExt.GromacsTopologyParser import GromacsTopologyParser
            print "Reading in Gromacs topology '{0}'...".format(filename)
            if not GromacsTopologyParser._GroTopParser:
                GromacsTopologyParser._GroTopParser = GromacsTopologyParser()
            GromacsTopologyParser._GroTopParser.parseTopology(filename)
            print "Topology loaded\n"

        elif extension == '.cms':
            from intermol.DesmondExt.DesmondParser import DesmondParser as DesmondParser
            DesmondParser = DesmondParser()
            print "Reading in Desmond structure '{0}'..."
            DesmondParser.readFile(filename)
            print "Sructure loaded\n"

        else:
            pdb.set_trace()
            raise Exception("{0} is not a supported file format".format(extension))

def write(*files):
    """
    """
    for f in files:
        try:
            filename = str(f)
        except:
            raise Exception("Arguments should be string filenames")

        extension = splitext(filename)[1].lower()

        if extension == '.gro':
            import intermol.GromacsExt.GromacsStructureParser as GromacsStructureParser
            print "Writing Gromacs structure file..."
            GromacsStructureParser.writeStructure(filename)

        elif extension == '.top':
            from intermol.GromacsExt.GromacsTopologyParser import GromacsTopologyParser
            print "Writing Gromacs topology file..."
            GromacsTopologyParser._GroTopParser.writeTopology(filename)

        elif extension == '.cms':
            from intermol.DesmondExt.DesmondParser import DesmondParser
            print "Writing Desmond topology file..."
            DesmondParser = DesmondParser()
            DesmondParser.writeFile(filename)

        elif extension == '.lmp':
            import intermol.LammpsExt.LammpsParser as LammpsParser
            print "Writing LAMMPS data & input files..."
            LammpsParser.writeData(filename)
            input_name = splitext(filename)[0] + '.input'
            print "Finished writing '{0}'".format(input_name)

        else:
            raise Exception("{0} is not a supported file format".format(extension))

        print "Finished writing '{0}'".format(filename)
