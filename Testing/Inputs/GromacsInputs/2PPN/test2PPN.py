import ctools.Driver as Driver
import pdb

Driver.initSystem("Solvated 2PPN")
Driver.loadTopology("2PPN.top")
Driver.loadStructure("2PPN.gro")
Driver.writeStructure("2PPNout.gro")
Driver.writeTopology("2PPNout.top")
