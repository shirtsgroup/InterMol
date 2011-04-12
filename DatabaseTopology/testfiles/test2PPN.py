import Topology.Driver as Driver
import pdb

Driver.initSystem("Solvated 2PPN")
Driver.loadStructure("2PPN.gro")
Driver.loadTopology("2PPN.top")
Driver.writeStructure("2PPNout.gro")
Driver.writeTopology("2PPNout.top")
