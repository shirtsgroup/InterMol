import Topology.Driver as Driver
import pdb
Driver.initSystem("Solvated GMX")

Driver.loadTopology("system_GMX.top")
Driver.loadStructure("system_GMX.gro")

Driver.writeStructure("system_GMX_out.gro")
Driver.writeTopology("system_GMX_out.top")


