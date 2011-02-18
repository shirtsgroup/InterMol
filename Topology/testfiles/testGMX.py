import Topology.Driver as Driver
import pdb

Driver.initSystem("Solvated GMX")
Driver.loadStructure("system_GMX.gro")
Driver.loadTopology("system_GMX.top")
pdb.set_trace()
Driver.writeStructure("system_GMX_out.gro")
Driver.writeTopology("system_GMX_out.top")


