import Topology.Driver as Driver
import clearDB
import pdb

Driver.initSystem("Solvated GMX", "stardock.cs.virginia.edu", "spring2011")
Driver.loadTopology("system2_GMX.top")
#pdb.set_trace()
#Driver.loadStructure("system2_GMX.gro")
#Driver.writeStructure("system2_GMX_out.gro")
#Driver.writeTopology("system2_GMX_out.top")

