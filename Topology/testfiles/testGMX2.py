import Topology.Driver as Driver
import pdb

Driver.initSystem("Solvated GMX")
Driver.loadStructure("system2_GMX.gro")
Driver.loadTopology("system2_GMX.top")
Driver.writeStructure("system2_GMX_out.gro")
Driver.writeTopology("system2_GMX_out.top")

