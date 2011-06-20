import ctools.Driver as Driver
import pdb

Driver.initSystem("Solvated GMX")
Driver.loadTopology("system2_GMX.top")
#pdb.set_trace()
Driver.loadStructure("system2_GMX.gro")
Driver.writeStructure("system2_GMX_out.gro")
Driver.writeTopology("system2_GMX_out.top")

