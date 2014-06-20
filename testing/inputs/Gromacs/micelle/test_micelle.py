import ctools.Driver as Driver
import pdb
Driver.initSystem("Solvated Micelle")

Driver.loadTopology("micelle.top")
Driver.loadStructure("micelle.gro")

Driver.writeStructure("micelle_out.gro")
Driver.write_topology("micelle_out.top")


