import ctools.Driver as Driver
import pdb
Driver.initSystem("Solvated Micelle")

Driver.loadTopology("micelle.top")
Driver.loadStructure("micelle.gro")

Driver.writeStructure("micelle_out.gro")
Driver.writeTopology("micelle_out.top")


