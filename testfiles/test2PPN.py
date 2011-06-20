import ctools.Driver as Driver
import pdb

Driver.initSystem("Solvated 2PPN")
Driver.loadStructure("2PPN.gro")
Driver.loadctools("2PPN.top")
Driver.writeStructure("2PPNout.gro")
Driver.writectools("2PPNout.top")
