Example
=======

>>> Driver.initSystem("Solvated GMX")
>>> # Reading
>>> Driver.loadctools("system_GMX.top")
>>> Driver.loadStructure("system_GMX.gro")
>>> # Writing
>>> Driver.writeStructure("system_GMX_out.gro")
>>> Driver.writectools("system_GMX_out.top")
