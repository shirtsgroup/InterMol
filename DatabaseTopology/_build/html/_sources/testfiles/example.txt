Example
=======

>>> Driver.initSystem("Solvated GMX")
>>> # Reading
>>> Driver.loadTopology("system_GMX.top")
>>> Driver.loadStructure("system_GMX.gro")
>>> # Writing
>>> Driver.writeStructure("system_GMX_out.gro")
>>> Driver.writeTopology("system_GMX_out.top")
