import copy
import sys
import pdb

from groupy.gbb import *
from groupy.system import *
from groupy.mdio import *
from groupy.general import *

names = ['system2_GMX']
for name in names:
    sys = Gbb()
    sys.load_lammps_data('{0}.lmp'.format(name))

    write_xyz(sys.xyz, sys.types, '{0}.xyz'.format(name))

