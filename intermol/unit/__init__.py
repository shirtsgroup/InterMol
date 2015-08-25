"""
Physical quantities with units for dimensional analysis and automatic unit
conversion.
"""
from __future__ import division

__docformat__ = "epytext en"

__author__ = "Christopher M. Bruns"
__copyright__ = "Copyright 2010, Stanford University and Christopher M. Bruns"
__credits__ = []
__license__ = "MIT"
__maintainer__ = "Christopher M. Bruns"
__email__ = "cmbruns@stanford.edu"

# This code is copied from the `simtk.unit` package distributed as a standalone
# package (https://pypi.python.org/pypi/simtk.unit/) and as part of OpenMM
# (https://simtk.org/home/openmm) and as part of ParmEd
# (https://github.com/ParmEd/ParmEd)

# When parmed can be imported, the unit package will be taken from there.
# Otherwise, the implementation here will be used. This way, the
# `intermol.unit` package can be used interchangeably with ParmEd

try:
    from parmed.unit import *
except ImportError:
    from intermol.unit.unit import Unit, is_unit
    from intermol.unit.quantity import Quantity, is_quantity
    from intermol.unit.unit_math import *
    from intermol.unit.unit_definitions import *
    from intermol.unit.constants import *
