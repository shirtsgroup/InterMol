import warnings
from intermol.utils import which

for exe in ['lammps', 'lmp_mpi', 'lmp_serial', 'lmp_openmpi',
            'lmp_mac_mpi']:
    if which(exe):
        LMP_PATH = exe
        break
    else:
        warnings.warn('Found no LAMMPS executable.')

import intermol.lammps.lammps_driver
