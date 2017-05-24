from collections import OrderedDict
import logging
import os

import parmed.unit as units

from intermol.utils import which, run_subprocess

CRM_PATH = ''

logger = logging.getLogger('InterMolLog')


to_canonical = {
    'BONDs': 'bond',

    'ANGLes': 'angle',
    'UREY-b': ['angle','urey-bradley'],
    
    'DIHEdrals': ['dihedral', 'proper'],
    'IMPRopers': ['dihedral', 'improper'],
    'CMAP': ['dihedral','cmap'],

    'HBONds': ['h-bond'],
    'IMHBnd': ['h-bond'],
    'VDWaals': ['vdw total', 'vdw (SR)'],
    'IMNBvdw': ['vdw total', 'vdw (LR)'],
    'EVDW': ['vdw total', 'vdw (SR)'],
    'RXNField': ['coulomb total','coulomb (SR)'],
    'ELEC': ['coulomb total', 'coulomb (SR)'],
    'IMELec': ['coulomb total', 'coulomb (LR)'],
    'EWKSum': ['coulomb total', 'coulomb (LR)'],
    'EWSElf': ['coulomb total', 'coulomb (LR)'],
    'EWEXcl': ['coulomb total', 'coulomb (LR)'],

    'ENERgy': 'potential'
}


def pick_crystal_type(box):
    ''' take a box vector and determine the crystal type (string output).
    has not been thoroughly debugged because of lack of examples
    '''
    a = box[0]
    b = box[1]
    c = box[2]
    alpha = box[3]
    beta = box[4]
    gamma = box[5]
    rectangular = (alpha == 90.0 and beta == 90.0 and gamma == 90.0)

    if rectangular:
        if a == b:
            if b == c:
                boxtype = 'cubic'
            else:
                boxtype = 'tetragonal'
        else:
            boxtype = 'orthorombic'

    elif alpha == gamma and alpha == 90:
        if alpha == beta:
            boxtype = 'orthorhombic'
        else:
            boxtype = 'monoclinic'

    elif a == b:
        if alpha == beta and alpha == 90.0 and gamma == 120.0:
            boxtype = 'hexagonal'
        elif a == c and alpha == 109.4712206344907 and alpha == beta and alpha == gamma:
            boxtype = 'octahedral' # this will be hard to match because of the degree of precision of alpha
        elif a == c and alpha == 60.0 and gamma == alpha and beta == 90.0:
            boxtype = 'rhdo'
        elif a == c and alpha == beta and gamma == beta and alpha < 120.0:
            boxtype = 'rhombohedral'
    else:
        boxtype = 'triclinic'

    return boxtype


def write_input_file(inpfile, psffile, rtfs, prms, strms,
                     boxtype, boxvecs, crdfile, charmm_settings, ignore_warnings=False):

    #annoyingly we need to write the charmm input file, which containes nonbonded interaction parameters,
    #and files thatare used.
    counter = 10
    increment = 10

    with open(inpfile, 'w') as charmm_inp:
        # will use relative paths because of length of line issues in charmm
        charmm_inp.write('! CHARMM Energy for %s\n' % os.path.relpath(inpfile))
        if (ignore_warnings):
            charmm_inp.write("BOMLEV  -1\n")

        for r in rtfs:
            charmm_inp.write('open read card unit %d name \"%s\"\nread rtf card unit %d\n' % (counter, os.path.relpath(r), counter))
            counter = counter + increment
        for p in prms:
            charmm_inp.write('open read card unit %d name \"%s\"\nread para card unit %d\n' % (counter, os.path.relpath(p), counter))
            counter = counter + increment
        for s in strms:
            charmm_inp.write('stream \"%s\"\n' % (os.path.relpath(s)))
        charmm_inp.write('read psf card name \"%s\"\n' % (os.path.relpath(psffile)))
        charmm_inp.write('read coor card name \"%s\"\nread coor card comp name \"%s\"\n' % (os.path.relpath(crdfile), os.path.relpath(crdfile)))
        charmm_inp.write('crystal define %s %.8f %.8f %.8f %.8f %.8f %.8f\ncrystal build noper 0\n' %  (boxtype,
               boxvecs[0], boxvecs[1], boxvecs[2], boxvecs[3], boxvecs[4], boxvecs[5]))
        # ! These segments are used for water and ions in bulk solvent
        #define bulks sele .not. (segid A .or. segid B) end
        # ! You may need to change these depending on how you plan to do recentering
        #image byseg sele .not. resname tip3 .and. .not. bulks end
        #image byres sele resname tip3 .or. bulks end
        charmm_inp.write("%s\n" % (charmm_settings))
        charmm_inp.write("energy\nstop")


def energies(inpfile, crm_path):
    """Compute single-point energies using CHARMM.

    Args:
        inpfile (str)
        crm_path (str):

    Returns:
        e_out:
        ener_xvg:
    """

    logger.info('Evaluating energy of {0}'.format(inpfile))

    # find the directory the input file is in.
    directory, _ = os.path.split(os.path.abspath(inpfile))

    # create files for output
    stdout_path = os.path.join(directory, 'charmm_stdout.txt')
    stderr_path = os.path.join(directory, 'charmm_stderr.txt')

    # delete previous stdout and stderr
    if os.path.isfile(stdout_path):
        os.remove(stdout_path)
    if os.path.isfile(stderr_path):
        os.remove(stderr_path)

    if not which(crm_path):
        raise IOError('Unable to find CHARMM executable (charmm).')
    if os.path.isdir(crm_path):
        charmm_bin = os.path.join(crm_path, 'charmm')
    else:
        charmm_bin = crm_path

    # run charmm - this assumes all files required in the input file are present
    cmd = [charmm_bin, '-i', inpfile]
    proc = run_subprocess(cmd, 'charmm', stdout_path, stderr_path)
    if proc.returncode != 0:
        logger.error('charmm failed. See %s' % stderr_path)

    # Extract energies from charmm output

    return _group_energy_terms(stdout_path)


def _group_energy_terms(mdout):
    """Parse CHARMM output file to extract and group the energy terms in a dict. """

    with open(mdout) as f:
        all_lines = f.readlines()

    startline = -1
    # find where the energy information starts
    for i, line in enumerate(all_lines):
        if line[0:9] == 'ENER ENR:':
            startline = i
            break

    if startline == -1:
        logger.error('No energy data (ENER ENR line) found in %s' % mdout)
        return

    energy_types = []
    energy_values = []

    for line in all_lines[startline:]:
        if line[0:4] == 'ENER':
            if '>' in line:
                # unforunately, we can' just split; we have to do it by fixed width
                startcol = 15
                colwidth = 13
                klim = (len(line)-startcol) // colwidth
                for k in range(klim):
                    startc = startcol+k*colwidth-1
                    endc = startc + colwidth
                    v = line[startc:endc]
                    energy_values.append(float(v)*units.kilocalories_per_mole)
            else:
                names = line.split()
                for n in names[2:]:
                    if n != "Eval#":
                        energy_types.append(n)

    e_out = OrderedDict(zip(energy_types, energy_values))
    return e_out, mdout
