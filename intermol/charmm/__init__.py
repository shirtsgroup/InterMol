from collections import OrderedDict
import logging
import os

import simtk.unit as units

from intermol.utils import which, run_subprocess

CRM_PATH = ''

logger = logging.getLogger('InterMolLog')

# --------- energy evaluation methods ---------- #
# this dict hasn't been upated yet:
key_dict = {'ENERgy': 'Potential',
            'BONDs': 'Bond',
            'ANGLes': 'Angle',
            'DIHEdrals': 'Proper Dih.',
            'IMPRopers': 'Improper Dih.',
            'EEL': 'Coulomb',
            'VDWaals': 'LJ (SR)',
            'EVDW': 'Disper. corr.'
            }

def remove_zero_terms(sdict):
    # remove entries that are all zero
    for key in sdict:
        if sdict[key]._value == 0.0:
            sdict.pop(key)

def standardize_key(sdict,in_key):
    if in_key in key_dict:
        sdict[key_dict[in_key]] = sdict[in_key]
        sdict.pop(in_key)

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
                     boxtype, boxvecs, crdfile, charmm_settings):

    #annoyingly we need to write the charmm input file, which containes nonbonded interaction parameters,
    #and files thatare used.
    counter = 10
    increment = 10
    with open(inpfile, 'w') as charmm_inp:
        # will use relative paths because of length of line issues in charmm
        charmm_inp.write('! CHARMM Energy for %s\n' % os.path.relpath(inpfile))
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


def charmm_energies(inpfile, crm_path):
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

    # find where the energy information starts
    i = 0
    for line in all_lines:
        if line[0:9] == 'ENER ENR:':
            startline = i
            break
        i+=1

    energy_types = []
    energy_values = []

    for line in all_lines[startline:]:
        if line[0:4] == 'ENER':
            if '>' in line:
                # unforunately, we can' just split; we have to do it by fixed width
                startcol = 15
                colwidth = 13
                klim = (len(line)-startcol)/colwidth
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
    # remove zero terms from the comparison
    remove_zero_terms(e_out)

    # remove components that are not energy
    nonenergykeys = ['GRMS','VIRI']
    for key in nonenergykeys:
        if key in e_out:
            e_out.pop(key)

    # rename energy terms to standardize
    for key in e_out:
        key = standardize_key(e_out, key)

    # sum up terms to components we can jointly report
    # this will likely need to change because the precise terms reported by CHARMM will vary.
    vanderwaals = ['LJ (SR)', 'IMNBvdw', 'Disper. corr.']
    electrostatic = ['ELEC', 'IMELec', 'EWKSum', 'EWSElf', 'EWEXcl']
    dihedrals = ['Proper Dih.', 'Improper Dih.']
    bonded = ['Bond', 'Angle', 'All dihedrals']
    nonbonded = ['Electrostatic', 'van der Waals'] # must come last, since is a sum of summed terms
    sumterms = [vanderwaals, electrostatic, dihedrals, bonded, nonbonded]
    newkeys = ['van der Waals','Electrostatic', 'All dihedrals', 'Bonded', 'Nonbonded']
    for k, key in enumerate(newkeys):
        e_out[key] =  0 * units.kilocalories_per_mole
        for group in sumterms[k]:
            if group in e_out:
                e_out[key] += e_out[group]
    return e_out, mdout
