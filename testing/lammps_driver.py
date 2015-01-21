import os
import subprocess
import intermol.unit as units
import pdb
import logging
from intermol.lammps_extension.lammps_parser import LammpsParser

logger = logging.getLogger('InterMolLog')

def readFile(infile):
    logger.info('Reading LAMMPS files {0}'.format(infile))
    parser = LammpsParser()
    parser.read_system(infile)
    logger.info('Structure loaded')

def writeFile(outfile):
    logger.info('Writing LAMMPS file {0}'.format(outfile))
    parser = LammpsParser()
    parser.write(outfile)
    logger.info('Write complete')

def lammps_energies(input_file, lmppath='lmp_openmpi'):
    """Evaluate energies of LAMMPS files

    Args:
        input_file = path to input file (expects data file in same folder)
        lmppath = path to LAMMPS binaries
    """
    logger.info('Evaluating energy of {0}'.format(input_file))

    directory, input_file = os.path.split(os.path.abspath(input_file))

    # mdrunin'
    saved_path = os.getcwd()
    os.chdir(directory)

    cmd = "{lmppath} < {input_file}".format(
            lmppath=lmppath, input_file=input_file)
    logger.debug('Running LAMMPS with command:\n    %s' % cmd)
    with open('lammps_stdout.txt', 'w') as out, open('lammps_stderr.txt', 'w') as err:
        exit = subprocess.call(cmd, stdout=out, stderr=err, shell=True)
    os.chdir(saved_path)
    if exit:
        logger.error('Energy evaluation failed. See %s/lammps_stderr.txt' % directory)
        raise Exception('Energy evaluation failed for {0}'.format(input_file))

    # energizin'
    proc = subprocess.Popen(["awk '/E_bond/{getline; print}' %s/lammps_stdout.txt" % (directory)],
            stdout=subprocess.PIPE, shell=True)
    (energies, err) = proc.communicate()
    if not energies:
        raise Exception('Unable to read LAMMPS energy output')


    # give everything units
    data = map(float, energies.split())
    data = [value * units.kilocalories_per_mole for value in data]

    # pack it all up in a dictionary
    types = ['Bond', 'Angle', 'Proper Dih.', 'Improper Dih.', 'Non-bonded',
            'Dispersive', 'Electrostatic', 'Coul. recip.', 'Disper. corr.',
            'Potential']
    e_out = dict(zip(types, data))

    # groupings
    e_out['Electrostatic'] += e_out['Coul. recip.']
    e_out['All dihedrals'] = e_out['Proper Dih.'] + e_out['Improper Dih.']
    return e_out, '%s/lammps_stdout.txt' % directory
