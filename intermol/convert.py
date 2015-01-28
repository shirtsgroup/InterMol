import argparse
import logging
import os
import sys
import warnings

import numpy as np

from intermol.gromacs import gromacs_driver
from intermol.lammps import lammps_driver
import intermol.tests


# Make a global logging object.
from intermol.tests.testing_tools import which

logger = logging.getLogger('InterMolLog')


def parse_args(args):
    parser = argparse.ArgumentParser(description='Perform a file conversion')

    # Input arguments.
    group_in = parser.add_argument_group('Choose input conversion format')
    group_in.add_argument('--des_in', nargs=1, metavar='file',
            help='.cms file for conversion from DESMOND file format')
    group_in.add_argument('--gro_in', nargs=2, metavar='file',
            help='.gro and .top file for conversion from GROMACS file format')
    group_in.add_argument('--lmp_in', nargs=1, metavar='file',
            help='input file for conversion from LAMMPS file format (expects'
                 ' data file in same directory and a read_data call)')

    # Output arguments.
    group_out = parser.add_argument_group('Choose output conversion format(s)')
    group_out.add_argument('--desmond', action='store_true',
            help='convert to DESMOND')
    group_out.add_argument('--gromacs', action='store_true',
            help='convert to GROMACS')
    group_out.add_argument('--lammps', action='store_true',
            help='convert to LAMMPS')

    # Other options.
    group_misc = parser.add_argument_group('Other optional arguments')
    group_misc.add_argument('--odir', metavar='directory', default='.',
            help='specification of output directory (default: ./)')
    group_misc.add_argument('--oname', metavar='prefix', default='',
            help='specification of prefix for output filenames '
                 '(default: inputprefix_converted)')
    group_misc.add_argument('-e', '--energy', dest='energy', action='store_true',
            help='evaluate energy of input and output files for comparison')
    group_misc.add_argument('--efile', default='',
            help='optional run file for energy evaluation (e.g. .cfg, .mdp, .input)')
    group_misc.add_argument('-dp', '--despath', dest='despath',
            metavar='path', default='',
            help='path for DESMOND binary, needed for energy evaluation')
    group_misc.add_argument('-gp', '--gropath', dest='gropath',
            metavar='path', default='',
            help='path for GROMACS binary, needed for energy evaluation')
    group_misc.add_argument('-lp', '--lmppath', dest='lmppath',
            metavar='path', default='lmp_openmpi',
            help='path for LAMMPS binary, needed for energy evaluation')
    group_misc.add_argument('-f', '--force', dest='force', action='store_true',
            help='ignore warnings (NOTE: may lead to partially correct topologies)')
    group_misc.add_argument('-v', '--verbose', dest='verbose', action='store_true',
            help='high verbosity, includes DEBUG level output')

    # Prints help if no arguments given.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(args)


def main(args=None):
    logger.info('Beginning InterMol conversion')
    if not args:
        args = vars(parse_args(args))

    if args.get('verbose'):
        h.setLevel(logging.DEBUG)
    # Print warnings.
    warnings.simplefilter("always")
    if not args.get('force'):
        # Warnings will be treated as exceptions unless force flag is used.
        warnings.simplefilter("error")

    gropath = args.get('gropath')
    if not gropath:
        gropath = ''
    lmppath = args.get('lmppath')
    if not lmppath:
        for exe in ['lmp_mpi', 'lmp_openmpi']:
            if which(exe):
                lmppath = exe
                break
        else:
            logger.exception('Found no LAMMPS executable.')

    # --------------- PROCESS INPUTS ----------------- #
    if args.get('gro_in'):
        gromacs_files = args['gro_in']

        prefix = os.path.splitext(os.path.basename(gromacs_files[0]))[0]
        # Find the top file since order of inputs is not enforced.
        top_in = [x for x in gromacs_files if x.endswith('.top')]
        assert(len(top_in) == 1)
        top_in = top_in[0]

        # Find the gro file since order of inputs is not enforced.
        gro_in = [x for x in gromacs_files if x.endswith('.gro')]
        assert(len(gro_in) == 1)
        gro_in = gro_in[0]

        system = gromacs_driver.read_file(top_in, gro_in, gropath)
    elif args.get('lmp_in'):
        lammps_file = args['lmp_in']
        prefix = os.path.splitext(os.path.basename(lammps_file))[0]
        system = lammps_driver.read_file(in_file=lammps_file)
    else:
        logger.error('No input file')
        sys.exit(1)

    # --------------- WRITE OUTPUTS ----------------- #
    if not args.get('oname'):
        oname = '{0}_converted'.format(prefix)
    else:
        oname = args['oname']
    oname = os.path.join(args['odir'], oname)  # Prepend output directory to oname.

    output_status = dict()
    # TODO: factor out exception handling
    if args.get('gromacs'):
        try:
            gromacs_driver.write_file(system,
                '{0}.top'.format(oname), '{0}.gro'.format(oname))
        except Exception as e:
            logger.exception(e)
            output_status['gromacs'] = e
        else:
            output_status['gromacs'] = 0

    if args.get('lammps'):
        try:
            lammps_driver.write_file('{0}.input'.format(oname), system)
        except Exception as e:
            logger.exception(e)
            output_status['lammps'] = e
        else:
            output_status['lammps'] = 0

    # --------------- ENERGY EVALUATION ----------------- #
    if args.get('energy'):
        # Run control file paths.
        tests_path = os.path.dirname(intermol.tests.__file__)

        # Required for all gromacs conversions.
        mdp_path = os.path.join(tests_path, 'gromacs', 'grompp.mdp')

        # Evaluate input energies.
        if args.get('gro_in'):
            input_type = 'gromacs'
            tests_path = os.path.dirname(intermol.tests.__file__)
            mdp_path = os.path.join(tests_path, 'gromacs', 'grompp.mdp')
            e_in, e_infile = gromacs_driver.gromacs_energies(top_in, gro_in,
                    mdp_path, gropath, '')
        elif args.get('lmp_in'):
            input_type = 'lammps'
            e_in, e_infile = lammps_driver.lammps_energies(lammps_file,
                                                           lmppath=lmppath)
        else:
            logger.warn('Code should have never made it here!')

        # Evaluate output energies.
        output_type = []
        e_outfile = []
        e_out = []
        if args.get('gromacs') and output_status['gromacs'] == 0:
            output_type.append('gromacs')
            try:
                out, outfile = gromacs_driver.gromacs_energies(
                        '{0}.top'.format(oname), '{0}.gro'.format(oname),
                        mdp_path, gropath, '')
            except Exception as e:
                output_status['gromacs'] = e
                logger.exception(e)
                e_out.append(-1)
                e_outfile.append(-1)
            else:
                output_status['gromacs'] = get_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)

        if args.get('lammps') and output_status['lammps'] == 0:
            output_type.append('lammps')
            try:
                out, outfile = lammps_driver.lammps_energies(
                        '{0}.input'.format(oname), lmppath)
            except Exception as e:
                output_status['lammps'] = e
                logger.exception(e)
                e_out.append(-1)
                e_outfile.append(-1)
            else:
                output_status['lammps'] = get_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)

        # Display energy comparison results.
        out = ['InterMol Conversion Energy Comparison Results','']
        out.append('{0} input energy file: {1}'.format(input_type, e_infile))
        for out_type, file in zip(output_type, e_outfile):
            out.append('{0} output energy file: {1}'.format(out_type, file))
        out += summarize_energy_results(e_in, e_out, input_type, output_type)
        logger.info('\n'.join(out))
    logger.info('Finished!')
    return output_status


def get_diff(e_in, e_out):
    """Returns difference in potential energy.

    arguments:
        e_in  - dictionary of energy groups from input file
        e_out - dictionary of energy groups from output file

    returns:
        potential energy difference in units of the input
    """
    type = 'Potential'
    input = e_in[type]
    diff = e_out[type].in_units_of(input.unit) - input
    return diff._value


def find_match(key, dictionary, unit):
    """Helper function for `summarize_energy_results`. """
    if key in dictionary:
        return dictionary[key].value_in_unit(unit)
    else:
        return np.nan


def summarize_energy_results(energy_input, energy_outputs, input_type, output_types):
    """Creates a table comparing input and output energy groups.

    Args:
        energy_input (dict): energy groups from input file
        energy_output(list): containing dictionary of energy groups or -1 for
            each output file
        input_type (str): input engine
        output_types (list): containing output formats

    Returns:
        out (list of strings): which forms a summary table using "\n".join(out)
    """
    out = []
    # remove failed evaluations (-1 in energy_outputs)
    failed_i = [i for i, x in enumerate(energy_outputs) if x == -1]
    failed = [output_types[i] for i in failed_i]
    output_types = [x for i, x in enumerate(output_types) if i not in failed_i]
    energy_outputs = [x for x in energy_outputs if x != -1]

    # find all distinct labels
    labels = set(energy_input.keys())
    for e_out in energy_outputs:
        for key, value in e_out.iteritems():
            labels.add(key)

    # set up energy comparison table
    labels = list(labels)
    unit = energy_input[energy_input.keys()[0]].unit
    energy_all = [energy_input] + energy_outputs
    data = np.empty((len(labels), len(energy_all)))
    for i in range(len(data)):
        for j in range(len(energy_all)):
            data[i, j] = find_match(labels[i], energy_all[j], unit)

    # TODO: sort table
    out.append('')
    out.append('Energy group summary')
    out.append('=======================================================================')
    header = '%20s %18s ' % ('type', 'input (%s)' % input_type)
    for otype in output_types:
        header += '%37s' %('output (%s) diff (%s)' % (otype, otype))
    out.append(header)
    for i in range(len(data)):
        line = '%20s ' % labels[i]
        line += '%18.8f ' % data[i][0]
        for j in range(1, len(data[i])):
            line += '%18.8f %18.8f' % (data[i][j], data[i][j]-data[i][0])
        out.append(line)
    out.append('')
    # get differences in potential energy
    i = labels.index('Potential')
    diff = data[i, 1::] - data[i, 0]
    for d, otype in zip(diff, output_types):
        out.append('difference in potential energy from %s=>%s conversion: %18.8f'
                    % (input_type, otype, d))
    for f in failed:
        out.append('energy comparison for {0} output failed'.format(f))
    out.append('=======================================================================')
    return out


if __name__ == '__main__':
    # Specifies lowest severity log messages to handle.
    logger.setLevel(logging.DEBUG)
    h = logging.StreamHandler()
    h.setLevel(logging.INFO)  # Ignores DEBUG level for now.
    f = logging.Formatter("%(levelname)s %(asctime)s %(message)s",
                          "%Y-%m-%d %H:%M:%S")
    h.setFormatter(f)
    logger.addHandler(h)

    # Redirect warnings module messages to logging system.
    logging.captureWarnings(True)
    warning_logger = logging.getLogger('py.warnings')
    warning_logger.addHandler(h)

    main()
