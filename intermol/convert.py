import argparse
import logging
import os
import sys
import warnings

import numpy as np

from intermol.gromacs import gromacs_driver
from intermol.lammps import lammps_driver
from intermol.desmond import desmond_driver
from intermol.amber import amber_driver
import intermol.tests
from intermol.tests.testing_tools import which


# Make a global logging object.
logger = logging.getLogger('InterMolLog')
logger.setLevel(logging.DEBUG)
logging.captureWarnings(True)
warning_logger = logging.getLogger('py.warnings')

def record_exception(runtype, logger, e_out, e_outfile, e):
    output_status[runtype] = e
    logger.exception(e)
    e_out.append(-1)
    e_outfile.append(-1)

def parse_args(args):
    parser = argparse.ArgumentParser(description='Perform a file conversion')

    # Input arguments.
    group_in = parser.add_argument_group('Choose input conversion format')
    group_in.add_argument('--des_in', metavar='file',
            help='.cms file for conversion from DESMOND file format')
    group_in.add_argument('--gro_in', nargs=2, metavar='file',
            help='.gro and .top file for conversion from GROMACS file format')
    group_in.add_argument('--lmp_in', metavar='file',
            help='input file for conversion from LAMMPS file format (expects'
                 ' data file in same directory and a read_data call)')
    group_in.add_argument('--amb_in', nargs=2, metavar='file',
            help='input files (.prmtop) and (.crd or equivalent) for conversion from AMBER file format')

    # Output arguments.
    group_out = parser.add_argument_group('Choose output conversion format(s)')
    group_out.add_argument('--desmond', action='store_true',
            help='convert to DESMOND')
    group_out.add_argument('--gromacs', action='store_true',
            help='convert to GROMACS')
    group_out.add_argument('--lammps', action='store_true',
            help='convert to LAMMPS')
    group_out.add_argument('--amber', action='store_true',
            help='convert to AMBER')

    # Other options.
    group_misc = parser.add_argument_group('Other optional arguments')
    group_misc.add_argument('--odir', metavar='directory', default='.',
            help='specification of output directory (default: ./)')
    group_misc.add_argument('--oname', metavar='prefix', default='',
            help='specification of prefix for output filenames '
                 '(default: inputprefix_converted)')

    group_misc.add_argument('-e', '--energy', dest='energy', action='store_true',
            help='evaluate energy of input and output files for comparison')
    group_misc.add_argument('--inefile', default='',
            help='optional run file for input energy evaluation (e.g. .cfg, .mdp, .input)')

    # desmond settings 
    group_misc.add_argument('-dp', '--despath', dest='des_path',
            metavar='path', default='',
            help='path for DESMOND binary, needed for energy evaluation')
    group_misc.add_argument('-ds', '--desmondsettings', dest='desset',
            metavar='settings', default='',
            help='Desmond .cfg settings file used for energy evaluation')

    # gromacs settings
    group_misc.add_argument('-gp', '--gropath', dest='gro_path',
            metavar='path', default='',
            help='path for GROMACS binary, needed for energy evaluation')
    group_misc.add_argument('-gs', '--gromacssettings', dest='gro_set',
            metavar='settings', default='',
            help='Gromacs .mdp settings file used for energy evaluation')

    # lammps settings
    group_misc.add_argument('-lp', '--lmppath', dest='lmp_path',
            metavar='path', default='',
            help='path for LAMMPS binary, needed for energy evaluation')
    group_misc.add_argument('-ls', '--lammpssettings', dest='lmpset',
            metavar='settings', default='',
            help='LAMMPS .input settings file used for energy evaluation. This is simply the pair_style and long range coulomb type to use')

    # amber settings
    group_misc.add_argument('-ap', '--amberpath', dest='amberpath',
            metavar='path', default='',
            help='path for AMBER binary, needed for energy evaluation')
    group_misc.add_argument('-as', '--ambersettings', dest='ambset',
            metavar='settings', default='',
            help='Amber .in settings file used for energy evaluation')

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

    # Paths to simulator executables.
    # GROMACS
    gro_path = args.get('gro_path')
    if not gro_path:
        gro_path = ''
    # LAMMPS
    lpm_path = args.get('lpm_path')
    if not lpm_path:
        for exe in ['lmp_mpi', 'lmp_serial', 'lmp_openmpi']:
            if which(exe):
                lpm_path = exe
                break
        else:
            logger.exception('Found no LAMMPS executable.')
    # DESMOND
    des_path = args.get('des_path')

    #AMBERPATH
    amberpath = args.get('amberpath')

    # --------------- PROCESS INPUTS ----------------- #
    if args.get('gro_in'):
        gromacs_files = args['gro_in']

        prefix = os.path.splitext(os.path.basename(gromacs_files[0]))[0]
        # Find the top file since order of inputs is not enforced.
        top_in = [x for x in gromacs_files if x.endswith('.top')]
        assert(len(top_in) == 1)
        top_in = os.path.abspath(top_in[0])

        # Find the gro file since order of inputs is not enforced.
        gro_in = [x for x in gromacs_files if x.endswith('.gro')]
        assert(len(gro_in) == 1)
        gro_in = os.path.abspath(gro_in[0])
        system = gromacs_driver.read_file(top_in, gro_in, gro_path)

    elif args.get('des_in'):
        cms_file = args['des_in']
        prefix = os.path.splitext(os.path.basename(cms_file))[0]
        system = desmond_driver.read_file(cms_file=cms_file)

    elif args.get('lmp_in'):
        lammps_file = args['lmp_in']
        prefix = os.path.splitext(os.path.basename(lammps_file))[0]
        system = lammps_driver.read_file(in_file=lammps_file)

    elif args.get('amb_in'):
        amber_files = args['amb_in']
        prefix = os.path.splitext(os.path.basename(amber_files[0]))[0]

        # Find the prmtop file since order of inputs is not enforced.
        prmtop_in = [x for x in amber_files if x.endswith('.prmtop')]
        assert(len(prmtop_in) == 1)
        prmtop_in = os.path.abspath(prmtop_in[0])

        # Find the crd file since order of inputs is not enforced.
        crd_in = [x for x in amber_files if (x.endswith('.rst7') or x.endswith('.crd') or x.endswith('.rst') or x.endswith('.inpcrd'))]
        assert(len(crd_in) == 1)
        crd_in = os.path.abspath(crd_in[0])

        # next, convert to gromacs with ParmEd
        try:
            import parmed
        except Exception as e:
            record_exception('amber',logger, e_out, e_outfile, e)

        structure = parmed.amber.AmberParm(prmtop_in,crd_in)
        #Make GROMACS topology
        parmed_system = parmed.gromacs.GromacsTopologyFile.from_structure(structure)
        # write out the files
        fromamber_top_in = prefix + '_from_amber.top'
        fromamber_gro_in = prefix + '_from_amber.gro'
        parmed.gromacs.GromacsTopologyFile.write(parmed_system, fromamber_top_in)
        parmed.gromacs.GromacsGroFile.write(parmed_system, fromamber_gro_in, precision = 8)

        # now, read in using gromacs
        system = gromacs_driver.read_file(fromamber_top_in, fromamber_gro_in, gro_path)
    else:
        logger.error('No input file')
        sys.exit(1)

    # --------------- WRITE OUTPUTS ----------------- #
    if not args.get('oname'):
        oname = '{0}_converted'.format(prefix)
    else:
        oname = args['oname']
    oname = os.path.abspath(os.path.join(args['odir'], oname))  # Prepend output directory to oname.

    output_status = dict()
    # TODO: factor out exception handling
    if args.get('gromacs'):
        try:
            gromacs_driver.write_file(system, '{0}.top'.format(oname), '{0}.gro'.format(oname))
        except Exception as e:
            logger.exception(e)
            output_status['gromacs'] = e
        else:
            output_status['gromacs'] = 'Converted'

    if args.get('lammps'):
        try:
            lammps_driver.write_file('{0}.input'.format(oname), system)
        except Exception as e:
            logger.exception(e)
            output_status['lammps'] = e
        else:
            output_status['lammps'] = 'Converted'

    if args.get('desmond'):
        try:
            desmond_driver.write_file('{0}.cms'.format(oname), system)
        except Exception as e:
            logger.exception(e)
            output_status['desmond'] = e
        else:
            output_status['desmond'] = 'Converted'

    if args.get('amber'):
        try:
            import parmed
        except Exception as e:
            record_exception('amber', logger, e_out, e_outfile, e)

        # first, check if the gro files exit from writing
        gro_out = oname + '.gro'
        top_out = oname + '.top'
        top = None
        e = None
        if (os.path.isfile(gro_out) and os.path.isfile(top_out)):
            # if so, use these files.  Load them into ParmEd
            try:
                top = parmed.load_file(top_out, xyz=gro_out)
            except Exception as e:
                output_status['amber'] = e
            if top:
                try:
                    parm = parmed.amber.AmberParm.from_structure(top)
                    prmtop_out = oname + '.prmtop'
                    crd_out = oname + '.rst7'
                    try:
                        parm.write_parm(prmtop_out)
                    except Exception as e:
                        output_status['amber'] = e
                    try:
                        parm.write_rst7(crd_out)
                    except Exception as e:
                        output_status['amber'] = e
                    if e == None:
                        output_status['amber'] = 'Converted'
                except Exception as e:
                    output_status['amber'] = e
        else:
            print "Can't convert to AMBER unless gromacs is also selected"


    # --------------- ENERGY EVALUATION ----------------- #
    if args.get('energy'):
        # Run control file paths.
        tests_path = os.path.abspath(os.path.dirname(intermol.tests.__file__))

        mdp_path = os.path.abspath(os.path.join(tests_path, 'gromacs', 'grompp.mdp'))
        cfg_path = os.path.abspath(os.path.join(tests_path, 'desmond', 'onepoint.cfg'))
        lmp_path = os.path.abspath(os.path.join(tests_path, 'lammps', 'lammps.input_include'))
        amb_path = os.path.abspath(os.path.join(tests_path, 'amber', 'amber.in'))

        # Evaluate input energies.
        if args.get('gro_in'):
            if args.get('gro_set'):
                mdp_path = args['gro_set']
            else:
                if gro_in.endswith('_vacuum.gro'):
                    mdp_path = os.path.abspath(os.path.join(tests_path, 'gromacs', 'grompp_vacuum.mdp'))
            input_type = 'gromacs'
            e_in, e_infile = gromacs_driver.gromacs_energies(top_in, gro_in, mdp_path, gro_path, '')
        elif args.get('lmp_in'):
            if args.get('lmp_set'):
                lmp_path = args['lmp_set']
            input_type = 'lammps'
            # add the insertion of pair types
            e_in, e_infile = lammps_driver.lammps_energies(lammps_file, lpm_path=lpm_path)

        elif args.get('des_in'):
            if args.get('des_set'):
                cfg_path = args['des_set']
            else:
                if cms_file.endswith('_vacuum.cms'):
                    cfg_path = os.path.abspath(os.path.join(tests_path, 'desmond', 'onepoint_vacuum.cfg'))
            input_type = 'desmond'
            e_in, e_infile = desmond_driver.desmond_energies(cms_file, cfg_path, des_path=des_path)

        if args.get('amb_in'):
            if args.get('amb_set'):
                ambin_path = args['amb_set']
            else:
                ambin_path = os.path.abspath(os.path.join(tests_path, 'amber', 'min.in'))
            input_type = 'amber'
            e_in = amber_driver.amber_energies(prmtop_in, crd_in, ambin_path, gro_path, '')
        else:
            logger.warn('Code should have never made it here!')

        # Evaluate output energies.
        output_type = []
        e_outfile = []
        e_out = []
        if args.get('gromacs') and output_status['gromacs'] == 'Converted':
            output_type.append('gromacs')
            try:
                out, outfile = gromacs_driver.gromacs_energies(
                    '{0}.top'.format(oname), '{0}.gro'.format(oname), mdp_path, gro_path, '')
            except Exception as e:
                record_exception('gromacs',logger, e_out, e_outfile, e)
            else:
                output_status['gromacs'] = potential_energy_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)

        if args.get('lammps') and output_status['lammps'] == 'Converted':
            output_type.append('lammps')
            try:
                out, outfile = lammps_driver.lammps_energies(
                    '{0}.input'.format(oname), lpm_path=lpm_path)
            except Exception as e:
                record_exception('lammps',logger, e_out, e_outfile, e)
            else:
                output_status['lammps'] = potential_energy_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)

        if args.get('desmond') and output_status['desmond'] == 'Converted':
            output_type.append('desmond')
            try:
                out, outfile = desmond_driver.desmond_energies(
                    '{0}.cms'.format(oname), cfg_path, des_path)
            except Exception as e:
                record_exception('desmond',logger, e_out, e_outfile, e)
            else:
                output_status['desmond'] = potential_energy_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)

        if args.get('amber') and output_status['amber'] == 'Converted':
            output_type.append('amber')
            try:
                out, outfile = amber_driver.amber_energies(
                    '{0}.cms'.format(oname), amb_path, ambpath)
            except Exception as e:
                record_exception('amber',logger, e_out, e_outfile, e)
            else:
                output_status['amber'] = potential_energy_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)

        # Display energy comparison results.
        out = ['InterMol Conversion Energy Comparison Results', '',
               '{0} input energy file: {1}'.format(input_type, e_infile)]
        for out_type, file in zip(output_type, e_outfile):
            out.append('{0} output energy file: {1}'.format(out_type, file))
        out += summarize_energy_results(e_in, e_out, input_type, output_type)
        logger.info('\n'.join(out))
    logger.info('Finished!')
    return output_status


def potential_energy_diff(e_in, e_out):
    """Returns difference in potential energy.

    arguments:
        e_in  - dictionary of energy groups from input file
        e_out - dictionary of energy groups from output file

    returns:
        potential energy difference in units of the input
    """
    energy_type = 'Potential'
    input_energy = e_in[energy_type]
    diff = e_out[energy_type].in_units_of(input_energy.unit) - input_energy
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
    # Remove failed evaluations (-1 in energy_outputs)
    failed_i = [i for i, x in enumerate(energy_outputs) if x == -1]
    failed = [output_types[i] for i in failed_i]
    output_types = [x for i, x in enumerate(output_types) if i not in failed_i]
    energy_outputs = [x for x in energy_outputs if x != -1]

    # Find all distinct labels
    labels = set(energy_input.keys())
    for e_out in energy_outputs:
        for key, value in e_out.items():
            labels.add(key)

    # Set up energy comparison table
    labels = list(labels)
    unit = energy_input[list(energy_input.keys())[0]].unit
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
        header += '%37s' % ('output (%s) diff (%s)' % (otype, otype))
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
    for fail in failed:
        out.append('energy comparison for {0} output failed'.format(fail))
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
