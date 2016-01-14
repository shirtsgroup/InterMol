import argparse
from argparse import RawTextHelpFormatter
import logging
import os
import sys

import intermol.gromacs as gmx
import intermol.lammps as lmp
import intermol.desmond as des
import intermol.amber as amb
import intermol.tests
from intermol.utils import potential_energy_diff, summarize_energy_results


# Make a global logging object.
logger = logging.getLogger('InterMolLog')
logger.setLevel(logging.DEBUG)
logging.captureWarnings(True)
warning_logger = logging.getLogger('py.warnings')


def parse_args(args):
    parser = argparse.ArgumentParser(description='Perform a file conversion',
                                     formatter_class=RawTextHelpFormatter)

    # Input arguments.
    group_in = parser.add_mutually_exclusive_group(required=True)
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
    group_misc.add_argument('--inefile', dest='inefile', default='',
            help='optional run settings file for input energy evaluation (e.g. .cfg, .mdp, .input)')

    # desmond settings 
    group_misc.add_argument('-dp', '--despath', dest='desmond_path',
            metavar='path', default='',
            help='path for DESMOND binary, needed for energy evaluation')
    group_misc.add_argument('-ds', '--desmondsettings', dest='desmond_set',
            metavar='settings', default=None,
            help='Desmond .cfg settings file used for energy evaluation')

    # gromacs settings
    group_misc.add_argument('-gp', '--gropath', dest='gromacs_path',
            metavar='path', default='',
            help='path for GROMACS binary, needed for energy evaluation')
    group_misc.add_argument('-gs', '--gromacssettings', dest='gromacs_set',
            metavar='settings', default=None,
            help='Gromacs .mdp settings file used for energy evaluation')

    # lammps settings
    group_misc.add_argument('-lp', '--lmppath', dest='lammps_path',
            metavar='path', default='',
            help='path for LAMMPS binary, needed for energy evaluation')
    group_misc.add_argument('-ls', '--lammpssettings', dest='lmp_style',
            #metavar='settings', default="pair_style lj/cut/coul/long 15.0 15.0\npair_modify tail yes\nkspace_style pppm 1e-8\n\n",
            metavar='settings', default="pair_style lj/cut/coul/long 9.0 9.0\npair_modify tail yes\nkspace_style pppm 1e-8\n\n",
            help='pair_style string to use in the output file. Default is a periodic Ewald simulation')

    # amber settings
    group_misc.add_argument('-ap', '--amberpath', dest='amber_path',
            metavar='path', default='',
            help='path for AMBER binary, needed for energy evaluation')
    group_misc.add_argument('-as', '--ambersettings', dest='amber_set',
            metavar='settings', default=None,
            help='Amber .in settings file used for energy evaluation')

    # Prints help if no arguments given.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(args)


def main(args=None):
    logger.info('Beginning InterMol conversion')
    if not args:
        args = vars(parse_args(args))

    if args.get('gromacs_path'):
        gmx.GMX_PATH = args['gromacs_path']
    if args.get('lammps_path'):
        lmp.LMP_PATH = args['lammps_path']
    if args.get('desmond_path'):
        des.DES_PATH = args['desmond_path']
    if args.get('amber_path'):
        amb.AMB_PATH = args['amber_path']

    # --------------- PROCESS INPUTS ----------------- #
    output_status = dict()
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
        system = gmx.load(top_in, gro_in)

    elif args.get('des_in'):
        cms_file = args['des_in']
        prefix = os.path.splitext(os.path.basename(cms_file))[0]
        system = des.load(cms_file)

    elif args.get('lmp_in'):
        lammps_file = args['lmp_in']
        prefix = os.path.splitext(os.path.basename(lammps_file))[0]
        system = lmp.load(in_file=lammps_file)

    elif args.get('amb_in'):

        amber_files = args['amb_in']
        prefix = os.path.splitext(os.path.basename(amber_files[0]))[0]

        # Find the prmtop file since order of inputs is not enforced.
        prmtop_in = [x for x in amber_files if x.endswith('.prmtop')]
        assert(len(prmtop_in) == 1)
        prmtop_in = os.path.abspath(prmtop_in[0])

        # Find the crd file since order of inputs is not enforced, not is suffix
        crd_in = [x for x in amber_files if (x.endswith('.rst7') or x.endswith('.crd') or x.endswith('.rst') or x.endswith('.inpcrd'))]
        assert(len(crd_in) == 1)
        crd_in = os.path.abspath(crd_in[0])

        # next, convert to gromacs with ParmEd
        try:
            import parmed
        except ImportError as e:
            output_status['amber'] = e
            raise ImportError('ParmEd is required for AMBER conversions.')

        structure = parmed.amber.AmberParm(prmtop_in, crd_in)
        #Make GROMACS topology
        try:  # TODO: this is just to skirt a pytest issue
            parmed_system = parmed.gromacs.GromacsTopologyFile.from_structure(structure)

            # write out the files.  Should write them out in the proper directory (the one reading in)
            pathprefix = os.path.dirname(prmtop_in)
            fromamber_top_in = os.path.join(pathprefix, prefix + '_from_amber.top')
            fromamber_gro_in = os.path.join(pathprefix, prefix + '_from_amber.gro')
            parmed.gromacs.GromacsTopologyFile.write(parmed_system, fromamber_top_in)
            parmed.gromacs.GromacsGroFile.write(parmed_system, fromamber_gro_in, precision=8)

            # now, read in using gromacs
            system = gmx.load(fromamber_top_in, fromamber_gro_in)
        except OSError as e:
            logger.exception(e)
    else:
        logger.error('No input file')
        sys.exit(1)

    # --------------- WRITE OUTPUTS ----------------- #
    if not args.get('oname'):
        oname = '{0}_converted'.format(prefix)
    else:
        oname = args['oname']
    oname = os.path.abspath(os.path.join(args['odir'], oname))  # Prepend output directory to oname.

    # TODO: factor out exception handling
    if args.get('gromacs'):
        try:
            gmx.save('{0}.top'.format(oname), '{0}.gro'.format(oname), system)
        except Exception as e:
            logger.exception(e)
            output_status['gromacs'] = e
        else:
            output_status['gromacs'] = 'Converted'

    if args.get('lammps'):
        try:
            lmp.save('{0}.input'.format(oname), system, nonbonded_style=args.get('lmp_style'))
        except Exception as e:
            logger.exception(e)
            output_status['lammps'] = e
        else:
            output_status['lammps'] = 'Converted'

    if args.get('desmond'):
        try:
            des.save('{0}.cms'.format(oname), system)
        except Exception as e:
            logger.exception(e)
            output_status['desmond'] = e
        else:
            output_status['desmond'] = 'Converted'

    if args.get('amber'):
        try:
            import parmed
        except ImportError as e:
            output_status['amber'] = e
            raise ImportError('ParmEd is required for AMBER conversions.')

        # first, check if the gro files exit from writing
        gro_out = oname + '.gro'
        top_out = oname + '.top'
        output_status['amber'] = 'Converted'
        if os.path.isfile(gro_out) and os.path.isfile(top_out):
            # if so, use these files.  Load them into ParmEd
            try:
                top = parmed.load_file(top_out, xyz=gro_out)
                try:
                    top.save(oname + '.prmtop', overwrite=True)
                except Exception as e:
                    output_status['amber'] = e
                try:
                    top.save(oname + '.rst7', overwrite=True)
                except Exception as e:
                    output_status['amber'] = e
            except Exception as e:
                output_status['amber'] = e
        else:
            print("Can't convert to AMBER unless gromacs is also selected")

    # --------------- ENERGY EVALUATION ----------------- #
    if args.get('energy'):
        # Run control file paths.
        tests_path = os.path.abspath(os.path.dirname(intermol.tests.__file__))

        # default locations of setting control files

        mdp_in_default = os.path.abspath(os.path.join(tests_path, 'gromacs', 'grompp.mdp'))
        cfg_in_default = os.path.abspath(os.path.join(tests_path, 'desmond', 'onepoint.cfg'))
        in_in_default = os.path.abspath(os.path.join(tests_path, 'amber', 'min.in'))

        # this hardcoding of which input energy files to use should not be here; it should be in the testing files.
        # Evaluate input energies.
        if args.get('gro_in'):
            if args.get('inefile'):
                if os.path.splitext(args.get('inefile'))[-1] != '.mdp':
                    logger.warn("GROMACS energy settings file does not end with .mdp")
                mdp_in = args['inefile']
            else:
                mdp_in = mdp_in_default
            input_type = 'gromacs'
            try:
                e_in, e_infile = gmx.energies(top_in, gro_in, mdp_in, gmx.GMX_PATH)
            except Exception as e:
                output_status['gromacs'] = e
                logger.exception(e)

        elif args.get('lmp_in'):
            if args.get('inefile'):
                logger.warn("LAMMPS energy settings should not require a separate infile")
            input_type = 'lammps'
            try:
                e_in, e_infile = lmp.energies(lammps_file, lmp.LMP_PATH)
            except Exception as e:
                output_status['lammps'] = e
                logger.exception(e)

        elif args.get('des_in'):
            if args.get('inefile'):
                if os.path.splitext(args.get('inefile'))[-1] != '.cfg':
                    logger.warn("DESMOND energy settings file does not end with .cfg")
                cfg_in = args['inefile']
            else:
                cfg_in = cfg_in_default
            input_type = 'desmond'
            try:
                e_in, e_infile = des.energies(cms_file, cfg_in, des.DES_PATH)
            except Exception as e:
                output_status['desmond'] = e
                logger.exception(e)

        elif args.get('amb_in'):
            if args.get('inefile'):
                if os.path.splitext(args.get('inefile'))[-1] != '.in':
                    logger.warn("AMBER energy settings file does not end with .in")
                in_in = args['inefile']
            else:
                in_in = in_in_default
            input_type = 'amber'
            try:
                e_in, e_infile = amb.energies(prmtop_in, crd_in, in_in, amb.AMB_PATH)
            except Exception as e:
                output_status['amber'] = e
                logger.exception(e)
        else:
            logger.warn('No format for input files identified! Code should have never made it here!')

        # Evaluate output energies.
        output_type = []
        e_outfile = []
        e_out = []
        if args.get('gromacs') and output_status['gromacs'] == 'Converted':
            output_type.append('gromacs')
            if args.get('gromacs_set'):
                mdp = args.get('gromacs_set')
            else:
                mdp = mdp_in_default

            try:
                out, outfile = gmx.energies('{0}.top'.format(oname),
                                            '{0}.gro'.format(oname),
                                            mdp, gmx.GMX_PATH)
            except Exception as e:
                output_status['gromacs'] = e
                logger.exception(e)
            else:
                output_status['gromacs'] = potential_energy_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)

        if args.get('lammps') and output_status['lammps'] == 'Converted':
            output_type.append('lammps')

            try:
                out, outfile = lmp.energies('{0}.input'.format(oname), lmp.LMP_PATH)
            except Exception as e:
                output_status['lammps'] = e
                logger.exception(e)
            else:
                output_status['lammps'] = potential_energy_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)

        if args.get('desmond') and output_status['desmond'] == 'Converted':
            output_type.append('desmond')
            if args.get('desmond_set'):
                cfg = args.get('desmond_set')
            else:
                cfg = cfg_in_default

            try:
                out, outfile = des.energies('{0}.cms'.format(oname), cfg, des.DES_PATH)
            except Exception as e:
                output_status['desmond'] = e
                logger.exception(e)
            else:
                output_status['desmond'] = potential_energy_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)

        if args.get('amber') and output_status['amber'] == 'Converted':
            output_type.append('amber')
            if args.get('amber_set'):
                in_amber = args.get('amber_set')
            else:
                in_amber = in_in_default

            try:
                out, outfile = amb.energies('{0}.prmtop'.format(oname),
                                            '{0}.rst7'.format(oname),
                                            in_amber, amb.AMB_PATH)
            except Exception as e:
                output_status['amber'] = e
                logger.exception(e)
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
