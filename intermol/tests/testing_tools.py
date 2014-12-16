import argparse
import logging
import os
import sys

#ENGINES = ['gromacs', 'lammps', 'desmond']
ENGINES = ['gromacs']

INFO_LOG = 'info.log'
DEBUG_LOG = 'debug.log'

# # get reference to loggers used in convert.py
# logger = logging.getLogger('InterMolLog')
# logger.setLevel(logging.DEBUG)
# logging.captureWarnings(True) # redirect warnings module messages to logging system
# warning_logger = logging.getLogger('py.warnings')
from intermol.convert import warning_logger, logger

# make a logger for OLDtest_all.py
testing_logger = logging.getLogger('systems_test')
testing_logger.setLevel(logging.DEBUG)
h = logging.StreamHandler()
h.setLevel(logging.INFO) # ignores DEBUG level for now
f = logging.Formatter("%(levelname)s %(asctime)s %(message)s", "%Y-%m-%d %H:%M:%S")
h.setFormatter(f)
testing_logger.addHandler(h)


def parse_args():
    parser = argparse.ArgumentParser(prog='PROG',
     formatter_class=argparse.RawDescriptionHelpFormatter,
     description="""
         InterMol Systems Testing Script
         --------------------------------
            After specifying --unit and input type X, this script will
            convert files found in X/unit_tests/ to all file formats.
            All output files will be found in ./unit_test_output.

            After specifying --stress and input type X, this script will
            convert files found in X/stress_tests/ to all file formats.
            All output files will be found in ./stress_test_output.

            Additional systems can be tested simply by adding files to the
            appropriate path specified above.
         """)

    group_type = parser.add_argument_group('Choose unit tests and/or stress tests')
    group_type.add_argument('--unit', action='store_true',
            help='run unit tests found in ****/unit_tests/')
    group_type.add_argument('--stress', action='store_true',
            help='run stress tests found in ****/unit_tests/')

    group_in = parser.add_argument_group('Choose one or more of the following input format(s)')
    group_in.add_argument('-d', '--desmond', action='store_true',
            help='test conversion of DESMOND files found in desmond/****_tests/')
    group_in.add_argument('-g', '--gromacs', action='store_true',
            help='test conversion of GROMACS files found in gromacs/****_tests/')
    group_in.add_argument('-l', '--lammps', action='store_true',
            help='test conversion of LAMMPS files found in lammps/****_tests/')

    group_misc = parser.add_argument_group('Other optional arguments')
    group_misc.add_argument('-e', '--energy', dest='energy', action='store_true',
            help='evaluate energy of input and output files for comparison')
    group_misc.add_argument('-dp', '--despath', dest='despath',
            metavar='path', default='',
            help='path for DESMOND binary, needed for energy evaluation')
    group_misc.add_argument('-gp', '--gropath', dest='gropath',
            metavar='path', default='',
            help='path for GROMACS binary, needed for energy evaluation')
    group_misc.add_argument('-lp', '--lmppath', dest='lmppath',
            metavar='path', default='lmp_openmpi',
            help='path for LAMMPS binary, needed for energy evaluation')
    group_misc.add_argument('-v', '--verbose', dest='verbose', action='count',
            help='print conversion output to console, -v for INFO level, -vv for DEBUG level')

    if len(sys.argv)==1: # if no arguments, print help
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    if not args.unit and not args.stress:
        print 'Please choose a testing type (--unit/--stress)'
        sys.exit(1)
    if not args.desmond and not args.gromacs and not args.lammps:
        print 'Please choose one or more input formats (-d/-g/-l)'
        sys.exit(1)
    return args


def add_flags(args, flags):
    if args.energy:
        flags.append('-e')
    if args.despath:
        flags.append('-dp')
        flags.append(args.despath)
    if args.gropath:
        flags.append('-gp')
        flags.append(args.gropath)
    if args.lmppath:
        flags.append('-lp')
        flags.append(args.lmppath)
    return flags


def add_handler(dir):
    """Adds two FileHandlers to the global logger object.

    Args:
        dir: path to output directory of the conversion
    Returns:
        h1: logs all >=INFO-level messages to dir/INFO_LOG
        h2: logs all >=DEBUG-level messages to dir/DEBUG_LOG
            also includes function and line number for each message
    """
    h1 = logging.FileHandler('{dir}/{info}'.format(dir=dir,info=INFO_LOG), mode='w')
    f1 = logging.Formatter("%(levelname)s %(asctime)s %(message)s",
                           "%Y-%m-%d %H:%M:%S")
    h1.setFormatter(f1)
    h1.setLevel(logging.INFO)
    logger.addHandler(h1)
    warning_logger.addHandler(h1)

    h2 = logging.FileHandler('{dir}/{debug}'.format(dir=dir,debug=DEBUG_LOG), mode='w')
    f2 = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s",
                           "%Y-%m-%d %H:%M:%S")
    h2.setFormatter(f2)
    h2.setLevel(logging.DEBUG)
    logger.addHandler(h2)
    warning_logger.addHandler(h2)
    return h1, h2


def remove_handler(h1, h2):
    """Removes the filehandlers h1, h2.

    Discontiniues messages being logged to their respective files.
    """
    logger.removeHandler(h1)
    logger.removeHandler(h2)
    warning_logger.removeHandler(h1)
    warning_logger.removeHandler(h2)


def summarize_results(input_type, results, outdir):
    col1_width = max(len(x) for x in results[input_type])
    # 10 is length of Input File
    col1_width = max(col1_width, 10)
    col2_width = max(len(str(x)) for x in results[input_type].itervalues())
    # 28 is length of Status/Potential Energy Diff.
    col2_width = max(col2_width, 28)
    total_width = col1_width + col2_width + 3

    with open('all_results.txt', 'w') as out:
        for engine in ENGINES:
            result = results[engine]

            out.write('\n')
            out.write('{:^{}}'.format('Results for {0} to {1} Conversion\n'.format(
                input_type.upper(), engine), total_width))
            out.write('='*total_width + '\n')
            out.write('{:<{}}   {:>{}}\n'.format('Input File', col1_width, 'Status/Potential Energy Diff', col2_width))
            out.write('-'*total_width + '\n')
            for name, res in result.iteritems():
                out.write('{:{}}   {!s:>{}}\n'.format(name, col1_width, res, col2_width))
            out.write('\n')

        out.write('For the standard output of each conversion, see {dir}/[system name]/{log}\n'
                .format(dir=outdir, log=INFO_LOG))
        out.write( 'For a detailed DEBUG-level log of each conversion, see {dir}/[system name]/{log}\n'
                .format(dir=outdir, log=DEBUG_LOG))
        out.write('\n')

    with open('all_results.txt') as out:
        print(out.read())


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

# def initialize_logger(verbose=1):
#     # Adds a StreamHandler to convert.py loggers so that.
#     # Their messages will be printed to console during conversion.
#     h = logging.StreamHandler()
#     if verbose == 2:
#         h.setLevel(logging.DEBUG)  # more verbosity
#     else:
#         h.setLevel(logging.INFO)  # no DEBUG-level messages
#     f = logging.Formatter("%(levelname)s %(asctime)s %(message)s",
#                           "%Y-%m-%d %H:%M:%S")
#     h.setFormatter(f)
#     logger.addHandler(h)
#     warning_logger.addHandler(h)


def main():
    args = vars(parse_args())


if __name__ == '__main__':
    main()