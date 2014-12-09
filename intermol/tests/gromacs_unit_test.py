import pdb
import sys
import argparse
import glob
import os
import logging
from intermol import convert

# any hardcoded variables
UNIT_DES_IN = './inputs/Desmond/UnitTest'
#UNIT_GRO_IN = './inputs/Gromacs/UnitTest'
UNIT_GRO_IN = './gromacs/unit_tests'
UNIT_LMP_IN = './inputs/Lammps/UnitTest'
UNIT_OUTPUT_DIR = 'unit_test_output'
STRESS_DES_IN = './inputs/Desmond/StressTest'
STRESS_GRO_IN = './inputs/Gromacs/StressTest'
STRESS_LMP_IN = './inputs/Lammps/StressTest'
STRESS_OUTPUT_DIR = 'stress_test_output'
INFO_LOG = 'info.log'
DEBUG_LOG = 'debug.log'
N_FORMATS = 3

# get reference to loggers used in convert.py 
logger = logging.getLogger('InterMolLog')
logger.setLevel(logging.DEBUG)
logging.captureWarnings(True) # redirect warnings module messages to logging system
warning_logger = logging.getLogger('py.warnings') 

# make a logger for systems_test.py
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
     description='''
         InterMol Systems Testing Script
         --------------------------------
            After specifying --unit and input type X, this script will 
            convert files found in ./inputs/X/UnitTest/ to all file formats. 
            All output files will be found in ./unit_test_output.
             
            After specifying --stress and input type X, this script will 
            convert files found in ./inputs/X/StressTest/ to all file formats. 
            All output files will be found in ./stress_test_output.

            Additional systems can be tested simply by adding files to the
            appropriate path specified above.
         ''')

    group_type = parser.add_argument_group('Choose unit tests and/or stress tests')
    group_type.add_argument('--unit', action='store_true',
            help='run unit tests found in inputs/****/UnitTest/')
    group_type.add_argument('--stress', action='store_true',
            help='run stress tests found in inputs/****/StressTest/')

    group_in = parser.add_argument_group('Choose one or more of the following input format(s)')
    group_in.add_argument('-d', '--desmond', action='store_true',
            help='tests conversion of DESMOND files found in inputs/Desmond/****Test/')
    group_in.add_argument('-g', '--gromacs', action='store_true',
            help='tests conversion of GROMACS files found in inputs/Gromacs/****Test/')
    group_in.add_argument('-l', '--lammps', action='store_true',
            help='tests conversion of LAMMPS files found in inputs/Lammps/****Test/')

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
    ''' 
    Adds two FileHandlers to the global logger object
    args:
        dir: path to output directory of the conversion
    return: 
        h1, h2 -- FileHandlers so that they can be removed later

    h1: logs all >=INFO-level messages to dir/INFO_LOG
    h2: logs all >=DEBUG-level messages to dir/DEBUG_LOG
        also includes function and line number for each message 
    '''
    h1 = logging.FileHandler('{dir}/{info}'.format(dir=dir,info=INFO_LOG), mode='w') # don't append to existing file
    f1 = logging.Formatter("%(levelname)s %(asctime)s %(message)s", "%Y-%m-%d %H:%M:%S")
    h1.setFormatter(f1)
    h1.setLevel(logging.INFO)
    logger.addHandler(h1)
    warning_logger.addHandler(h1)

    h2 = logging.FileHandler('{dir}/{debug}'.format(dir=dir,debug=DEBUG_LOG), mode='w')
    f2 = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s", "%Y-%m-%d %H:%M:%S")
    h2.setFormatter(f2)
    h2.setLevel(logging.DEBUG)
    logger.addHandler(h2)
    warning_logger.addHandler(h2)
    return h1, h2

def remove_handler(h1, h2):
    '''Removes the filehandlers h1, h2 so that messages do not continue to be logged to their respective files'''
    logger.removeHandler(h1)
    logger.removeHandler(h2)
    warning_logger.removeHandler(h1)
    warning_logger.removeHandler(h2)

def test_desmond(args, indir, outdir):
    files = sorted(glob.glob('%s/*/*.cms' % indir)) # return list of files that match the string
    files = [x for x in files if not x.endswith('-out.cms')] 
    results = []

    basedir = '%s/FromDesmond' % outdir
    if not os.path.isdir(basedir):
        os.mkdir(basedir)

    for f in files:
        if not args.verbose:
            testing_logger.info('Converting {f}'.format(f=f))
        prefix = f[f.rfind('/') + 1:-4]
        odir = '%s/%s' % (basedir, prefix)
        if not os.path.isdir(odir):
            os.mkdir(odir)
        h1, h2 = add_handler(odir)
        flags = ['--des_in', f, '--desmond', '--gromacs', '--lammps', '--odir', odir]
        flags = add_flags(args, flags)
        logger.info('Converting %s with command:\n    python convert.py %s' % (f,' '.join(flags)))
        try:
            diff = convert.main(flags) # reuses code from convert.py
            assert(len(diff) == N_FORMATS)
            results += diff
        except Exception as e:
            logger.exception(e)
            results += [e]*N_FORMATS
        remove_handler(h1, h2)
    return files, results

def test_gromacs(args, indir, outdir):
    gro_files = sorted(glob.glob('%s/*/*.gro' % indir)) # return list of files that match the string
    gro_files = [x for x in gro_files if not x.endswith('out.gro')] 
    top_files = sorted(glob.glob('%s/*/*.top' % indir)) # return list of files that match the string
    results = []
    
    basedir = '%s/FromGromacs' % outdir
    if not os.path.isdir(basedir):
        os.mkdir(basedir)

    for g, t in zip(gro_files, top_files):
        if not args.verbose:
            testing_logger.info('Converting {g}'.format(g=g))
        prefix = g[g.rfind('/') + 1:-4]
        odir = '%s/%s' % (basedir, prefix)
        if not os.path.isdir(odir):
            os.mkdir(odir)
        h1, h2 = add_handler(odir)
        flags = ['--gro_in', g, t, '--desmond', '--gromacs', '--lammps', '--odir', odir]
        flags = add_flags(args, flags)
        logger.info('Converting %s, %s with command:\n    python convert.py %s' 
                    % (g, t,' '.join(flags)))
        try:
            diff = convert.main(flags) # reuses code from convert.py
            assert(len(diff) == N_FORMATS)
            results += diff
        except Exception as e:
            logger.exception(e)
            results += [e]*N_FORMATS
        remove_handler(h1, h2)
    return gro_files, results

def test_lammps(args, indir, outdir):
    files = sorted(glob.glob('%s/*/*.input' % indir)) # return list of files that match the string
    results = []

    basedir = '%s/FromLammps' % outdir
    if not os.path.isdir(basedir):
        os.mkdir(basedir)

    for f in files:
        if not args.verbose:
            testing_logger.info('Converting {f}'.format(f=f))
        prefix = f[f.rfind('/') + 1:-4]
        odir = '%s/%s' % (basedir, prefix)
        if not os.path.isdir(odir):
            os.mkdir(odir)
        h1, h2 = add_handler(odir)
        flags = ['--lmp_in', f, '--desmond', '--gromacs', '--lammps', '--odir', odir]
        flags = add_flags(args, flags)
        logger.info('Converting %s with command:\n    python convert.py %s' % (f,' '.join(flags)))
        try:
            diff = convert.main(flags) # reuses code from convert.py
            assert(len(diff) == N_FORMATS)
            results += diff
        except Exception as e:
            logger.exception(e)
            results += [e]*N_FORMATS
        remove_handler(h1, h2)
    return files, results

def summarize_results(input_type, files, results, outdir):
    files = [x.split('/')[-1] for x in files] # get filename, not full path
    col1_width = max(len(x) for x in files)
    col1_width = max(col1_width, 10) # 10 is length of Input File
    col2_width = max(len(str(x)) for x in results)
    col2_width = max(col2_width, 28) # 28 is length of Status/Potential Energy Diff
    total_width = col1_width + col2_width + 3
    n = len(files)

    des_res = results[::N_FORMATS]
    print ''
    print '{:^{}}'.format('Results for {0} to DESMOND Conversion'.format(input_type.upper()), total_width)
    print '='*total_width
    print '{:<{}}   {:>{}}'.format('Input File', col1_width, 'Status/Potential Energy Diff', col2_width)
    print '-'*total_width
    for f,r in zip(files, des_res):
        print '{:{}}   {!s:>{}}'.format(f, col1_width, r, col2_width)
    print ''

    gro_res = results[1:][::N_FORMATS]
    print ''
    print '{:^{}}'.format('Results for {0} to GROMACS Conversion'.format(input_type.upper()), total_width)
    print '='*total_width
    print '{:<{}}   {:>{}}'.format('Input File', col1_width, 'Status/Potential Energy Diff', col2_width)
    print '-'*total_width
    for f,r in zip(files, gro_res):
        print '{:{}}   {!s:>{}}'.format(f, col1_width, r, col2_width)
    print ''

    lmp_res = results[2:][::N_FORMATS]
    print ''
    print '{:^{}}'.format('Results for {0} to LAMMPS Conversion'.format(input_type.upper()), total_width)
    print '='*total_width
    print '{:<{}}   {:>{}}'.format('Input File', col1_width, 'Status/Potential Energy Diff', col2_width)
    print '-'*total_width
    for f,r in zip(files, lmp_res):
        print '{:{}}   {!s:>{}}'.format(f, col1_width, r, col2_width)
    print ''
    print 'For the standard output of each conversion, see {dir}/From{type}/[system name]/{log}' \
            .format(dir=outdir, type=input_type, log=INFO_LOG)
    print 'For a detailed DEBUG-level log of each conversion, see {dir}/From{type}/[system name]/{log}' \
            .format(dir=outdir, type=input_type, log=DEBUG_LOG)
    print ''

def main():
    args = parse_args()

    if args.verbose:
        # adds a StreamHandler to convert.py loggers so that
        #  their messages will be printed to console during conversion
        h = logging.StreamHandler()
        if args.verbose == 2:
            h.setLevel(logging.DEBUG) # more verbosity
        else: # args.verbose == 1
            h.setLevel(logging.INFO) # no DEBUG-level messages
        f = logging.Formatter("%(levelname)s %(asctime)s %(message)s", "%Y-%m-%d %H:%M:%S")
        h.setFormatter(f)
        logger.addHandler(h)
        warning_logger.addHandler(h)

    if args.unit:
        testing_logger.info('Running unit tests')
        if not os.path.isdir(UNIT_OUTPUT_DIR):
            os.mkdir(UNIT_OUTPUT_DIR)
        if args.desmond:
            des_input_files, results_des = test_desmond(args, UNIT_DES_IN, UNIT_OUTPUT_DIR)
        if args.gromacs:
            gro_input_files, results_gro = test_gromacs(args, UNIT_GRO_IN, UNIT_OUTPUT_DIR)
        if args.lammps:
            lmp_input_files, results_lmp = test_lammps(args, UNIT_LMP_IN, UNIT_OUTPUT_DIR)
        
        if args.desmond:
            summarize_results('Desmond', des_input_files, results_des, UNIT_OUTPUT_DIR)
        if args.gromacs:
            summarize_results('Gromacs', gro_input_files, results_gro, UNIT_OUTPUT_DIR)
        if args.lammps:
            summarize_results('Lammps', lmp_input_files, results_lmp, UNIT_OUTPUT_DIR)

    if args.stress:
        testing_logger.info('Running stress tests')
        if not os.path.isdir(STRESS_OUTPUT_DIR):
            os.mkdir(STRESS_OUTPUT_DIR)
    
        if args.desmond:
            des_input_files, results_des = test_desmond(args, STRESS_DES_IN, STRESS_OUTPUT_DIR)
        if args.gromacs:
            gro_input_files, results_gro = test_gromacs(args, STRESS_GRO_IN, STRESS_OUTPUT_DIR)
        if args.lammps:
            lmp_input_files, results_lmp = test_lammps(args, STRESS_LMP_IN, STRESS_OUTPUT_DIR)
        
        if args.desmond:
            summarize_results('Desmond', des_input_files, results_des, STRESS_OUTPUT_DIR)
        if args.gromacs:
            summarize_results('Gromacs', gro_input_files, results_gro, STRESS_OUTPUT_DIR)
        if args.lammps:
            summarize_results('Lammps', lmp_input_files, results_lmp, STRESS_OUTPUT_DIR)


if __name__ == '__main__':
    main()
