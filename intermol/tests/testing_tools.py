import logging
import os

#ENGINES = ['gromacs', 'lammps', 'desmond']
ENGINES = ['gromacs']

# Log filenames which will be written for each system tested.
INFO_LOG = 'info.log'
DEBUG_LOG = 'debug.log'

# Get reference to loggers used in convert.py.
logger = logging.getLogger('InterMolLog')
logger.setLevel(logging.DEBUG)
logging.captureWarnings(True)
warning_logger = logging.getLogger('py.warnings')

# Make a logger for unit and system tests.
testing_logger = logging.getLogger('testing')
if not testing_logger.handlers:
    testing_logger.setLevel(logging.DEBUG)
    h = logging.StreamHandler()
    h.setLevel(logging.INFO)  # ignores DEBUG level for now
    f = logging.Formatter("%(levelname)s %(asctime)s %(message)s", "%Y-%m-%d %H:%M:%S")
    h.setFormatter(f)
    testing_logger.addHandler(h)


def add_handler(directory):
    """Adds two FileHandlers to the global logger object.

    Args:
        directory: path to output directory of the conversion
    Returns:
        h1: logs all >=INFO-level messages to directory/INFO_LOG
        h2: logs all >=DEBUG-level messages to directory/DEBUG_LOG
            also includes function and line number for each message
    """
    h1 = logging.FileHandler('{directory}/{info}'.format(directory=directory, info=INFO_LOG), mode='w')
    f1 = logging.Formatter("%(levelname)s %(asctime)s %(message)s",
                           "%Y-%m-%d %H:%M:%S")
    h1.setFormatter(f1)
    h1.setLevel(logging.INFO)
    logger.addHandler(h1)
    warning_logger.addHandler(h1)

    h2 = logging.FileHandler('{directory}/{debug}'.format(directory=directory, debug=DEBUG_LOG), mode='w')
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
