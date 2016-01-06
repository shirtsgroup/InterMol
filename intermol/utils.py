import logging
import os
from six import string_types
from subprocess import PIPE, Popen


# Log filenames which will be written for each system tested.
INFO_LOG = 'info.log'
DEBUG_LOG = 'debug.log'

# Get reference to various loggers.
logger = logging.getLogger('InterMolLog')   # From convert.py
warning_logger = logging.getLogger('py.warnings')  # From convert.py
testing_logger = logging.getLogger('testing')  # From test_all.py


def command_line_flags(flags):
    """Convert a dict of flags to a string for use on the command line. """
    cmd_line_equivalent = []
    for flag, flag_value in flags.items():
        if isinstance(flag_value, list):
            in_files = ' '.join(flag_value)
            arg = '--{0} {1}'.format(flag, in_files)
        elif not isinstance(flag_value, string_types):
            # E.g. {'gromacs': True}
            arg = '--{0}'.format(flag)
        else:
            arg = '--{0} {1}'.format(flag, flag_value)
        cmd_line_equivalent.append(arg)
    return cmd_line_equivalent


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


def run_subprocess(cmd, engine, stdout_path, stderr_path, stdin=None):
    """Run a subprocess and log the stdout and stderr. """
    logger.debug('Running {0} with command:\n    {1}'.format(engine.upper(), ' '.join(cmd)))

    proc = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input=stdin)
    with open(stdout_path, 'a') as stdout, open(stderr_path, 'a') as stderr:
        stdout.write(out)
        stderr.write(err)
    return proc


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
    h1.close()
    h2.close()

