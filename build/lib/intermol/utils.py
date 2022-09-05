import logging
import os
from six import string_types
from subprocess import PIPE, Popen

import numpy as np


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
        if np.isnan(data[i][0]):
            line += '%18s' % 'n/a'
        else:
            line += '%18.8f' % (data[i][0])
        for j in range(1, len(data[i])):
            if np.isnan(data[i][j]):
                line += '%18s' % 'n/a'
            else:
                line += '%18.8f' % (data[i][j])
            if np.isnan(data[i][j]) or np.isnan(data[i][0]):
                line += '%18s' % 'n/a'
            else:
                line += '%18.8f' % (data[i][j]-data[i][0])
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


