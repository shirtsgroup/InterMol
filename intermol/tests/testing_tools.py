from collections import OrderedDict
from copy import copy
from glob import glob
import logging
import os
from pkg_resources import resource_filename
from subprocess import Popen, PIPE

import numpy as np
from six import string_types

from intermol.exceptions import MultipleValidationErrors

ENGINES = ['gromacs', 'lammps', 'desmond']

# Log filenames which will be written for each system tested.
INFO_LOG = 'info.log'
DEBUG_LOG = 'debug.log'

# Get reference to various loggers.
logger = logging.getLogger('InterMolLog')   # From convert.py
warning_logger = logging.getLogger('py.warnings')  # From convert.py
testing_logger = logging.getLogger('testing')  # From test_all.py




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


def summarize_results(input_type, results, outdir):
    col1_width = max(len(x) for x in results[input_type])
    # 10 is length of Input File
    col1_width = max(col1_width, 10)
    col2_width = max(len(str(x)) for x in results[input_type].values())
    # 28 is length of Status/Potential Energy Diff.
    col2_width = max(col2_width, 28)
    total_width = col1_width + col2_width + 3

    with open('all_results.txt', 'w') as out:
        for engine in ENGINES:
            result = results[engine]

            out.write('\n')
            out.write('{:^{}}'.format('Results for {0} to {1} Conversion\n'.format(
                input_type.upper(), engine.upper()), total_width))
            out.write('='*total_width + '\n')
            out.write('{:<{}}   {:>{}}\n'.format('Input File', col1_width, 'Status/Potential Energy Diff', col2_width))
            out.write('-'*total_width + '\n')
            for name, res in result.items():
                out.write('{:{}}   {!s:>{}}\n'.format(name, col1_width, res, col2_width))
            out.write('\n')

        out.write('For the standard output of each conversion,'
                  ' see {dir}/[system name]/{log}\n'.format(dir=outdir, log=INFO_LOG))
        out.write('For a detailed DEBUG-level log of each conversion,'
                  ' see {dir}/[system name]/{log}\n'.format(dir=outdir, log=DEBUG_LOG))
        out.write('\n')

    with open('all_results.txt') as out:
        print(out.read())


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


def convert_one_to_all(input_engine, test_type, energy, test_tolerance=1e-4):
    """Convert all tests of a type from one engine to all others. """
    flags = {'energy': energy,
             input_engine: True}

    testing_logger.info('Running {} {} tests'.format(input_engine.upper(), test_type))
    output_dir = resource_filename('intermol', 'tests/{}_test_outputs'.format(test_type))
    try:
        os.makedirs(output_dir)
    except OSError:
        if not os.path.isdir(output_dir):
            raise

    results = _convert_from_engine(input_engine, flags, test_type=test_type)
    if not energy:
        return  # No need to compare energies.

    exceptions = list()
    for output_engine, tests in results.items():
        # All non-numeric results are assumed to be intended error messages such
        # as unconvertible functional forms for specific packages.
        numeric_results = list()
        for val in tests.values():
            try:
                float_val = float(val)
            except (ValueError, TypeError):
                continue
            else:
                numeric_results.append(float_val)
        if not numeric_results:
            exceptions.append('No {0} tests were successfully converted from {1} to {2}'.format(
                test_type, input_engine.upper(), output_engine.upper()))
            continue
        zeros = np.zeros(shape=(len(numeric_results)))
        if np.allclose(numeric_results, zeros, atol=test_tolerance):
            print('All successfully converted {0} tests for {1} to {2} match within {3:.1e} kJ/mol.'.format(
                test_type, input_engine.upper(), output_engine.upper(), test_tolerance))
        else:
            # Continue making all the comparisons before throwing the error.
            exceptions.append('{0} to {1} {2} tests do not match within {3:.1e} kJ/mol.'.format(
                input_engine.upper(), output_engine.upper(), test_type, test_tolerance))

    if exceptions:
        raise MultipleValidationErrors(*exceptions)


def _convert_from_engine(input_engine, flags, test_type='unit'):
    """Run all tests of a particular type from one engine to all others.

    Args:
        engine (str): The input engine.
        flags (dict): Flags passed to `convert.py`.
        test_type (str, optional): The test suite to run. Defaults to 'unit'.
    Returns:
        results (dict of dicts): The results of all conversions.
            {'gromacs': {'bond1: result, 'bond2: results...},
            'lammps': {'bond1: result, 'bond2: results...},
            ...}

    """
    from intermol import convert
    test_dir = resource_filename('intermol', 'tests/{0}/{1}_tests'.format(
        input_engine, test_type))

    print('LOOKING FOR TEST FILES IN {}'.format(test_dir))
    get_test_files = eval('_get_{}_test_files'.format(input_engine))
    test_files, names = get_test_files(test_dir)
    print('FOUND {} TESTFILES FOR ENGINE {}'.format(len(test_files), input_engine))
    # The results of all conversions are stored in nested dictionaries:
    # results = {'gromacs': {'bond1: result, 'bond2: result...},
    #            'lammps': {'bond1: result, 'bond2: result...},
    #            ...}
    per_file_results = OrderedDict((k, None) for k in names)
    results = OrderedDict((engine, copy(per_file_results)) for engine in ENGINES)

    base_output_dir = resource_filename('intermol', 'tests/{0}_test_outputs/from_{1}'.format(test_type, input_engine))
    try:
        os.makedirs(base_output_dir)
    except OSError:
        if not os.path.isdir(base_output_dir):
            raise

    for test_file, name in zip(test_files, names):
        testing_logger.info('Converting {0}'.format(name))
        odir = '{0}/{1}'.format(base_output_dir, name)
        try:
            os.makedirs(odir)
        except OSError:
            if not os.path.isdir(odir):
                raise

        if input_engine == 'gromacs':
            gro, top = test_file
            flags['gro_in'] = [gro, top]
        elif input_engine == 'lammps':
            flags['lmp_in'] = test_file
        elif input_engine == 'desmond':
            flags['des_in'] = test_file

        flags['odir'] = odir
        for engine in ENGINES:
            flags[engine] = True

        h1, h2 = add_handler(odir)
        cmd_line_equivalent = command_line_flags(flags)

        logger.info('Converting {0} with command:\n'.format(test_file))
        logger.info('    python convert.py {0}'.format(
            ' '.join(cmd_line_equivalent)))
        diff = convert.main(flags)
        for engine, result in diff.items():
            results[engine][name] = result
        remove_handler(h1, h2)

    summarize_results(input_engine, results, base_output_dir)
    return results


def _get_gromacs_test_files(test_dir):
    gro_files = sorted(glob(os.path.join(test_dir, '*/*.gro')))
    gro_files = [x for x in gro_files if not x.endswith('out.gro')]
    top_files = sorted(glob(os.path.join(test_dir, '*/*.top')))
    names = [os.path.splitext(os.path.basename(gro))[0] for gro in gro_files]
    return zip(gro_files, top_files), names


def _get_lammps_test_files(test_dir):
    input_files = sorted(glob(os.path.join(test_dir, '*/*.input')))
    names = [os.path.splitext(os.path.basename(inp))[0] for inp in input_files]
    return input_files, names


def _get_desmond_test_files(test_dir):
    cms_files = sorted(glob(os.path.join(test_dir, '*/*.cms')))
    names = [os.path.splitext(os.path.basename(cms))[0] for cms in cms_files]
    return cms_files, names

