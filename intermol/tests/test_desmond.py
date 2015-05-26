from collections import OrderedDict
from copy import copy
from glob import glob
import os
from pkg_resources import resource_filename

import logging
import numpy as np
from six import string_types
import sys

from intermol import convert
from intermol.tests.testing_tools import (add_handler, remove_handler,
                                          summarize_results, ENGINES)

logger = logging.getLogger('InterMolLog')
testing_logger = logging.getLogger('testing')
if not testing_logger.handlers:
    testing_logger.setLevel(logging.DEBUG)
    h = logging.StreamHandler()
    h.setLevel(logging.INFO)  # ignores DEBUG level for now
    f = logging.Formatter("%(levelname)s %(asctime)s %(message)s", "%Y-%m-%d %H:%M:%S")
    h.setFormatter(f)
    testing_logger.addHandler(h)


def test_desmond_unit():
    """Run the Desmond stress tests. """
    flags = {'unit': True,
             'energy': True,
             'desmond': True}

    testing_logger.info('Running unit tests')
    output_dir = resource_filename('intermol', 'tests/unit_test_outputs')
    try:
        os.makedirs(output_dir)
    except OSError:
        if not os.path.isdir(output_dir):
            raise
    _run_desmond_and_compare(flags, test_tolerance=1e-4, test_type='unit')


def test_desmond_stress():
    """Run the Desmond stress tests. """
    flags = {'stress': True,
             'energy': True,
             'desmond': True}

    testing_logger.info('Running stress tests')
    output_dir = resource_filename('intermol', 'tests/stress_test_outputs')
    try:
        os.makedirs(output_dir)
    except OSError:
        if not os.path.isdir(output_dir):
            raise
    _run_desmond_and_compare(flags, test_tolerance=1e-4, test_type='stress')


def _run_desmond_and_compare(flags, test_tolerance, test_type='unit'):
    results = _convert_from_desmond(flags, test_type=test_type)
    zeros = np.zeros(shape=(len(results['desmond'])))
    for engine, tests in results.items():
        tests = np.array(tests.values())
        try:
            passed = np.allclose(tests, zeros, atol=test_tolerance)
        except TypeError:
            raise Exception('Found non-numeric result for Desmond to {0}. This '
                            'probably means an error occured somewhere in the '
                            'conversion.\n\n'.format(engine.upper()))
        if passed:
            print('All {0} tests for Desmond to {1} match within {2:.1e} kJ/mol.'.format(test_type, engine.upper(), test_tolerance))
        else:
            raise Exception('{0} tests do not match within {1:.1e} kJ/mol.'.format(test_type, test_tolerance))


def _convert_from_desmond(flags, test_type='unit'):
    """Run all Desmond tests of a particular type.

    Args:
        flags (dict): Flags passed to `convert.py`.
        test_type (str, optional): The test suite to run. Defaults to 'unit'.
    Returns:
        results (dict of dicts): The results of all conversions.
            {'gromacs': {'bond1: result, 'bond2: results...},
            'lammps': {'bond1: result, 'bond2: results...},
            ...}

    """
    test_dir = resource_filename('intermol', 'tests/desmond/{0}_tests'.format(test_type))

    gro_files = sorted(glob(os.path.join(test_dir, '*/*.gro')))
    gro_files = [x for x in gro_files if not x.endswith('out.gro')]
    top_files = sorted(glob(os.path.join(test_dir, '*/*.top')))
    names = [os.path.splitext(os.path.basename(gro))[0] for gro in gro_files]

    # The results of all conversions are stored in nested dictionaries:
    # results = {'gromacs': {'bond1: result, 'bond2: result...},
    #            'lammps': {'bond1: result, 'bond2: result...},
    #            ...}
    per_file_results = OrderedDict((k, None) for k in names)
    results = {engine: copy(per_file_results) for engine in ENGINES}

    base_output_dir = resource_filename('intermol', 'tests/{0}_test_outputs/from_desmond'.format(test_type))
    try:
        os.makedirs(base_output_dir)
    except OSError:
        if not os.path.isdir(base_output_dir):
            raise

    for gro, top, name in zip(gro_files, top_files, names):
        testing_logger.info('Converting {0}'.format(name))
        odir = '{0}/{1}'.format(base_output_dir, name)
        try:
            os.makedirs(odir)
        except OSError:
            if not os.path.isdir(odir):
                raise
        h1, h2 = add_handler(odir)

        flags['gro_in'] = [gro, top]
        flags['odir'] = odir
        for engine in ENGINES:
            flags[engine] = True

        cmd_line_equivalent = _command_line_flags(flags)

        logger.info('Converting {0}, {1} with command:\n'.format(gro, top))
        logger.info('    python convert.py {0}'.format(
            ' '.join(cmd_line_equivalent)))

        diff = convert.main(flags)
        for engine, result in diff.items():
            results[engine][name] = result
        remove_handler(h1, h2)

    summarize_results('desmond', results, base_output_dir)
    return results


def _command_line_flags(flags):
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


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Run the Desmond tests.')

    type_of_test = parser.add_argument('-t', '--type', metavar='test_type',
            default='unit', help="The type of tests to run: 'unit', 'stress' or 'all'.")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # TODO: Rewrite assertions or testing calls so that you can run 'all' but
    #       still have working py.test.
    args = vars(parser.parse_args())
    if args['type'] in ['unit', 'all']:
        test_desmond_unit()
    if args['type'] in ['stress', 'all']:
        test_desmond_stress()














