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
from testing_tools import (add_handler, remove_handler, summarize_results,
                           ENGINES)

logger = logging.getLogger('InterMolLog')
testing_logger = logging.getLogger('testing')
if not testing_logger.handlers:
    testing_logger.setLevel(logging.DEBUG)
    h = logging.StreamHandler()
    h.setLevel(logging.INFO)  # ignores DEBUG level for now
    f = logging.Formatter("%(levelname)s %(asctime)s %(message)s",
                          "%Y-%m-%d %H:%M:%S")
    h.setFormatter(f)
    testing_logger.addHandler(h)


def desmond(flags, test_type='unit'):
    """Run all DESMOND tests of a particular type.

    Args:
        flags (dict): Flags passed to `convert.py`.
        test_type (str, optional): The test suite to run. Defaults to 'unit'.
    Returns:
        results (dict of dicts): The results of all conversions.
            {'desmond': {'bond1: result, 'bond2: results...},
            'lammps': {'bond1: result, 'bond2: results...},
            ...}

    """
    resource_dir = resource_filename('intermol',
                                     'tests/desmond/{0}_tests'.format(test_type))

    cms_files = sorted(glob(os.path.join(resource_dir, '*/*.cms')))
    names = [os.path.splitext(os.path.basename(cms))[0] for cms in cms_files]

    # The results of all conversions are stored in nested dictionaries:
    # results = {'desmond': {'bond1: result, 'bond2: results...},
    #            'lammps': {'bond1: result, 'bond2: results...},
    #            ...}
    per_file_results = OrderedDict((k, None) for k in names)
    results = {engine: copy(per_file_results) for engine in ENGINES}

    unit_test_outputs = '{0}_test_outputs/from_desmond'.format(test_type)
    basedir = os.path.join(os.path.dirname(__file__), unit_test_outputs)
    if not os.path.isdir(basedir):
        os.mkdir(basedir)

    for cms, name in zip(cms_files, names):
        testing_logger.info('Converting {0}'.format(name))
        odir = '{0}/{1}'.format(basedir, name)
        if not os.path.isdir(odir):
            os.mkdir(odir)
        h1, h2 = add_handler(odir)

        flags['des_in'] = cms
        flags['odir'] = odir

        for engine in ENGINES:
            flags[engine] = True

        cmd_line_equivalent = []
        for flag, flag_value in flags.iteritems():
            if not isinstance(flag_value, string_types):
                # E.g. {'desmond': True}
                arg = '--{0}'.format(flag)
            else:
                arg = '--{0} {1}'.format(flag, flag_value)
            cmd_line_equivalent.append(arg)

        logger.info('Converting {0} with command:\n'.format(cms))
        logger.info('    python convert.py {0}'.format(
            ' '.join(cmd_line_equivalent)))

        diff = convert.main(flags)
        for engine, result in diff.iteritems():
            results[engine][name] = result
        remove_handler(h1, h2)

    summarize_results('desmond', results, basedir)
    return results


def test_desmond_unit():
    """Run the LAMMPS stress tests. """
    unit_test_tolerance = 1e-4
    flags = {'unit': True,
             'energy': True,
             'desmond': True}

    testing_logger.info('Running unit tests')

    output_dir = os.path.join(os.path.dirname(__file__), 'unit_test_outputs')
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    results = desmond(flags, test_type='unit')
    zeros = np.zeros(shape=(len(results['desmond'])))
    for engine, tests in results.iteritems():
        tests = np.array(tests.values())
        try:
            passed = np.allclose(tests, zeros, atol=unit_test_tolerance)
        except:
            raise ValueError('Found non-numeric result for DESMOND to {0}. This'
                             ' probably means an error occured somewhere in the'
                             ' conversion.'.format(engine.upper()))
        assert passed, 'Unit tests do not match within {0:.1e} kJ/mol.'.format(
                unit_test_tolerance)
        print('All unit tests for DESMOND to {0} match within {1:.1e} kJ/mol.'.format(
                engine.upper(), unit_test_tolerance))


def test_desmond_stress():
    """Run the DESMOND stress tests. """
    stress_test_tolerance = 1e-4
    flags = {'stress': True,
             'energy': True,
             'desmond': True}

    testing_logger.info('Running stress tests')

    output_dir = os.path.join(os.path.dirname(__file__), 'stress_test_outputs')
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    results = desmond(flags, test_type='stress')
    zeros = np.zeros(shape=(len(results['desmond'])))
    for engine, tests in results.iteritems():
        tests = np.array(tests.values())
        try:
            passed = np.allclose(tests, zeros, atol=stress_test_tolerance)
        except:
            raise ValueError('Found non-numeric result. This probably means an '
                             'error occured somewhere in the conversion.')
        assert passed, 'Stress tests do not match within {0:.1e} kJ/mol.'.format(
                stress_test_tolerance)
        print('All unit tests for DESMOND to {1} match within {0:.1e} kJ/mol.'.format(
                stress_test_tolerance, engine.upper()))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Run the DESMOND tests.')

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














