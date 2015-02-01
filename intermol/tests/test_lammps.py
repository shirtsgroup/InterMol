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
    f = logging.Formatter("%(levelname)s %(asctime)s %(message)s", "%Y-%m-%d %H:%M:%S")
    h.setFormatter(f)
    testing_logger.addHandler(h)


def lammps(flags, test_type='unit'):
    """Run all LAMMPS tests of a particular type.

    Args:
        flags (dict): Flags passed to `convert.py`.
        test_type (str, optional): The test suite to run. Defaults to 'unit'.
    Returns:
        results (dict of dicts): The results of all conversions.
            {'gromacs': {'bond1: result, 'bond2: results...},
            'lammps': {'bond1: result, 'bond2: results...},
            ...}

    """
    resource_dir = resource_filename('intermol',
                                     'tests/lammps/{0}_tests'.format(test_type))

    input_files = sorted(glob(os.path.join(resource_dir, '*/*.input')))
    names = [os.path.splitext(os.path.basename(inp))[0] for inp in input_files]

    per_file_results = {k: None for k in names}
    results = {engine: per_file_results for engine in ENGINES}

    unit_test_outputs = '{0}_test_outputs/from_lammps'.format(test_type)
    basedir = os.path.join(os.path.dirname(__file__), unit_test_outputs)
    if not os.path.isdir(basedir):
        os.mkdir(basedir)

    for input_file, name in zip(input_files, names):
        testing_logger.info('Converting {0}'.format(name))
        odir = '{0}/{1}'.format(basedir, name)
        if not os.path.isdir(odir):
            os.mkdir(odir)
        h1, h2 = add_handler(odir)

        flags['lmp_in'] = input_file
        flags['odir'] = odir

        for engine in ENGINES:
            flags[engine] = True

        cmd_line_equivalent = []
        for flag, flag_value in flags.iteritems():
            if not isinstance(flag_value, string_types):
                # E.g. {'lammps': True}
                arg = '--{0}'.format(flag)
            else:
                arg = '--{0} {1}'.format(flag, flag_value)
            cmd_line_equivalent.append(arg)

        logger.info('Converting {0} with command:\n'.format(input_file))
        logger.info('    python convert.py {0}'.format(
            ' '.join(cmd_line_equivalent)))

        try:
            diff = convert.main(flags)
        except Exception as e:
            logger.exception(e)
            for engine in ENGINES:
                results[engine][name] = e
        else:
            for engine, result in diff.iteritems():
                results[engine][name] = result
        remove_handler(h1, h2)

    summarize_results('lammps', results, basedir)
    return results


def test_lammps_unit():
    """Run the LAMMPS unit tests. """
    unit_test_tolerance = 1e-4
    flags = {'unit': True,
             'energy': True,
             'lammps': True}

    testing_logger.info('Running unit tests')

    output_dir = os.path.join(os.path.dirname(__file__), 'unit_test_outputs')
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    results = lammps(flags, test_type='unit')
    zeros = np.zeros(shape=(len(results['lammps'])))
    for engine, tests in results.iteritems():
        tests = np.array(tests.values())
        try:
            passed = np.allclose(tests, zeros, atol=unit_test_tolerance)
        except:
            raise ValueError('Found non-numeric result. This probably means an'
                             'error occured somewhere in the conversion.')
        assert passed, 'Unit tests do not match within {0:.1e} kJ/mol.'.format(
                unit_test_tolerance)
        print('All unit tests match within {0:.1e} kJ/mol.'.format(unit_test_tolerance))


def test_lammps_stress():
    """Run the LAMMPS stress tests. """
    stress_test_tolerance = 1e-4
    flags = {'stress': True,
             'energy': True,
             'lammps': True}

    testing_logger.info('Running stress tests')

    output_dir = os.path.join(os.path.dirname(__file__), 'stress_test_outputs')
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    results = lammps(flags, test_type='stress')
    zeros = np.zeros(shape=(len(results['lammps'])))
    for engine, tests in results.iteritems():
        tests = np.array(tests.values())
        try:
            passed = np.allclose(tests, zeros, atol=stress_test_tolerance)
        except:
            raise ValueError('Found non-numeric result. This probably means an '
                             'error occured somewhere in the conversion.')
        assert passed, 'Stress tests do not match within {0:.1e} kJ/mol.'.format(
                stress_test_tolerance)
        print('All stress tests match within {0:.1e} kJ/mol.'.format(stress_test_tolerance))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Run the LAMMPS tests.')

    type_of_test = parser.add_argument('-t', '--type', metavar='test_type',
            default='unit', help="The type of tests to run: 'unit', 'stress' or 'all'.")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = vars(parser.parse_args())
    if args['type'] in ['unit', 'all']:
        test_lammps_unit()
    if args['type'] in ['stress', 'all']:
        test_lammps_stress()














