from glob import glob
import os
from pkg_resources import resource_filename

import logging
import numpy as np
import pytest
from six import string_types

from intermol import convert
from testing_tools import (add_handler, remove_handler, summarize_results,
                           ENGINES)

testing_logger = logging.getLogger('testing')

def gromacs(flags, test_type='unit'):
    """

    Args:
        args:
    Returns:

    """
    resource_dir = resource_filename('intermol',
                                     'tests/gromacs/{0}_tests'.format(test_type))

    gro_files = sorted(glob(os.path.join(resource_dir, '*/*.gro')))
    gro_files = [x for x in gro_files if not x.endswith('out.gro')]
    top_files = sorted(glob(os.path.join(resource_dir, '*/*.top')))
    names = [os.path.splitext(os.path.basename(gro))[0] for gro in gro_files]

    # The results of all conversions are stored in nested dictionaries:
    # results = {'gromacs': {'bond1: result, 'bond2: results...},
    #            'lammps': {'bond1: result, 'bond2: results...},
    #            ...}
    per_file_results = {k: None for k in names}
    results = {engine: per_file_results for engine in ENGINES}

    unit_test_outputs = '{0}_test_outputs/from_gromacs'.format(test_type)
    basedir = os.path.join(os.path.dirname(__file__), unit_test_outputs)
    if not os.path.isdir(basedir):
        os.mkdir(basedir)

    for gro, top, name in zip(gro_files, top_files, names):
        odir = '{0}/{1}'.format(basedir, name)
        if not os.path.isdir(odir):
            os.mkdir(odir)
        h1, h2 = add_handler(odir)

        flags['gro_in'] = [gro, top]
        flags['odir'] = odir

        for engine in ENGINES:
            flags[engine] = True

        cmd_line_equivalent = []
        for k, v in flags.iteritems():
            if isinstance(v, list):
                in_files = ' '.join(v)
                arg = '--{0} {1}'.format(k, in_files)
            elif not isinstance(v, string_types):
                arg = '--{0}'.format(k)
            else:
                arg = '--{0} {1}'.format(k, v)
            cmd_line_equivalent.append(arg)

        testing_logger.info('Converting {0}, {1} with command:\n'.format(gro, top))
        testing_logger.info('    python convert.py {0}'.format(' '.join(cmd_line_equivalent)))

        try:
            diff = convert.main(flags)
        except Exception as e:
            testing_logger.exception(e)
            for engine in ENGINES:
                results[engine][name] = e
        else:
            for engine, result in diff.iteritems():
                results[engine][name] = result
        remove_handler(h1, h2)

    summarize_results('gromacs', results, basedir)
    return results


def test_gromacs_unit():
    """

    Args:
        gromacs:
    Returns:

    """
    flags = {'unit': True,
             'energy': True,
             'gromacs': True}

    testing_logger.info('Running unit tests')

    output_dir = os.path.join(os.path.dirname(__file__), 'unit_test_outputs')
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    results = gromacs(flags, test_type='unit')
    zeros = np.zeros(shape=(len(results['gromacs'])))
    for engine, tests in results.iteritems():
        tests = np.array(tests.values())
        assert np.allclose(tests, zeros, atol=1e-4)


# def test_gromacs_stress():
#     """
#
#     Args:
#         gromacs:
#     Returns:
#
#     """
#     flags = {'stress': True,
#              'energy': True,
#              'gromacs': True}
#
#     testing_logger.info('Running stress tests')
#
#     output_dir = os.path.join(os.path.dirname(__file__), 'stress_test_outputs')
#     if not os.path.isdir(output_dir):
#         os.mkdir(output_dir)
#
#     results = gromacs(flags, test_type='stress')
#     zeros = np.zeros(shape=(len(results['gromacs'])))
#     for engine, tests in results.iteritems():
#         tests = np.array(tests.values())
#         assert np.allclose(tests, zeros, atol=1e-4)

if __name__ == "__main__":
    import pdb
    #test_gromacs_unit()
    #test_gromacs_stress()
    pass














