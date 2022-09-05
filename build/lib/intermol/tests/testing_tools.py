from collections import OrderedDict
from copy import copy
from glob import glob
import logging
import os
from pkg_resources import resource_filename

import numpy as np

from intermol.exceptions import MultipleValidationErrors
from intermol.utils import add_handler, remove_handler, command_line_flags


ENGINES = ['gromacs', 'lammps', 'desmond', 'amber']

set_suffix = {'gromacs': '.mdp',
              'lammps': None,
              'desmond': '.cfg',
              'amber': '.in'
              }

set_prefix = {'gromacs': 'grompp',
              'lammps': None,
              'desmond': 'onepoint',
              'amber': 'min'
            }




# Log filenames which will be written for each system tested.
INFO_LOG = 'info.log'
DEBUG_LOG = 'debug.log'

# Get reference to various loggers.
logger = logging.getLogger('InterMolLog')   # From convert.py
warning_logger = logging.getLogger('py.warnings')  # From convert.py
testing_logger = logging.getLogger('testing')  # From test_all.py


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


def convert_one_to_all(input_engine, test_type, energy, output_dir,
                       test_tolerance=1e-4):
    """Convert all tests of a type from one engine to all others. """
    flags = {'energy': energy,
             input_engine: True}

    testing_logger.info('Running {} {} tests'.format(input_engine.upper(), test_type))
    output_dir = os.path.join(output_dir, '{}_test_outputs'.format(test_type))

    try:
        os.makedirs(output_dir)
    except OSError:
        if not os.path.isdir(output_dir):
            raise

    results = _convert_from_engine(input_engine, flags, test_type=test_type,
                                   output_dir=output_dir)
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


def _convert_from_engine(input_engine, flags, output_dir, test_type='unit'):
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
    test_dir = resource_filename('intermol', 'tests/{}/{}_tests'.format(
        input_engine, test_type))

    get_test_files = test_finders[input_engine]

    test_files, names = get_test_files(test_dir)
    test_files = list(test_files)
    if len(test_files) < 1:
        testing_logger.info('No {} tests found for {}.'.format(
            test_type, input_engine.upper()))
        return
    # The results of all conversions are stored in nested dictionaries:
    # results = {'gromacs': {'bond1: result, 'bond2: result...},
    #            'lammps': {'bond1: result, 'bond2: result...},
    #            ...}
    per_file_results = OrderedDict((k, None) for k in names)
    results = OrderedDict((engine, copy(per_file_results)) for engine in ENGINES)

    base_output_dir = os.path.join(output_dir, 'from_{}'.format(input_engine))
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

        test_input_dir = os.path.abspath('')

        # is it a vacuum file or not.
        s = ''
        vacstring = '_vacuum'
        if vacstring in test_file[0]:
            s = vacstring

        if input_engine == 'gromacs':
            gro, top = test_file
            flags['gro_in'] = [gro, top]
            flags['inefile'] = os.path.join(test_input_dir, 'gromacs' , 'grompp' + s + '.mdp')

        elif input_engine == 'lammps':
            flags['lmp_in'] = test_file

        elif input_engine == 'desmond':
            flags['des_in'] = test_file
            flags['inefile'] = os.path.join(test_input_dir, 'desmond', 'onepoint' + s + '.cfg')

        elif input_engine == 'amber':
            prmtop, rst = test_file
            flags['amb_in'] = [prmtop, rst]
            flags['inefile'] = os.path.join(test_input_dir, 'amber', 'min' + s + '.in')

        flags['odir'] = odir
        for engine in ENGINES:
            flags[engine] = True
            if set_prefix[engine]:  # set the input file to use for the output
                flags[engine+'_set'] = os.path.join(test_input_dir, engine, set_prefix[engine] + s + set_suffix[engine])

        if flags['lammps']:
            if '_vacuum' in test_file[0]:
                warning_logger.info("I'm seeing an input file (%s) with 'vacuum' in the title; I'm assuming a nonperiodic system with no cutoffs for the LAMMPS input file." % (test_file[0]))
                flags['lmp_settings'] = "pair_style lj/cut/coul/cut 20.0 20.0\nkspace_style none\n\n"
            elif 'lj3_bulk' in test_file[0]:  # annoying workaround since kspace won't work with zero charges.
                flags['lmp_settings'] = "pair_style lj/cut 9.0 \npair_modify tail yes\nkspace_style none\n\n"
            else:
                flags['lmp_settings'] = "pair_style lj/cut/coul/long 15.0 15.0\npair_modify tail yes\nkspace_style pppm 1e-8\n\n"  # Ewald defaults
                #flags['lmp_settings'] = "pair_style lj/cut/coul/long 9.0 9.0\npair_modify tail yes\nkspace_style pppm 1e-6\n\n"  # Ewald defaults
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


def _get_amber_test_files(test_dir):
    prmtop_files = sorted(glob(os.path.join(test_dir, '*/*.prmtop')))
    prmtop_files = [x for x in prmtop_files if not x.endswith('out.prmtop')]
    suffix = ['rst','rst7','crd','inpcrd']
    crd_files = []
    for s in suffix:
        crd_files += glob(os.path.join(test_dir, '*/*.' + s))    # generalize to other coordinate files
    crd_files = sorted(crd_files)    
    names = [os.path.splitext(os.path.basename(prmtop))[0] for prmtop in prmtop_files]
    return zip(prmtop_files, crd_files), names

test_finders = {'gromacs': _get_gromacs_test_files,
                'lammps': _get_lammps_test_files,
                'amber': _get_amber_test_files,
                'desmond': _get_desmond_test_files,
                }
