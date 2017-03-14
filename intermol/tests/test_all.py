import logging
import os
import pytest
import sys

from intermol.tests.testing_tools import convert_one_to_all, ENGINES

logger = logging.getLogger('InterMolLog')

# Make a logger for unit and system tests.
testing_logger = logging.getLogger('testing')
if not testing_logger.handlers:
    testing_logger.setLevel(logging.DEBUG)
    h = logging.StreamHandler()
    h.setLevel(logging.INFO)  # ignores DEBUG level for now
    f = logging.Formatter("%(levelname)s %(asctime)s %(message)s",
                          "%Y-%m-%d %H:%M:%S")
    h.setFormatter(f)
    testing_logger.addHandler(h)


def test_gromacs_unit(energy=False, output_dir=os.getcwd()):
    convert_one_to_all(input_engine='gromacs', test_type='unit', energy=energy,
                       output_dir=output_dir)


def test_gromacs_stress(energy=False, output_dir=os.getcwd()):
    convert_one_to_all(input_engine='gromacs', test_type='stress', energy=energy,
                       output_dir=output_dir)


def test_lammps_unit(energy=False, output_dir=os.getcwd()):
    convert_one_to_all(input_engine='lammps', test_type='unit', energy=energy,
                       output_dir=output_dir)


def test_lammps_stress(energy=False, output_dir=os.getcwd()):
    convert_one_to_all(input_engine='lammps', test_type='stress', energy=energy,
                       output_dir=output_dir)


def test_desmond_unit(energy=False, output_dir=os.getcwd()):
    convert_one_to_all(input_engine='desmond', test_type='unit', energy=energy,
                       output_dir=output_dir)


@pytest.mark.skipif(os.getenv('CI') == 'true',
                    reason='Desmond stress tests take too long to run')
def test_desmond_stress(energy=False, output_dir=os.getcwd()):
    convert_one_to_all(input_engine='desmond', test_type='stress', energy=energy,
                       output_dir=output_dir)


def test_amber_unit(energy=False, output_dir=os.getcwd()):
    convert_one_to_all(input_engine='amber', test_type='unit', energy=energy,
                       output_dir=output_dir)


def test_amber_stress(energy=False, output_dir=os.getcwd()):
    convert_one_to_all(input_engine='amber', test_type='stress', energy=energy,
                       output_dir=output_dir)

find_test_func = {('unit', 'gromacs'): test_gromacs_unit,
                  ('unit', 'amber'): test_amber_unit,
                  ('unit', 'desmond'): test_desmond_unit,
                  ('unit', 'lammps'): test_lammps_unit,

                  ('stress', 'gromacs'): test_gromacs_stress,
                  ('stress', 'amber'): test_amber_stress,
                  ('stress', 'desmond'): test_desmond_stress,
                  ('stress', 'lammps'): test_lammps_stress,
}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert tests from one engine to all others.')

    engine = parser.add_argument('-p', '--program', metavar='engine',
            default='gromacs', help="The engine to convert from: {}".format(', '.join(ENGINES)))
    type_of_test = parser.add_argument('-t', '--type', metavar='test_type',
            default='unit', help="The type of tests to run: unit, stress.")
    compute_energies = parser.add_argument('-e', '--energy', dest='compute_energies',
            action='store_true', help="Compute and compare the input and output energies.")
    output_dir = parser.add_argument('-o', '--output', dest='output_dir',
            default=os.getcwd(), help="Output directory.")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = vars(parser.parse_args())
    testing_function = find_test_func[(args['type'].lower(),
                                       args['program'].lower())]
    testing_function(args['compute_energies'], args['output_dir'])
