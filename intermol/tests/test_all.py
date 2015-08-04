import logging
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


def test_gromacs_unit(energy=True):
    convert_one_to_all(input_engine='gromacs', test_type='unit', energy=energy)


def test_gromacs_stress(energy=True):
    convert_one_to_all(input_engine='gromacs', test_type='stress', energy=energy)


def test_lammps_unit(energy=True):
    convert_one_to_all(input_engine='lammps', test_type='unit', energy=energy)


def test_lammps_stress(energy=True):
    convert_one_to_all(input_engine='lammps', test_type='stress', energy=energy)


def test_desmond_unit(energy=True):
    convert_one_to_all(input_engine='desmond', test_type='unit', energy=energy)


def test_desmond_stress(energy=True):
    convert_one_to_all(input_engine='desmond', test_type='stress', energy=energy)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert tests from one engine to all others.')

    engine = parser.add_argument('-p', '--program', metavar='engine',
            default='gromacs', help="The engine to convert from: {}".format(', '.join(ENGINES)))
    type_of_test = parser.add_argument('-t', '--type', metavar='test_type',
            default='unit', help="The type of tests to run: 'unit' or 'stress'.")
    compute_energies = parser.add_argument('-e', '--energy', dest='compute_energies',
            action='store_true', help="Compute and compare the input and output energies.")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = vars(parser.parse_args())
    func_name = 'test_{}_{}'.format(args['program'].lower(), args['type'])
    testing_function = eval(func_name)
    testing_function(args['compute_energies'])
