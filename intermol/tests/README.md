### Running the tests
In order to cover the large number of functional forms for different energy
terms, we use a testing suite that converts unit tests, e.g. a tiny system with
harmonic bonds between, and stress tests, e.g., a large solvated protein or
micelle, between Gromacs, Lammps and Desmond.

To run ALL tests (this may take a while), use:
```bash
py.test -v
```

To run conversions from one package to all other packages, use `test_all.py`
```bash
python test_all.py -h
usage: test_all.py [-h] [-p engine] [-t test_type] [-e]

Convert tests from one engine to all others.

optional arguments:
  -h, --help            show this help message and exit
  -p engine, --program engine
                        The engine to convert from: gromacs, lammps, desmond
  -t test_type, --type test_type
                        The type of tests to run: unit, stress.
  -e, --energy          Compute and compare the input and output energies.
```

For example, to convert all unit tests from gromacs to all other packages and
compare energies between the inputs and outputs, run:
```bash
python test_all.py --type unit --program gromacs --energy
```