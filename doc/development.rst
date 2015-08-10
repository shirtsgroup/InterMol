============
Contributing
============
Contributions are welcome and they are greatly appreciated! Every little bit
helps and credit will always be given.

Bug Reports
-----------
If you find a bug, please file a detailed issue `here
<https://github.com/shirtsgroup/InterMol/issues>`_ and **provide input files**
so that we can try to reproduce the problem.

Code Style
----------
`PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ describes the style guide
for python code. Please try to follow it when reasonable, in particular the
naming conventions (e.g. public methods should use ``lower_case_with_underscores``
instead of ``camelCase``).

Two helpful tools for checking your code are ``flake8`` and ``pylint``. The latter
is a lot pickier but generally, your code should pass most ``flake8`` tests::

 $ pip install flake8 pylint

Also some IDE's like `PyCharm <https://www.jetbrains.com/pycharm/>`_ will
notify you of PEP8 violations as you code.

.. note:: Some older chunks of the code do not yet adhere to PEP8. If you find
          something easily fixable, please do so!

Running our tests
-----------------

In order to cover the large number of functional forms for different energy
terms, we use a testing suite that converts unit tests, e.g. a tiny system with
harmonic bonds between, and stress tests, e.g., a large solvated protein or
micelle, between Gromacs, Lammps and Desmond.

To run ALL tests (this may take a while), use::

    $ py.test -v

To run conversions from one package to all other packages, use `test_all.py`::

    $ python test_all.py -h
    usage: test_all.py [-h] [-p engine] [-t test_type] [-e]

    Convert tests from one engine to all others.

    optional arguments:
      -h, --help            show this help message and exit
      -p engine, --program engine
                            The engine to convert from: gromacs, lammps, desmond
      -t test_type, --type test_type
                        The type of tests to run: unit, stress.
      -e, --energy          Compute and compare the input and output energies.

For example, to convert all unit tests from gromacs to all other packages and
compare energies between the inputs and outputs, run::

    $ python test_all.py --type unit --program gromacs --energy

