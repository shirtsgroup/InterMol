InterMol: a conversion tool for molecular dynamics simulations
==============================================================

[![Linux Build Status](https://travis-ci.org/shirtsgroup/InterMol.svg?branch=develop)](https://travis-ci.org/shirtsgroup/InterMol)
[![Coverage Status](https://coveralls.io/repos/shirtsgroup/InterMol/badge.svg?branch=develop)](https://coveralls.io/r/shirtsgroup/InterMol)

We are currently in beta testing phase. Desmond<=>Gromacs<=>Lammps conversions are carried out natively in InterMol.  AMBER->X is carried out by converting AMBER to GROMACS, then to other programs using ParmEd. AMBER->CHARMM is carried out by ParmEd directly.

### [Installation instructions](doc/installation.rst)


### Basic usage
To check out how it works, use the ````convert.py```` script found in the ````intermol```` directory:

```bash
$ python convert.py -h
usage: convert.py [-h] [--des_in file] [--gro_in file file] [--lmp_in file] [--amb_in file file] [-crm_in file] 
                  [--desmond] [--gromacs] [--lammps] [--amber] [--charmm] [--odir directory]
                  [--oname prefix] [-e] [--efile EFILE] 
                  [-dp, --despath] [-ds, --desmondsettings]
                  [-gp, --gropath path] [-gs, --gromacssettings]
                  [-lp, --lmppath] [-ls, --lammpssettings]
                  [-ap, --amberpath] [-as, --ambersettings]
                  [-cp, --charmmpath path] [-cs, --charmmsettings]
                  [-l path] [-f] [-v]

Perform a file conversion

optional arguments:
  -h, --help            show this help message and exit

Choose input conversion format:
  --amb_in file         input files (.prmtop) and (.crd or equivalent) for
                        conversion from AMBER file format
  -crm_in file          input file (.imp) for conversion from CHARMM file format
  --des_in file         .cms file for conversion from DESMOND file format
  --gro_in file file    .gro and .top file for conversion from GROMACS file
                        format
  --lmp_in file         input file for conversion from LAMMPS file format
                        (expects data file in same directory and a read_data call)
  
Choose output conversion format(s):
  --amber               convert to AMBER
  --charmm              convert to CHARMM
  --desmond             convert to DESMOND
  --gromacs             convert to GROMACS
  --lammps              convert to LAMMPS
  
  
Other optional arguments:
  --odir directory      specification of output directory (default: ./)
  --oname prefix        specification of prefix for output filenames (default:
                        inputprefix_converted)
  -e, --energy          evaluate energy of input and output files for
                        comparison
  --efile EFILE         optional run file for energy evaluation (e.g. .cfg,
                        .mdp, .input)
  -dp, --despath path
                        path for DESMOND binary, needed for energy evaluation
  -ds, --desmondsettings settings
                        Desmond .cfg settings file used for energy evaluation
  -gp, --gropath path
                        path for GROMACS binary, needed for energy evaluation
  -gs, --gromacssettings settings
                        Gromacs .mdp settings file used for energy evaluation
  -lp, --lmppath path
                        path for LAMMPS binary, needed for energy evaluation
  -ls, --lammpssettings settings
                        pair_style string to use in the output file. Default
                        is a periodic Ewald simulation
  -ap, --amberpath path
                        path for AMBER binary, needed for energy evaluation
  -as, --ambersettings settings
                        AMBER .in file used for energy evaluation
  -cp, --charmmpath path
                        path for CHARMM binary, needed for energy evaluation
  -cs, --charmmsettings settings
                        CHARMM .in settings used for energy evaluation.
                        Because of the unified input file, settings lines
                        rather than a file are part


  -f, --force           ignore warnings
  -v, --verbose         high verbosity, includes DEBUG level output
  -n, --noncanonical    print the 'noncanonical' enery terms, the ones that
                        only occur in 1 or 2 programs
````

For example, to convert from desmond to gromacs and evalutate the energy of the input and output files:

```bash
mkdir test_output
python convert.py --des_in validation/inputs/Desmond/UnitTest/frag_opls2001/frag_opls2001.cms --gromacs --odir test_output -e
```

Note that the program may not give the same ASCII output; the goal is to
produce the same energy output, and there are frequently multiple ways to
express the same molecular potential energy function using differently
formatted files.

We're also developing a comprehensive test suite to ensure that all InterMol
conversions happen correctly. For instructions on how to run this check out
the [testing README](tests/README.md).

### Citation [![Citing InterMol](https://img.shields.io/badge/DOI-10.1063%2F1.4880555-blue.svg)](http://dx.doi.org/10.1007/s10822-016-9977-1)

InterMol is research software. If you make use of InterMol in scientific
publications please cite the following reference:

```
@article{shirts2016lessons,
  title={Lessons learned from comparing molecular dynamics engines on the SAMPL5 dataset},
  author={Shirts, Michael R and Klein, Christoph and Swails, Jason M and Yin, Jian and Gilson, Michael K and Mobley, David L and Case, David A and Zhong, Ellen D},
  journal={Journal of computer-aided molecular design},
  pages={1--15},
  year={2016},
  doi={doi:10.1007/s10822-016-9977-1}
  publisher={Springer International Publishing}
}
```
