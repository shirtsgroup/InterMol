InterMol
========

Conversion tool for molecular simulations. 

We are currently in alpha testing phase, debugging Desmond<=>Gromacs<=>Lammps conversions. 

To check out how it works, use the ````convert.py```` script found in the ````intermol```` directory:

````
$ python convert.py -h
usage: convert.py [-h] [--des_in file] [--gro_in file file] [--lmp_in file]
                  [--desmond] [--gromacs] [--lammps] [--odir directory]
                  [--oname prefix] [-e] [--efile EFILE] [-d path] [-g path]
                  [-l path] [-f] [-v]

Perform a file conversion

optional arguments:
  -h, --help            show this help message and exit

Choose input conversion format:
  --des_in file         .cms file for conversion from DESMOND file format
  --gro_in file file    .gro and .top file for conversion from GROMACS file
                        format
  --lmp_in file         input file for conversion from LAMMPS file format
                        (expects data file in same directory and a read_data
                        call)

Choose output conversion format(s):
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
  -d path, --despath path
                        path for DESMOND binary, needed for energy evaluation
  -g path, --gropath path
                        path for GROMACS binary, needed for energy evaluation
  -l path, --lmppath path
                        path for LAMMPS binary, needed for energy evaluation
  -f, --force           ignore warnings
  -v, --verbose         high verbosity, includes DEBUG level output
````

For example, to convert from desmond to gromacs and evalutate the energy of the input and output files:

````
mkdir test_output
python convert.py --des_in validation/inputs/Desmond/UnitTest/frag_opls2001/frag_opls2001.cms --gromacs --odir test_output -e
````

Note that the program may not give the same ASCII output; the goal is to produce the same energy output, and there are frequently multiple ways to express the same molecular potential energy function using differently formatted files.
