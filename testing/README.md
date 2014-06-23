Use the following two scripts to run InterMol. 

convert.py
==========
````
$ python convert.py -h
usage: convert.py [-h] [--des_in file] [--gro_in file file] [--lmp_in file]
                  [--desmond] [--gromacs] [--lammps] [--odir directory]
                  [--oname prefix] [-e] [--efile EFILE] [-d path] [-g path]
                  [-l path] [-v]

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
  -v, --verbose         high verbosity, includes DEBUG level output
````

systems_test.py
===========

````
$ python systems_test.py -h
usage: PROG [-h] [--unit] [--stress] [-d] [-g] [-l] [-e] [-dp path] [-gp path]
            [-lp path] [-v]

         InterMol Systems Testing Script
         --------------------------------
            After specifying --unit and input type X, this script will 
            convert files found in ./inputs/X/UnitTest/ to all file formats. 
            All output files will be found in ./unit_test_output.
             
            After specifying --stress and input type X, this script will 
            convert files found in ./inputs/X/StressTest/ to all file formats. 
            All output files will be found in ./stress_test_output.

            Additional systems can be tested simply by adding files to the
            appropriate path specified above.
         

optional arguments:
  -h, --help            show this help message and exit

Choose unit tests and/or stress tests:
  --unit                run unit tests found in inputs/****/UnitTest/
  --stress              run stress tests found in inputs/****/StressTest/

Choose one or more of the following input format(s):
  -d, --desmond         test conversion of DESMOND files found in
                        inputs/Desmond/****Test/
  -g, --gromacs         test conversion of GROMACS files found in
                        inputs/Gromacs/****Test/
  -l, --lammps          test conversion of LAMMPS files found in
                        inputs/Lammps/****Test/

Other optional arguments:
  -e, --energy          evaluate energy of input and output files for
                        comparison
  -dp path, --despath path
                        path for DESMOND binary, needed for energy evaluation
  -gp path, --gropath path
                        path for GROMACS binary, needed for energy evaluation
  -lp path, --lmppath path
                        path for LAMMPS binary, needed for energy evaluation
  -v, --verbose         print conversion output to console, -v for INFO level,
                        -vv for DEBUG level
````
