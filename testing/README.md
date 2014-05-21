Use the following two scripts to run InterMol. 

convert.py
==========
````
$ python convert.py -h
usage: convert.py [-h] [--des_in file] [--gro_in file file] [--lmp_in file]
                  [--desmond] [--gromacs] [--lammps] [--odir directory]
                  [--oname prefix] [-e] [--efile EFILE] [-d path] [-g path]
                  [-l path] [-v] [-f]

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
  -v, --verbose         verbosity
  -f, --force           ignore warnings
````

UnitTest.py
===========

````
$ python UnitTest.py -h
usage: PROG [-h] [--desmond] [--gromacs] [--lammps] [-e] [-d path] [-g path]
            [-l path]

         InterMol Unit Testing Script
         --------------------------------
            After specifying input type X, this script will convert files
            found in ./Inputs/X/UnitTest/ to all file formats. All output
            files will be found in ./UnitTestOutput.  
             
         

optional arguments:
  -h, --help            show this help message and exit

Run unit test on the following input format(s):
  --desmond             test conversion of DESMOND files found in
                        Inputs/Desmond/UnitTest/
  --gromacs             test conversion of GROMACS files found in
                        Inputs/Gromacs/UnitTest/
  --lammps              test conversion of LAMMPS files found in
                        Inputs/Lammps/UnitTest/

Other optional arguments:
  -e, --energy          evaluate energy of input and output files for
                        comparison
  -d path, --despath path
                        path for DESMOND binary, needed for energy evaluation
  -g path, --gropath path
                        path for GROMACS binary, needed for energy evaluation
  -l path, --lmppath path
                        path for LAMMPS binary, needed for energy evaluation
````
