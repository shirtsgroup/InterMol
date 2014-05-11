import glob
import numpy
import convert
import os
import pdb


# This script runs a unit test on subdirectories found within Inputs/*/UnitTest/

DES_IN = './Inputs/Desmond/UnitTest'
GRO_IN = './Inputs/Gromacs/UnitTest' 
LMP_IN = './Inputs/Lammps/UnitTest'
OUTPUT_DIR = 'UnitTestOutput'

def parse_args():
    parser = argparse.ArgumentParser(description = 'Run InterMol unit test')

    group_out = parser.add_argument_group('Run unit test on the following input format(s)')
    group_out.add_argument('--desmond', action='store_true',
            help='test conversion of DESMOND input files (Inputs/Desmond/UnitTest/)')
    group_out.add_argument('--gromacs', action='store_true',
            help='test conversion of GROMACS input files (Inputs/Gromacs/UnitTest/)')
    group_out.add_argument('--lammps', action='store_true',
            help='test conversion of LAMMPS input files (Inputs/Lammps/UnitTest/)')

    group_misc = parser.add_argument_group('Other optional arguments')
    group_misc.add_argument('-e', '--energy', dest='energy', action='store_true',
            help='evaluate energy of input and output files for comparison')
    group_misc.add_argument('-d', '--despath', dest='despath',
            metavar='path', default='',
            help='path for DESMOND binary, needed for energy evaluation')
    group_misc.add_argument('-g', '--gropath', dest='gropath',
            metavar='path', default='',
            help='path for GROMACS binary, needed for energy evaluation')
    group_misc.add_argument('-l', '--lmppath', dest='lmppath',
            metavar='path', default='lmp_openmpi',
            help='path for LAMMPS binary, needed for energy evaluation')
    return parser.parse_args()

def test_desmond(args):
    # for now only attempts desmond to desmond conversion of input files
    # to do: need to add checks that right input .cms files are globbed
    files = glob.glob('%s/*/*.cms' % DES_IN) # return list of files that match the string
    results = []

    for f in files:
        try:
            args = ['--des_in', f, '--desmond', '--odir', OUTPUT_DIR, '-e']
            print 'converting %s to desmond file format with command' % f
            print 'python convert.py', ' '.join(args)
            rms = convert.main(args) # reuses code from convert.py
            results.append(rms) 
        except:
            results.append(-1)
    return results

def test_gromacs(args):
    return 1

def test_lammps(args):
    return 1

def main():
    args = parse_args()

    if not os.path.isdir(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)

    if args.desmond:
        results_des = test_desmond(args)
    if args.gromacs:
        results_gro = test_gromacs(args)
    if args.lammps:
        results_lmp = test_lampps(args)
    
    print results_des
    print results_gro
    print results_lmp

    return

    input_files = [glob.glob('%s/*' % x) for x in INPUT_DIR] 

    des_in = input_files[0]

    gro_in = input_files[1]

    lmp_in = input_files[2]



if __name__ == '__main__':
    main()




