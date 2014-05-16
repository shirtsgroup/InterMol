import argparse
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

def add_flags(args, flags):
    if args.energy:
        flags.append('-e')
    if args.despath:
        flags.append('-d')
        flags.append(args.despath)
    if args.gropath:
        flags.append('-g')
        flags.append(args.gropath)
    if args.lmppath:
        flags.append('-l')
        flags.append(args.lmppath)
    return flags

def test_desmond(args):
    files = glob.glob('%s/*/*.cms' % DES_IN) # return list of files that match the string
    files = [x for x in files if not x.endswith('-out.cms')] 
    results = []

    odir = '%s/DesmondToDesmond' % OUTPUT_DIR
    if not os.path.isdir(odir):
        os.mkdir(odir)
    for f in files:
        try:
            flags = ['--des_in', f, '--desmond', '--odir', odir]
            flags = add_flags(args, flags)
            print 'converting %s to desmond file format with command' % f
            print 'python convert.py', ' '.join(flags)
            rms = convert.main(flags) # reuses code from convert.py
            results.append(rms) 
        except:
            results.append(-1)

    odir = '%s/DesmondToGromacs' % OUTPUT_DIR
    if not os.path.isdir(odir):
        os.mkdir(odir)
    for f in files:
        try:
            flags = ['--des_in', f, '--gromacs', '--odir', odir]
            flags = add_flags(args, flags)
            print 'converting %s to gromacs file format with command' % f
            print 'python convert.py', ' '.join(flags)
            rms = convert.main(flags) # reuses code from convert.py
            results.append(rms) 
        except:
            results.append(-1)

    odir = '%s/DesmondToLammps' % OUTPUT_DIR
    if not os.path.isdir(odir):
        os.mkdir(odir)
    for f in files:
        try:
            flags = ['--des_in', f, '--lammps', '--odir', odir]
            flags = add_flags(args, flags)
            print 'converting %s to lammps file format with command' % f
            print 'python convert.py', ' '.join(flags)
            rms = convert.main(flags) # reuses code from convert.py
            results.append(rms) 
        except:
            results.append(-1)
    return files, results

def test_gromacs(args):
    gro_files = glob.glob('%s/*/*.gro' % GRO_IN) # return list of files that match the string
    top_files = glob.glob('%s/*/*.top' % GRO_IN) # return list of files that match the string
    results = []

    odir = '%s/GromacsToDesmond' % OUTPUT_DIR
    if not os.path.isdir(odir):
        os.mkdir(odir)
    for g, t in zip(gro_files, top_files):
        try:
            flags = ['--gro_in', g, t, '--desmond', '--odir', odir]
            flags = add_flags(args, flags)
            print 'converting {0}, {1} to desmond file format with command'.format(g,t)
            print 'python convert.py', ' '.join(flags)
            rms = convert.main(flags) # reuses code from convert.py
            results.append(rms) 
        except:
            results.append(-1)

    odir = '%s/GromacsToGromacs' % OUTPUT_DIR
    if not os.path.isdir(odir):
        os.mkdir(odir)
    for g, t in zip(gro_files, top_files):
        try:
            flags = ['--gro_in', g, t, '--gromacs', '--odir', odir]
            flags = add_flags(args, flags)
            print 'converting {0}, {1} to gromacs file format with command'.format(g,t) 
            print 'python convert.py', ' '.join(flags)
            rms = convert.main(flags) # reuses code from convert.py
            results.append(rms) 
        except:
            results.append(-1)

    odir = '%s/DesmondToLammps' % OUTPUT_DIR
    if not os.path.isdir(odir):
        os.mkdir(odir)
    for g, t in zip(gro_files, top_files):
        try:
            flags = ['--gro_in', g, t, '--lammps', '--odir', odir]
            flags = add_flags(args, flags)
            print 'converting {0}, {1} to lammps file format with command'.format(g,t)
            print 'python convert.py', ' '.join(flags)
            rms = convert.main(flags) # reuses code from convert.py
            results.append(rms) 
        except:
            results.append(-1)
    return gro_files, results

def test_lammps(args):
    files = glob.glob('%s/*/*.lmp' % LMP_IN) # return list of files that match the string
    results = []

    odir = '%s/LammpsToDesmond' % OUTPUT_DIR
    if not os.path.isdir(odir):
        os.mkdir(odir)
    for f in files:
        try:
            flags = ['--lmp_in', f, '--desmond', '--odir', odir]
            flags = add_flags(args, flags)
            print 'converting %s to desmond file format with command' % f
            print 'python convert.py', ' '.join(flags)
            rms = convert.main(flags) # reuses code from convert.py
            results.append(rms) 
        except:
            results.append(-1)

    odir = '%s/LammpsToGromacs' % OUTPUT_DIR
    if not os.path.isdir(odir):
        os.mkdir(odir)
    for f in files:
        try:
            flags = ['--lmp_in', f, '--gromacs', '--odir', odir]
            flags = add_flags(args, flags)
            print 'converting %s to gromacs file format with command' % f
            print 'python convert.py', ' '.join(flags)
            rms = convert.main(flags) # reuses code from convert.py
            results.append(rms) 
        except:
            results.append(-1)

    odir = '%s/LammpsToLammps' % OUTPUT_DIR
    if not os.path.isdir(odir):
        os.mkdir(odir)
    for f in files:
        try:
            flags = ['--lmp_in', f, '--lammps', '--odir', odir]
            flags = add_flags(args, flags)
            print 'converting %s to lammps file format with command' % f
            print 'python convert.py', ' '.join(flags)
            rms = convert.main(flags) # reuses code from convert.py
            results.append(rms) 
        except:
            results.append(-1)
    return files, results

def summarize_results(input_type, files, results):
    for i in range(len(results)):
        if results[i] == 0:
            results[i] = 'Converted'
        elif results[i] == 1:
            results[i] = 'Failed at reading input'
        elif results[i] == 2:
            results[i] = 'Failed at writing output'
        elif results[i] == 3:
            results[i] = 'Failed at evaluting energy of input'
        elif results[i] == 4:
            results[i] = 'Failed at evaluating energy of output'

    col1_width = max(len(x) for x in files)
    col2_width = max(len(str(x)) for x in results)
    col2_width = max(col2_width,9) # 9 is length of RMS/Error
    n = len(files)

    des_res = results[0:n]
    print ''
    print '*'*15 + 'Results for {0} to Desmond Conversion'.format(input_type) + '*'*15
    print '{:{}}   {:{}}'.format('File', col1_width, 'RMS/Error', col2_width)
    print '-'*(col1_width+3+col2_width)
    for f,r in zip(files, des_res):
        print '{:{}}   {:>{}}'.format(f, col1_width, r, col2_width)
    print ''

    gro_res = results[n:2*n]
    print ''
    print '*'*15 + 'Results for {0} to Gromacs Conversion'.format(input_type) + '*'*15  
    print '{:{}}   {:{}}'.format('File', col1_width, 'RMS/Error', col2_width)
    print '-'*(col1_width+3+col2_width)
    for f,r in zip(files, gro_res):
        print '{:{}}   {:>{}}'.format(f, col1_width, r, col2_width)
    print ''

    lmp_res = results[2*n::]
    print ''
    print '*'*15 + 'Results for {0} to Lammps Conversion'.format(input_type) + '*'*15
    print '{:{}}   {:{}}'.format('File', col1_width, 'RMS/Error', col2_width)
    print '-'*(col1_width+3+col2_width)
    for f,r in zip(files, lmp_res):
        print '{:{}}   {:>{}}'.format(f, col1_width, r, col2_width)
    print ''

def main():
    args = parse_args()

    if not os.path.isdir(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)

    if args.desmond:
        des_input_files, results_des = test_desmond(args)
    if args.gromacs:
        gro_input_files, results_gro = test_gromacs(args)
    if args.lammps:
        lmp_input_files, results_lmp = test_lammps(args)
    
    if args.desmond:
        summarize_results('Desmond', des_input_files, results_des)
    if args.gromacs:
        summarize_results('Gromacs', gro_input_files, results_gro)
    if args.lammps:
        summarize_results('Lammps', lmp_input_files, results_lmp)

if __name__ == '__main__':
    main()
