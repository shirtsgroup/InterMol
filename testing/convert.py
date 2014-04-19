import sys
import os
import argparse
from warnings import warn
import pdb

import numpy as np

import intermol.Driver as Driver
import evaluate
from gromacs_energies import gromacs_energies
from lammps_energies import lammps_energies
from helper_functions import combine_energy_results, print_energy_summary

def parse_args():
    # TODO:
    ############
    #   Priority: remove reliance on extensions, both here and in Driver
    ############
    #   move -h/--help flag to the bottom of the help message or into group_misc
    #   add --version flag

    parser = argparse.ArgumentParser(description = 'Perform a file conversion')

    # input arguments 
    #   (FYI not using argparse's mutually exclusive group since it
    #    doesn't support the description argument => ugly -h output)
    group_in = parser.add_argument_group('Choose input conversion format')
    group_in.add_argument('--des_in', nargs=1, metavar='file',
            help='.cms file for conversion from DESMOND file format')
    # TODO: ORDER OF GROMACS FILES MATTERS
    #       -need to somehow enforce that top gets parsed before gro
    group_in.add_argument('--gro_in', nargs=2, metavar='file',
            help='.gro and .top file for conversion from GROMACS file format')
    # TODO: change functionality so that input file is read and expects
    #       call to appropriate data file within itself
    group_in.add_argument('--lmp_in', nargs=1, metavar='file',
            help='.lmp file for conversion from LAMMPS file format '
                    + '(expects matching .input file)')

    # output arguments
    group_out = parser.add_argument_group('Choose output conversion format(s)')
    group_out.add_argument('--desmond', action='store_true',
            help='convert to DESMOND')
    group_out.add_argument('--gromacs', action='store_true',
            help='convert to GROMACS')
    group_out.add_argument('--lammps', action='store_true',
            help='convert to LAMMPS')

    # other
    group_misc = parser.add_argument_group('Other optional arguments')
    group_misc.add_argument('--odir', metavar='directory', default='.',
            help='specification of output directory (default: ./)')
    group_misc.add_argument('--oname', metavar='prefix', default='',
            help='specification of prefix for output filenames '
                    + '(default: inputprefix_converted)')
    group_misc.add_argument('-e', '--energy', dest='energy', action='store_true',
            help='evaluate energy of input and output files for comparison')
    group_misc.add_argument('--efile', default='',
            help='optional run file for energy evaluation (e.g. .cfg, .mdp, .input)')
    group_misc.add_argument('-d', '--despath', dest='despath',
            metavar='path', default='',
            help='path for DESMOND binary, needed for energy evaluation')
    group_misc.add_argument('-g', '--gropath', dest='gropath',
            metavar='path', default='', 
            help='path for GROMACS binary, needed for energy evaluation')
    group_misc.add_argument('-l', '--lmppath', dest='lmppath',
            metavar='path', default='lmp_openmpi',
            help='path for LAMMPS binary, needed for energy evaluation')
    group_misc.add_argument('-f', '--force', action='store_true',
            help='ignore warnings')
    
    return parser.parse_args()

def main():
    args = parse_args()

    print 'Performing InterMol conversion...'

    # check if inputs files exist
    # NOT checking extension type to allow flexibility in naming
    if args.des_in:
        if not os.path.isfile(args.des_in[0]):
            raise Exception('File not found: {0}'.format(args.des_in[0]))
        files = args.des_in
         # get prefix of file w/o full path
        prefix = args.des_in[0][args.des_in[0].rfind('/') + 1:-4]
    elif args.gro_in:
        for f in args.gro_in: # 2 arguments for top and gro file
            if not os.path.isfile(f):
                raise Exception('File not found: {0}'.format(f))
        files = args.gro_in
        prefix = args.gro_in[0][args.gro_in[0].rfind('/') + 1:-4]
    elif args.lmp_in:
        if not os.path.isfile(args.lmp_in[0]):
            raise Exception('File not found: {0}'.format(args.lmp_in[0]))
        files = args.lmp_in
        prefix = args.lmp_in[0][args.lmp_in[0].rfind('/') + 1:-4]
    else:
        print 'No input file'
        sys.exit(1)

    # check to make sure output directory exists before doing anything
    if not os.path.exists(args.odir):
        raise Exception('Output directory not found: {0}'.format(args.odir))

    # check to make sure a conversion type is selected before doing anything
    if not args.desmond and not args.gromacs and not args.lammps:
        print 'WARNING: no conversion type selected'
        print '         use -f to override'
        if not args.force:
            sys.exit(1)

    # initialize driver
    Driver.initSystem(prefix)  # what does this name do?

    # load files -- Driver should know which type based on extension... 
    # change this eventually?
    Driver.load(*files)

    # output
    if not args.oname: # default is to append _converted to the input prefix
        oname = '%s_converted' % prefix
    else:
        oname = args.oname
    oname = os.path.join(args.odir, oname) # prepend output directory to oname
    if args.desmond:
        print 'Converting to Desmond...writing %s.cms...' % oname
        Driver.write('%s.cms' % oname)
    if args.gromacs:
        print 'Converting to Gromacs...writing %s.gro, %s.top...' % (oname, oname)
        Driver.write('%s.gro' % oname)
        Driver.write('%s.top' % oname)
    if args.lammps:
        print 'Converting to Lammps...writing %s.lmp...' % oname
        Driver.write('%s.lmp' % oname)
    
    # calculate energy of input and output files to compare
    # to do: remove hardcoded parmeter filenames
    #        add lammps
    #        format printing of energy dictionary
    if args.energy:
        # input
        if args.des_in:
            e_in, e_infile = evaluate.desmond_energies(args.des_in[0],
                    'Inputs/Desmond/onepoint.cfg', args.despath)
        elif args.gro_in:
            top = [x for x in args.gro_in if x.endswith('.top')]
            gro = [x for x in args.gro_in if x.endswith('.gro')]
            e_in, e_infile = gromacs_energies(top[0], gro[0],
                    'Inputs/Gromacs/grompp.mdp', args.gropath, '')
        elif args.lmp_in:
            # TODO: fix this when --lmp_in gets changed to read input files
            temp = args.lmp_in[0].split('.')[0] + '.input'
            e_in, e_infile = lammps_energies(temp, args.lmppath) 
        else:
            warn('Something weird went wrong! Code should have never made it here.')
            pass # should never reach here
        
        # output
        if args.desmond:
            e_out, e_outfile = evaluate.desmond_energies('%s.cms' % oname,
                    'Inputs/Desmond/onepoint.cfg', args.despath) 
        if args.gromacs:
            e_out, e_outfile = gromacs_energies('%s.top' % oname,
                    '%s.gro' % oname, 'Inputs/Gromacs/grompp.mdp', args.gropath, '') 
        if args.lammps:
            e_out, e_outfile = lammps_energies('%s.input' % oname,
                    args.lmppath) 
        print 'Input energy file:', e_infile
        print 'Output energy file:', e_outfile
        results = combine_energy_results(e_in, e_out)
        print_energy_summary(results)

if __name__ == '__main__':
    main()
