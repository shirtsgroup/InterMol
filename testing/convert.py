import sys
import os
import argparse
import logging
import pdb
import numpy as np

from intermol.system import System
import helper_functions
import desmond_driver
import gromacs_driver
import lammps_driver

# global constants

DES_POS = 0
GRO_POS = 1
LMP_POS = 2
N_FORMATS = 3

# Make a global logging object.
logger = logging.getLogger('InterMolLog')
if __name__ == '__main__':
    logger.setLevel(logging.DEBUG) # specifies lowest severity log messages to handle
    h = logging.StreamHandler()
    h.setLevel(logging.INFO) # ignores DEBUG level for now
    f = logging.Formatter("%(levelname)s %(asctime)s %(message)s", "%Y-%m-%d %H:%M:%S")
    h.setFormatter(f)
    logger.addHandler(h)

def get_parser():
    parser = argparse.ArgumentParser(description = 'Perform a file conversion')

    # input arguments
    group_in = parser.add_argument_group('Choose input conversion format')
    group_in.add_argument('--des_in', nargs=1, metavar='file',
            help='.cms file for conversion from DESMOND file format')
    group_in.add_argument('--gro_in', nargs=2, metavar='file',
            help='.gro and .top file for conversion from GROMACS file format')
    group_in.add_argument('--lmp_in', nargs=1, metavar='file',
            help='input file for conversion from LAMMPS file format (expects'
                    + ' data file in same directory and a read_data call)')

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
    group_misc.add_argument('-v', '--verbose', dest='verbose', action='store_true',
            help='high verbosity, includes DEBUG level output')

    # prints help if no arguments given
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    return parser

def main(args=''):
    parser = get_parser() # returns parser object
    if not args: # args automatically obtained from sys.argv
        args = parser.parse_args()
    else: # flexibility to call main() from other scripts
        args = parser.parse_args(args)

    if args.verbose:
        h.setLevel(logging.DEBUG)
    System._sys = System('InterMol')

    # --------------- PROCESS INPUTS ----------------- #

    if args.des_in: 
        prefix = args.des_in[0][args.des_in[0].rfind('/') + 1:-4]
        desmond_driver.readFile(args.des_in[0])
    elif args.gro_in:
        prefix = args.gro_in[0][args.gro_in[0].rfind('/') + 1:-4]
        # find the top file since order of inputs is not enforced
        top_in = [x for x in args.gro_in if x.endswith('.top')]
        assert(len(top_in)==1)
        top_in = top_in[0] # now just string instead of list
        # find the gro file since order of inputs is not enforced
        gro_in = [x for x in args.gro_in if x.endswith('.gro')]
        assert(len(gro_in)==1)
        gro_in = gro_in[0]
        gromacs_driver.readFile(top_in, gro_in, args.gropath)
    elif args.lmp_in:
        prefix = args.lmp_in[0][args.lmp_in[0].rfind('/') + 1:-4]
        lammps_driver.readFile(args.lmp_in[0])
    else:
        logger.error('No input file')
        sys.exit(1)

    # --------------- WRITE OUTPUTS ----------------- #

    if not args.oname: # default is to append _converted to the input prefix
        oname = '%s_converted' % prefix
    else:
        oname = args.oname
    oname = os.path.join(args.odir, oname) # prepend output directory to oname

    status = [Exception('Failed on write')]*N_FORMATS # records status of conversion
    if args.desmond:
        try: 
            desmond_driver.writeFile('{0}.cms'.format(oname))
            status[DES_POS] = 0 
        except Exception as e: # other output types will still be attempted
            logger.exception(e) # reports the supressed exception
            status[DES_POS] = e
    if args.gromacs:
        try:
            gromacs_driver.writeFile('{0}.top'.format(oname), '{0}.gro'.format(oname))
            status[GRO_POS] = 0
        except Exception as e:
            logger.exception(e)
            status[GRO_POS] = e
    if args.lammps:
        try:
            lammps_driver.writeFile('{0}.lmp'.format(oname))
            status[LMP_POS] = 0
        except Exception as e:
            logger.exception(e)
            status[LMP_POS] = e

    if args.energy:
        # evaluate input energy
        if args.des_in:
            input_type = 'Desmond' # used for displaying results
            e_in, e_infile = desmond_driver.desmond_energies(args.des_in[0],
                    'inputs/Desmond/onepoint.cfg', args.despath)
        elif args.gro_in:
            input_type = 'Gromacs'
            e_in, e_infile = gromacs_driver.gromacs_energies(top_in, gro_in,
                    'inputs/Gromacs/grompp.mdp', args.gropath, '')
        elif args.lmp_in:
            input_type = 'Lammps'
            # TODO: fix this when --lmp_in gets changed to read input files
            temp = args.lmp_in[0].split('.')[0] + '.input'
            e_in, e_infile = lammps_driver.lammps_energies(temp, args.lmppath)
        else:
            logger.warn('Something weird went wrong! Code should have never made it here.')
            pass # should never reach here

        # evaluate output energy
        output_type = []
        e_outfile = []
        e_out = []
        if args.desmond and status[DES_POS] == 0:
            output_type.append('Desmond')
            try:
                out, outfile = desmond_driver.desmond_energies('%s.cms' % oname,
                        'inputs/Desmond/onepoint.cfg', args.despath)
                status[DES_POS] = helper_functions.get_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)
            except Exception as e:
                status[DES_POS] = e
                logger.exception(e) #reports the supressed exception
                e_out.append(-1)
                e_outfile.append(-1)

        if args.gromacs and status[GRO_POS] == 0:
            output_type.append('Gromacs')
            try:
                out, outfile = gromacs_driver.gromacs_energies('%s.top' % oname,
                        '%s.gro' % oname, 'inputs/Gromacs/grompp.mdp', args.gropath, '')
                status[GRO_POS] = helper_functions.get_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)
            except Exception as e:
                status[GRO_POS] = e
                logger.exception(e)
                e_out.append(-1)
                e_outfile.append(-1)
 
        if args.lammps and status[LMP_POS] == 0:
            output_type.append('Lammps')
            try:
                out, outfile = lammps_driver.lammps_energies('%s.input' % oname,
                        args.lmppath)
                status[LMP_POS] = helper_functions.get_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)
            except Exception as e:
                status[LMP_POS] = e
                logger.exception(e)
                e_out.append(-1)
                e_outfile.append(-1)

        # display energy comparison results
        out = ['InterMol Conversion Energy Comparison Results','']
        out.append('{0} input energy file: {1}'.format(input_type, e_infile))
        for type, file in zip(output_type, e_outfile):
            out.append('{0} output energy file: {1}'.format(type, file))
        out += helper_functions.summarize_energy_results(e_in, e_out, input_type, output_type)
        logger.info('\n'.join(out))
    status = ['Converted' if x is 0 else x for x in status] # changes 0 to 'Converted'
    return status

if __name__ == '__main__':
    logger.info('Beginning InterMol conversion')
    main()
    logger.info('Finished!')
