import sys
import os
import argparse
from warnings import warn
import logging
import traceback
import pdb
import numpy as np

from intermol.System import System
from desmond_energies import desmond_energies
from gromacs_energies import gromacs_energies
from lammps_energies import lammps_energies
import helper_functions
from user_exceptions import *
import desmond_driver
import gromacs_driver
import lammps_driver

# global constants

DES_POS = 0
GRO_POS = 1
LMP_POS = 2
N_FILETYPES = 3

# Make a global logging object.
logger = logging.getLogger('InterMolLog')
logger.setLevel(logging.INFO) # specifies lowest severity log messages to handle (DEBUG will be ignored)
h = logging.StreamHandler()
#f = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s")
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
            help='verbosity')

    # prints help if no arguments given
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    return parser

def main_simple():
    args = get_parser().parse_args()
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

    status = [-1]*N_FILETYPES # records status of conversion
    if args.desmond:
        try: 
            desmond_driver.writeFile('{0}.cms'.format(oname))
            status[DES_POS] = 0 
        except Exception as e: # other output types will still be attempted
            logger.exception(e) # reports the supressed exception
    if args.gromacs:
        try:
            gromacs_driver.writeFile('{0}.top'.format(oname), '{0}.gro'.format(oname))
            status[GRO_POS] = 0
        except Exception as e:
            logger.exception(e)
    if args.lammps:
        try:
            lammps_driver.writeFile('{0}.lmp'.format(oname))
            status[LMP_POS] = 0
        except Exception as e:
            logger.exception(e)

    if args.energy:
        # evaluate input energy
        if args.des_in:
            input_type = 'Desmond' # used for displaying results
            e_in, e_infile = desmond_driver.desmond_energies(args.des_in[0],
                    'Inputs/Desmond/onepoint.cfg', args.despath)
        elif args.gro_in:
            input_type = 'Gromacs'
            e_in, e_infile = gromacs_driver.gromacs_energies(top_in, gro_in,
                    'Inputs/Gromacs/grompp.mdp', args.gropath, '')
        elif args.lmp_in:
            input_type = 'Lammps'
            # TODO: fix this when --lmp_in gets changed to read input files
            temp = args.lmp_in[0].split('.')[0] + '.input'
            e_in, e_infile = lammps_driver.lammps_energies(temp, args.lmppath)
        else:
            warn('Something weird went wrong! Code should have never made it here.')
            pass # should never reach here

        # evaluate output energy
        output_type = []
        e_outfile = []
        e_out = []
        if args.desmond and status[DES_POS] == 0:
            output_type.append('Desmond')
            try:
                out, outfile = desmond_driver.desmond_energies('%s.cms' % oname,
                        'Inputs/Desmond/onepoint.cfg', args.despath)
                e_out.append(out)
                e_outfile.append(outfile)
            except Exception as e:
                logger.exception(e) #reports the supressed exception
                e_out.append(-1)
                e_outfile.append(-1)

        if args.gromacs and status[GRO_POS] == 0:
            output_type.append('Gromacs')
            try:
                out, outfile = gromacs_driver.gromacs_energies('%s.top' % oname,
                        '%s.gro' % oname, 'Inputs/Gromacs/grompp.mdp', args.gropath, '')
                e_out.append(out)
                e_outfile.append(outfile)
            except Exception as e:
                logger.exception(e)
                e_out.append(-1)
                e_outfile.append(-1)
 
        if args.lammps and status[LMP_POS] == 0:
            output_type.append('Lammps')
            try:
                out, outfile = lammps_driver.lammps_energies('%s.input' % oname,
                        args.lmppath)
                e_out.append(out)
                e_outfile.append(outfile)
            except Exception as e:
                logger.exception(e)
                e_out.append(-1)
                e_outfile.append(-1)

        # display energy comparison results
        print ''
        print '{0} input energy file: {1}'.format(input_type, e_infile)
        for type, file in zip(output_type, e_outfile):
            print '{0} output energy file: {1}'.format(type, file)
        helper_functions.print_multiple_energy_results(e_in, e_out, input_type, output_type)

def main(args=''):

    parser = get_parser() # returns parser object
    if not args: # args automatically obtained from sys.argv
        args = parser.parse_args()
    else: # flexibility to call main() from other scripts
        args = parser.parse_args(args)


    # check to make sure output directory exists before doing anything
    if not os.path.exists(args.odir):
        logger.error('Output directory not found: {0}'.format(args.odir))
        sys.exit(1)

    # check to make sure a conversion type is selected before doing anything
    if not args.desmond and not args.gromacs and not args.lammps:
        logger.error('no conversion type selected')
        sys.exit(1)

    System._sys = System('InterMol')

    # --------------- PROCESS INPUTS ----------------- #

    if args.des_in: # input type is desmond
        # check if input exists
        if not os.path.isfile(args.des_in[0]):
            logger.error('File not found: {0}'.format(args.des_in[0]))
            raise ReadError('File not found')
        prefix = args.des_in[0][args.des_in[0].rfind('/') + 1:-4]
        from intermol.DesmondExt.DesmondParser import DesmondParser as DesmondParser
        DesmondParser = DesmondParser()
        logger.info('Reading in Desmond structure {0}...'.format(args.des_in[0]))
        try:
            DesmondParser.readFile(args.des_in[0], args.verbose)
            logger.info('Sructure loaded')
        except Exception as err:
            logger.exception(err)
            raise ReadError('Failed on read')

    elif args.gro_in: # input type is gromacs
        # check if input exists
        for f in args.gro_in:
            if not os.path.isfile(f):
                logger.exception('File not found: {0}'.format(f))
                raise ReadError('File not found')
        # find the top file since order of inputs is not enforced
        top_in = [x for x in args.gro_in if x.endswith('.top')]
        assert(len(top_in)==1)
        top_in = top_in[0] # now just string instead of list
        # find the gro file since order of inputs is not enforced
        gro_in = [x for x in args.gro_in if x.endswith('.gro')]
        assert(len(gro_in)==1)
        gro_in = gro_in[0]
        prefix = args.gro_in[0][args.gro_in[0].rfind('/') + 1:-4]
        # ensure .gro and .top are a valid match
        try:
            gromacs_energies(top_in, gro_in,
                'Inputs/Gromacs/grompp.mdp', args.gropath, '',
                grompp_check=True)
        except Exception as err:
            logger.exception(err)
            raise ReadError('grompp check failed')

        from intermol.GromacsExt.GromacsTopologyParser import GromacsTopologyParser
        import intermol.GromacsExt.GromacsStructureParser as GromacsStructureParser
        logger.info('Reading in Gromacs topology {0}...'.format(top_in))
        GromacsTopologyParser._GroTopParser = GromacsTopologyParser()
        try:
            GromacsTopologyParser._GroTopParser.parseTopology(top_in)
            logger.info('Topology loaded')
            logger.info('Reading in Gromacs structure {0}...'.format(gro_in))
            GromacsStructureParser.readStructure(gro_in)
            logger.info('Structure loaded')
        except Exception as err:
            logger.exception(err)
            raise ReadError('Failed on read')

    elif args.lmp_in: # input type is lammps
        # check if file exists
        if not os.path.isfile(args.lmp_in[0]):
            logger.exception('File not found: {0}'.format(args.lmp_in[0]))
            raise ReadError('File not found')
        prefix = args.lmp_in[0][args.lmp_in[0].rfind('/') + 1:-4]
        from intermol.lammps_extension.lammps_parser import LammpsParser
        logger.info('Reading LAMMPS data & input files...')
        try:
            lammps_parser = LammpsParser()
            lammps_parser.read_system(args.lmp_in[0])
            logger.info('Structure loaded')
        except Exception as err:
            logger.exception(err)
            raise ReadError('Failed on read')

    else:
        logger.error('No input file')
        sys.exit(1)

    # ------------ WRITE OUTPUTS -------------- #

    results = [0]*N_FILETYPES # 
    if not args.oname: # default is to append _converted to the input prefix
        oname = '%s_converted' % prefix
    else:
        oname = args.oname
    oname = os.path.join(args.odir, oname) # prepend output directory to oname

    if args.desmond:  # Desmond OK so far
        results[DES_POS] == 2 # assume bad on write
        try:
            logger.info('Converting to Desmond...')
            logger.info('Writing %s.cms...' % oname)
            from intermol.DesmondExt.DesmondParser import DesmondParser
            DesmondParser = DesmondParser()
            DesmondParser.writeFile('%s.cms' % oname, args.verbose)
            results[DES_POS] = 0 # didn't fail on write, OK so far
        except Exception as err:
            print 'Failed on write for Desmond'
            print traceback.format_exc()

    if args.gromacs:  #Gromacs OK so far
        results[GRO_POS] == 2 # assume bad on write
        try:
            print 'Converting to Gromacs...writing %s.gro, %s.top...' % (oname, oname)
            import intermol.GromacsExt.GromacsStructureParser as GromacsStructureParser
            GromacsStructureParser.writeStructure('%s.gro' % oname)
            from intermol.GromacsExt.GromacsTopologyParser import GromacsTopologyParser
            if not GromacsTopologyParser._GroTopParser:
                GromacsTopologyParser._GroTopParser = GromacsTopologyParser()
            GromacsTopologyParser._GroTopParser.writeTopology('%s.top' % oname)
            results[GRO_POS] = 0
        except Exception as err:
            print 'Failed on write for Gromacs'
            print traceback.format_exc()

    if args.lammps:  # LAMMPS OK so far
        results[LMP_POS] == 2 # assume bad on write
        try:
            print 'Converting to Lammps...writing %s.lmp...' % oname
            from intermol.lammps_extension.lammps_parser import LammpsParser
            lammps_parser = LammpsParser()
            lammps_parser.write('%s.lmp' % oname)
            results[LMP_POS] = 0
        except Exception as err:
            print 'Failed on write for LAMMPS'
            print traceback.format_exc()

    # ------------ ENERGY COMPARISON ----------- #
    # to do: remove hardcoded parmeter filenames
    #        format printing of energy dictionary
    if args.energy:
        # evaluate input energy
        try:
            if args.des_in:
                input_type = 'Desmond' # used for displaying results
                e_in, e_infile = desmond_energies(args.des_in[0],
                        'Inputs/Desmond/onepoint.cfg', args.despath)
            elif args.gro_in:
                input_type = 'Gromacs'
                e_in, e_infile = gromacs_energies(top_in, gro_in,
                        'Inputs/Gromacs/grompp.mdp', args.gropath, '')
            elif args.lmp_in:
                input_type = 'Lammps'
                # TODO: fix this when --lmp_in gets changed to read input files
                temp = args.lmp_in[0].split('.')[0] + '.input'
                e_in, e_infile = lammps_energies(temp, args.lmppath)
            else:
                warn('Something weird went wrong! Code should have never made it here.')
                pass # should never reach here
        except Exception as err:
            print 'Failed at evaluating energy of input file'
            print traceback.format_exc()
            return [3]*N_FILETYPES # failed at input energy, used in UnitTest.py
                                   # input is poorly formatted, overwrite the other errors.
                                   # We can't expect conversion if the input energies are bad.

        # evaluate output energy
        output_type = []
        e_outfile = []
        e_out = []
        if args.desmond:
            output_type.append('Desmond')
            try:
                out, outfile = desmond_energies('%s.cms' % oname,
                        'Inputs/Desmond/onepoint.cfg', args.despath)
                e_out.append(out)
                e_outfile.append(outfile)
            except Exception as e:
                print 'Failed at evaluating energy of Desmond output file: {0}'.format(e)
                print traceback.format_exc()
                e_out.append(-1)
                e_outfile.append(-1)

        if args.gromacs:
            output_type.append('Gromacs')
            try:
                out, outfile = gromacs_energies('%s.top' % oname,
                        '%s.gro' % oname, 'Inputs/Gromacs/grompp.mdp', args.gropath, '')
                e_out.append(out)
                e_outfile.append(outfile)
            except Exception as e:
                print 'Failed at evaluating energy of Gromacs output file: {0}'.format(e)
                print traceback.format_exc()
                e_out.append(-1)
                e_outfile.append(-1)
 
        if args.lammps:
            output_type.append('Lammps')
            try:
                out, outfile = lammps_energies('%s.input' % oname,
                        args.lmppath)
                e_out.append(out)
                e_outfile.append(outfile)
            except Exception as e:
                print 'Failed at evaluating energy of output file: {0}'.format(e)
                print traceback.format_exc()
                e_out.append(-1)
                e_outfile.append(-1)

        if len(e_out) > 1 and e_out[1:] == e_out[:-1]: # all failed to evaluate (all are -1)
            return [4]*len(output_type) # 4 is errorcode for UnitTest.py script

        # display energy comparison results
        print ''
        print '{0} input energy file: {1}'.format(input_type, e_infile)
        for type, file in zip(output_type, e_outfile):
            print '{0} output energy file: {1}'.format(type, file)
        diff = helper_functions.print_multiple_energy_results(e_in, e_out, input_type, output_type)
        return diff # difference in potential energy for UnitTest.py script

    return results # converted sucessfully (no energy check)

if __name__ == '__main__':
    logger.info('Beginning InterMol conversion')
    main_simple()
    logger.info('Finished!')
