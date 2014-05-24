import sys
import os
import argparse
from warnings import warn
import traceback
import pdb
import numpy as np

from intermol.System import System
from desmond_energies import desmond_energies
from gromacs_energies import gromacs_energies
from lammps_energies import lammps_energies
import helper_functions

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

def main(args=''):

    #TODO: fix verbosity

    parser = get_parser() # returns parser object
    if not args: # args automatically obtained from sys.argv
        args = parser.parse_args()
    else: # flexibility to call main() from other scripts
        args = parser.parse_args(args)


    # check to make sure output directory exists before doing anything
    if not os.path.exists(args.odir):
        raise Exception('Output directory not found: {0}'.format(args.odir))

    # check to make sure a conversion type is selected before doing anything
    if not args.desmond and not args.gromacs and not args.lammps:
        print 'ERROR: no conversion type selected'
        sys.exit(1)

    print 'Performing InterMol conversion...'
    System._sys = System('InterMol')
    print 'System initialized'


    # --------------- PROCESS INPUTS ----------------- #

    if args.des_in: # input type is desmond
        # check if input exists
        if not os.path.isfile(args.des_in[0]):
            raise Exception('File not found: {0}'.format(args.des_in[0]))
        prefix = args.des_in[0][args.des_in[0].rfind('/') + 1:-4]
        from intermol.DesmondExt.DesmondParser import DesmondParser as DesmondParser
        DesmondParser = DesmondParser()
        print 'Reading in Desmond structure {0}...'.format(args.des_in[0])
        try:
            DesmondParser.readFile(args.des_in[0], args.verbose)
        except Exception as err:
            print 'Failed on read'
            print traceback.format_exc()
            return 1 # failed on read, used in UnitTest.py
        print 'Sructure loaded\n'

    elif args.gro_in: # input type is gromacs
        # check if input exists
        for f in args.gro_in:
            if not os.path.isfile(f):
                raise Exception('File not found: {0}'.format(f))
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
        _ = gromacs_energies(top_in, gro_in,
                'Inputs/Gromacs/grompp.mdp', args.gropath, '',
                grompp_check=True)

        from intermol.GromacsExt.GromacsTopologyParser import GromacsTopologyParser
        import intermol.GromacsExt.GromacsStructureParser as GromacsStructureParser
        print "Reading in Gromacs topology {0}...".format(top_in)
        if not GromacsTopologyParser._GroTopParser:
            GromacsTopologyParser._GroTopParser = GromacsTopologyParser()
        try:
            GromacsTopologyParser._GroTopParser.parseTopology(top_in)
            print "Topology loaded\n"
            print "Reading in Gromacs structure {0}...".format(gro_in)
            GromacsStructureParser.readStructure(gro_in)
        except Exception as err:
            print 'Failed on read'
            print traceback.format_exc()
            return 1 # failed on read, used in UnitTest.py
        print "Structure loaded\n"

    elif args.lmp_in: # input type is lammps
        # check if file exists
        if not os.path.isfile(args.lmp_in[0]):
            raise Exception('File not found: {0}'.format(args.lmp_in[0]))
        prefix = args.lmp_in[0][args.lmp_in[0].rfind('/') + 1:-4]
        from intermol.lammps_extension.lammps_parser import LammpsParser
        print "Reading LAMMPS data & input files..."
        try:
            lammps_parser = LammpsParser()
            lammps_parser.read_system(args.lmp_in[0])
        except Exception as err:
            print 'Failed on read'
            print traceback.format_exc()
            return 1 # failed on read, used in UnitTest.py
        print "Data loaded\n"

    else:
        print 'No input file'
        sys.exit(1)

    # ------------ WRITE OUTPUTS -------------- #

    if not args.oname: # default is to append _converted to the input prefix
        oname = '%s_converted' % prefix
    else:
        oname = args.oname
    oname = os.path.join(args.odir, oname) # prepend output directory to oname

    try:
        if args.desmond:
            print 'Converting to Desmond...writing %s.cms...' % oname
            from intermol.DesmondExt.DesmondParser import DesmondParser
            DesmondParser = DesmondParser()
            DesmondParser.writeFile('%s.cms' % oname, args.verbose)

        if args.gromacs:
            print 'Converting to Gromacs...writing %s.gro, %s.top...' % (oname, oname)
            import intermol.GromacsExt.GromacsStructureParser as GromacsStructureParser
            GromacsStructureParser.writeStructure('%s.gro' % oname)
            from intermol.GromacsExt.GromacsTopologyParser import GromacsTopologyParser
            if not GromacsTopologyParser._GroTopParser:
                GromacsTopologyParser._GroTopParser = GromacsTopologyParser()
            GromacsTopologyParser._GroTopParser.writeTopology('%s.top' % oname)

        if args.lammps:
            print 'Converting to Lammps...writing %s.lmp...' % oname
            from intermol.lammps_extension.lammps_parser import LammpsParser
            lammps_parser = LammpsParser()
            lammps_parser.write('%s.lmp' % oname)

    except Exception as err:
        print 'Failed on write'
        print traceback.format_exc()
        return 2 # failed on write, used in UnitTest.py


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
            return 3 # failed at input energy, used in UnitTest.py

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
        print '{0} input energy file: {1}'.format(input_type, e_infile)
        for type, file in zip(output_type, e_outfile):
            print '{0} output energy file: {1}'.format(type, file)
        diff = helper_functions.print_multiple_energy_results(e_in, e_out, input_type, output_type)
        return diff # difference in potential energy for UnitTest.py script

    return 0 # converted sucessfully (no energy check)

if __name__ == '__main__':
    main()
