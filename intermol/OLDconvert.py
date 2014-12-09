import sys
import os
import argparse
import logging
import warnings
import numpy as np
from collections import OrderedDict
#from intermol.desmond_extension import desmond_driver
from intermol.gromacs import gromacs_driver
#from intermol.lammps_extension import lammps_driver
from intermol.system import System

# for debugging
import pdb

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

    logging.captureWarnings(True) # redirect warnings module messages to logging system
    warning_logger = logging.getLogger('py.warnings') 
    warning_logger.addHandler(h)

def parse_args(args):
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
    group_misc.add_argument('-dp', '--despath', dest='despath',
            metavar='path', default='',
            help='path for DESMOND binary, needed for energy evaluation')
    group_misc.add_argument('-gp', '--gropath', dest='gropath',
            metavar='path', default='',
            help='path for GROMACS binary, needed for energy evaluation')
    group_misc.add_argument('-lp', '--lmppath', dest='lmppath',
            metavar='path', default='lmp_openmpi',
            help='path for LAMMPS binary, needed for energy evaluation')
    group_misc.add_argument('-f', '--force', dest='force', action='store_true',
            help='ignore warnings (NOTE: may lead to partially correct topologies)')
    group_misc.add_argument('-v', '--verbose', dest='verbose', action='store_true',
            help='high verbosity, includes DEBUG level output')

    # prints help if no arguments given
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    if args.verbose:
        h.setLevel(logging.DEBUG)
    warnings.simplefilter("always") # print warnings
    if not args.force:
        warnings.simplefilter("error") # warnings will be treated as exceptions unless force flag is used
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
        gromacs_driver.read_file(top_in, gro_in, args.gropath)
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
            temp = args.lmp_in[0].rsplit('.',1)[0] + '.input'   # split from the right, only one split to avoid leading './'
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
                status[DES_POS] = get_diff(e_in, out)
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
                status[GRO_POS] = get_diff(e_in, out)
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
                status[LMP_POS] = get_diff(e_in, out)
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
        out += summarize_energy_results(e_in, e_out, input_type, output_type)
        logger.info('\n'.join(out))
    status = ['Converted' if x is 0 else x for x in status] # changes 0 to 'Converted'
    return status # for systems_test.py

# helper functions

def get_diff(e_in, e_out):
    """returns difference in potential energy
    
    arguments:
        e_in  - dictionary of energy groups from input file
        e_out - dictionary of energy groups from output file
    
    returns:
        potential energy difference in units of the input
    """
    type = 'Potential' # getting difference in potential energy
    input = e_in[type]
    diff = e_out[type].in_units_of(input.unit) - input
    return diff._value

def find_match(key, dict, unit):
    """helper function for summarize_energy_results() """
    if key in dict:
        return dict[key].in_units_of(unit)._value
    else:
        return np.nan

def summarize_energy_results(energy_input, energy_outputs, input_type, output_types):
    """ creates a table comparing input and output energy groups

    args:
        energy_input  - dictionary of energy groups from input file
        energy_output - list containing dictionary of energy groups or -1 for each output file
        input_type    - input format string
        output_types  - list containing output formats
    
    returns:
        out - list of strings which forms a summary table using "\n".join(out) 
    """
    out = []
    # remove failed evaluations (-1 in energy_outputs)
    failed_i = [i for i,x in enumerate(energy_outputs) if x == -1]
    failed = [output_types[i] for i in failed_i]
    output_types = [x for i,x in enumerate(output_types) if i not in failed_i]
    energy_outputs = [x for x in energy_outputs if x != -1]
    # find all distinct labels
    labels = set(energy_input.keys())
    for e_out in energy_outputs:
        for key, value in e_out.iteritems():
            labels.add(key)
    # set up energy comparison table
    labels = list(labels)
    unit = energy_input[energy_input.keys()[0]].unit
    energy_all = [energy_input] + energy_outputs
    data = np.empty((len(labels),len(energy_all)))
    for i in range(len(data)):
        for j in range(len(energy_all)):
            data[i,j] = find_match(labels[i], energy_all[j], unit)
    # sort table <= to do
    # print table
    out.append('')
    out.append('energy group summary')
    out.append('=======================================================================')
    header = '%20s %18s ' % ('type', 'input (%s)' % input_type)
    for otype in output_types:
        header += '%37s' %('output (%s) diff (%s)' % (otype, otype))
    out.append(header)
    for i in range(len(data)):
        line = '%20s ' % labels[i]
        line += '%18.8f ' % data[i][0]
        for j in range(1,len(data[i])):
            line += '%18.8f %18.8f' % (data[i][j],data[i][j]-data[i][0])
        out.append(line)
    out.append('')
    # get differences in potential energy
    i = labels.index('Potential')
    diff = data[i,1::] - data[i,0]
    for d, otype in zip(diff, output_types):
        out.append('difference in potential energy from %s=>%s conversion: %18.8f' 
                    % (input_type, otype, d))
    for f in failed:
        out.append('energy comparison for {0} output failed'.format(f))
    out.append('=======================================================================')
    return out

if __name__ == '__main__':
    logger.info('Beginning InterMol conversion')
    main()
    logger.info('Finished!')
