import sys
import os
import pdb
import numpy as np
import argparse

import intermol.Driver as Driver
import evaluate

def parse_args():
    parser = argparse.ArgumentParser(description = 'Perform a file conversion')

    # input arguments
    group_in = parser.add_argument_group('Choose input conversion format')
    group_in.add_argument('--desmond', nargs=2, metavar='file', help='.cms and .cfg file for conversion from DESMOND file format')
    group_in.add_argument('--gromacs', nargs=2, metavar='file', help='.gro and .top file for conversion from GROMACS file format')
    group_in.add_argument('--lammps', nargs=2, metavar='file', help='.lmp and .input file for conversion from LAMMPS file format')

    # output arguments
    group_out = parser.add_argument_group('Choose output conversion format(s)')
    group_out.add_argument('--cms_out', dest='cms_out', metavar='.CMS', help='output .cms file (for DESMOND conversion)')
    group_out.add_argument('--gro_out', dest='gro_out', metavar='.GRO', help='output .gro file (for GROMACS conversion)')
    group_out.add_argument('--top_out', dest='top_out', metavar='.TOP', help='output .top file (for GROMACS conversion)')
    group_out.add_argument('--lmp_out', dest='lmp_out', metavar='.LMP', help='output .lmp file (for LAMMPS conversion)')

# 
#    parser_des_in = parser.add_argument_group('from_desmond', 'conversion from desmond file format')
#
#
#    parser_des_in = parser.add_argument_group('from_desmond', 'conversion from desmond file format')
#    parser_des_in.add_argument('--cms_in', dest='cms_in', help='input .cms file (for DESMOND)')
#    parser_des_in.add_argument('--cfg_in', dest='cfg_in', default='Inputs/Desmond/onepoint.cfg', help='input .cfg file (for DESMOND)')
#
#    parser_gro_in = parser.add_argument_group('from_gromacs', 'conversion from gromacs file format')
#    parser_gro_in.add_argument('--top_in', dest='top_in', help='input .top file (for GROMACS)')
#    parser_gro_in.add_argument('--gro_in', dest='gro_in', help='input .gro file (for GROMACS)')
#
#    parser_lam_in = parser.add_argument_group('from_lammps', 'conversion from lammps file format')
#    parser_lam_in.add_argument('--lmp_in', dest='lmp_in', help='input .lmp file (for LAMMPS)')
#    parser_lam_in.add_argument('--input_in', dest='input_in', help='input .input file (for LAMMPS)')
#
#    parser_des_out = parser.add_argument_group('to_desmond', 'conversion to desmond file format')
#    parser_des_out.add_argument('--cms_out', dest='cms_out', help='output .cms file (for DESMOND)')
#
#    parser_gro_out = parser.add_argument_group('to_gromacs', 'conversion to gromacs file format')
#    parser_gro_out.add_argument('--top_out', dest='top_out', help='output .top file (for GROMACS)')
#    parser_gro_out.add_argument('--gro_out', dest='gro_out', help='output .gro file (for GROMACS)')
#
#    parser_lam_out = parser.add_argument_group('to_lammps', 'conversion to lammps file format')
#    parser_lam_out.add_argument('--lmp_out', dest='lmp_out', help='output .lmp file (for LAMMPS)')
#

#    parser.add_argument('-f', '--from', dest='from', choices=['desmond','gromacs','lammps'])
#    parser.add_argument('-i', '--in', dest='infile', nargs='+', help='input file name')
#    parser.add_argument('-t', '--to', dest='to', choices=['desmond','gromacs','lammps'])
#    parser.add_argument('-o', '--out', dest='outfile', nargs='+', help='output file name')

    # other
    parser.add_argument('-e', dest='energy', action='store_true', help='evaluate energies',default=False)
    # is this argument needed? gromacs/lammps path missing right now
    parser.add_argument('-d', dest='despath', default='/opt/schrodinger2013/', help='path for DESMOND binary')

    return parser.parse_args()

def main():
    args = parse_args()

    # check if files exist
    if not os.path.isfile(args.infile):
	    raise Exception('File not found: {0}!'.format(args.infile))
    outdir = os.path.dirname(args.outfile)
    if not os.path.exists(outdir):
        raise Exception('Output directory not found: {0}!'.format(outdir))


    # perform conversion - only DtoD for now....
    try:
        print 'Performing InterMol conversion:'
        print '    From: {0}'.format(args.infile)
        print '    To: {0}'.format(args.outfile)

        name = 'system'
        Driver.initSystem(name) # what does this name do?

        Driver.load(args.infile)

        Driver.write(args.outfile)
    except:
        print 'Failed at InterMol conversion'
        sys.exit(1)

    # calc input energies
    if args.energy:
        e_in, e_infile = evaluate.desmond_energies(args.infile, args.cfgfile, args.despath)
        print 'Total energy from %s:' % e_infile
        print e_in
        e_out, e_outfile = evaluate.desmond_energies(args.outfile, args.cfgfile, args.despath)    
        print 'Total energy from %s:' % e_outfile
        print e_out

if __name__ == '__main__':
    main()
