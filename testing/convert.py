import sys
import os
import pdb
import numpy as np
import argparse

import intermol.Driver as Driver
import evaluate

def parse_args():
    parser = argparse.ArgumentParser(description = 'testing file conversions')
    parser.add_argument('-i', '--in', dest='infile', help="input struture/topology .cms file")
    parser.add_argument('-o', '--out', dest='outfile', help="output structure/topology .cms file")
    parser.add_argument('--cfg', dest='cfgfile', default='Inputs/Desmond/onepoint.cfg', help="input .cfg file (for DESMOND)")
    parser.add_argument('-d', dest='despath', default='/opt/schrodinger2013/', help="path for DESMOND binary")
    parser.add_argument('-e', dest='energy', action="store_true", help="evaluate energies",default=False)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    # check if files exist
    if not os.path.isfile(args.infile):
	    raise Exception("File not found: {0}!".format(args.infile))
    outdir = os.path.dirname(args.outfile)
    if not os.path.exists(outdir):
	    raise Exception("Output directory not found: {0}!".format(outdir))

    # perform conversion
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

if __name__ == "__main__":
    main()
