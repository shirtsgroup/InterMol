import glob
import numpy
import convert
import os
import pdb

INPUT_DIR = './Inputs/Desmond/RegTest'
OUTPUT_DIR = 'RegTest_Output'

def main():
    if not os.path.isdir(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)

    files = glob.glob('%s/*.cms' % INPUT_DIR) # return list of files that match the string
    results = []

    # for now only attempts desmond to desmond conversion of input files
    for f in files:
        try:
            args = ['--des_in', f, '--desmond', '--odir', OUTPUT_DIR, '-e']
            print 'converting %s to desmond file format with command' % f
            print 'python convert.py', ' '.join(args)
            rms = convert.main(args) # reuses code from convert.py
            results.append(rms) 
        except:
            results.append(-1)
    print 'RMS errors between input and output energies'
    print results

if __name__ == '__main__':
    main()




