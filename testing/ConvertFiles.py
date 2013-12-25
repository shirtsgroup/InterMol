import string
import os.path
from optparse import OptionParser
import intermol.Driver as Driver

parser = OptionParser()
parser.add_option('-i', type='str', dest='input', help = "Specify input program")
parser.add_option('-o', type='str', dest='output', help = "Specify output program")
parser.add_option('-f', type='str', dest='files', help = "Specify file to be converted")

(options, args) = parser.parse_args()
program_in = options.input
program_out = options.output
filename = options.files
args.append(filename)
    
print "Starting in %s and finishing in %s" %(program_in, program_out)
for name in args:
    
    Driver.initSystem(name)

    if  'Gromacs' == program_in:
        gro_in = os.path.join('Inputs/Gromacs/', name+'.gro')
        if not os.path.isfile(gro_in):
            raise Exception ("File not found: {0}!".format(gro_in))
        top_in = os.path.join('Inputs/Gromacs/', name+'.top')
        if not os.path.isfile(top_in):
            raise Exception ("File not found: {0}!".format(top_in))
        Driver.load(top_in, gro_in)    
    elif 'Desmond' == program_in:
        cms_in = os.path.join('Inputs/Desmond/', name, name+ '.cms')
        if not os.path.isfile(cms_in):
            raise Exception ("File not found: {0}!".format(cms_in))
        Driver.load(cms_in)
    
    if 'Gromacs' == program_out:
        gro_out = os.path.join('Outputs', program_in + 'to' + program_out,  name+ '_OUT.gro')
        gro_out = os.path.join('Outputs', program_in + 'to' + program_out,  name+ '_OUT.top')
        Driver.write(top_out, gro_out)
    elif 'Desmond' == program_out:
        cms_out = os.path.join('Outputs', program_in + 'to' + program_out, name, name+ '_OUT.cms')
        Driver.write(cms_out)
print "$s converted" %(args)
