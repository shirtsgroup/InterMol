import string
import os.path
import intermol.Driver as Driver

print "\nStarting Program? (Desmond or Gromacs)"
program_in = raw_input("")
print "\nConvert to? (Desmond or Gromacs)"
program_out = raw_input("")
print "\nSystem name?"
name = raw_input("")

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
