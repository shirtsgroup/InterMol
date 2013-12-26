import string
import sys
import os.path
from optparse import OptionParser
import intermol.Driver as Driver

#Error file
error = open("Errors.txt", "w")

#loop over all test cases for all x_to_y combinations
#List of Desmond files to be converted
filenames = ['a2a_dppc-out', 'protein', 'Rilpivirine', 'simulated_annealing_example', 'lck_Me_Cl_complex_12_1', 'lck_Me_Cl_solvent_12_1']

#Convert Desmond files to Desmond and Gromacs
for name in filenames:
    cms_in = os.path.join('Inputs/Desmond/', name, name+ '.cms')
    if not os.path.isfile(cms_in):
        #raise Exception ("File not found: {0}!".format(cms_in))
        error.write('\n(%s) -- File not found: {0}!' %(cms_in))
        pass
    try:
        Driver.load(cms_in)
    except:
        e = sys.exc_info()[0]
        error.write('\n(%s) -- %s' %(name, e))
        pass

    cms_out = os.path.join('Outputs', 'DesmondToDesmond', name, name+ '_OUT.cms')
    try:
        Driver.write(cms_out)
    except:
        e = sys.exc_info()[0]
        error.write('\n(%s) -- %s' %(name, e))
        pass

    gro_out = os.path.join('Outputs', 'DesmondToGromacs',  name+ '_OUT.gro')
    top_out = os.path.join('Outputs', 'DesmondToGromacs',  name+ '_OUT.top')
    try:
        Driver.write(top_out, gro_out)
    except:
        e = sys.exc_info()[0]
        error.write('\n(%s) -- %s' %(name, e))
        pass

#List of Gromacs files to be converted
filenames = ['2PPN', 'micelle']

#Converrt Gromacs files to Desmond and Gromacs
for name in filenames:
    gro_in = os.path.join('Inputs/Gromacs/', name, name+'.gro')
    if not os.path.isfile(gro_in):
        #raise Exception ("File not found: {0}!".format(gro_in))
        error.write('\n(%s) -- File not found: {0}!' %(cms_in))
        pass
    top_in = os.path.join('Inputs/Gromacs/', name, name+'.top')
    if not os.path.isfile(top_in):
        #raise Exception ("File not found: {0}!".format(top_in))
        error.write('\n(%s) -- File not found: {0}!' %(cms_in))
        pass
    try:
        Driver.load(top_in, gro_in)
    except:
        e = sys.exc_info()[0]
        error.write('\n(%s) -- %s' %(name, e))
        pass

    gro_out = os.path.join('Outputs', 'GromacsToGromacs', name, name+ '_OUT.gro')
    top_out = os.path.join('Outputs', 'GromacsToGromacs', name, name+ '_OUT.top')
    try:
        Driver.write(top_out, gro_out)
    except:
        e = sys.exc_info()[0]
        error.write('\n(%s) -- %s' %(name, e))
        pass

    cms_out = os.path.join('Outputs', 'GromacsToDesmond', name, name+ '_OUT.cms')
    try: 
        Driver.write(cms_out)
    except:
        e = sys.exc_info()[0]
        error.write('\n(%s) -- %s' %(name, e))
        pass

#Gromacs files whose directories and filenames aren't the same name
    #'complex' system located in micelle
complex_gro_in = os.path.join('Inputs/Gromacs/', 'micelle', 'complex.gro')
if not os.path.isfile(complex_gro_in):
    #raise Exception ("File not found: {0}!".format(complex_gro_in))
    error.write('\n(%s) -- File not found: {0}!' %(cms_in))
    pass
complex_top_in = os.path.join('Inputs/Gromacs/', 'micelle', 'complex.top')
if not os.path.isfile(complex_top_in):
    #raise Exception ("File not found: {0}!".format(complex_top_in))
    error.write('\n(%s) -- File not found: {0}!' %(cms_in))
    pass
try:
    Driver.load(complex_top_in, complex_gro_in)
except:
    e = sys.exc_info()[0]
    error.write('\n(%s) -- %s' %(name, e))
    pass

complex_gro_out = os.path.join('Outputs/GromacsToGromacs/micelle', 'complex_OUT.gro')
complex_top_out = os.path.join('Outputs/GromacsToGromacs/micelle', 'complex_OUT.top')
try:
    Driver.write(complex_top_out, complex_gro_out)
except:
    e = sys.exc_info()[0]
    error.write('\n(%s) -- %s' %(name, e))
    pass

complex_cms_out = os.path.join('Outputs/GromacsToDesmond/micelle', 'complex_OUT.cms')
try:
    Driver.write(complex_cms_out)
except:
    e = sys.exc_info()[0]
    error.write('\n(%s) -- %s' %(name, e))
    pass

#system_GMX located in system
GMX_gro_in = os.path.join('Inputs/Gromacs/system','system_GMX.gro')
if not os.path.isfile(GMX_gro_in):
    #raise Exception ("File not found: {0}!".format(GMX_gro_in))
    error.write('\n(%s) -- File not found: {0}!' %(cms_in))
    pass
GMX_top_in = os.path.join('Inputs/Gromacs/system','system_GMX.top')
if not os.path.isfile(GMX_top_in):
    #raise Exception ("File not found: {0}!".format(GMX_top_in))
    error.write('\n(%s) -- File not found: {0}!' %(cms_in))
    pass
try:
    Driver.load(GMX_top_in, GMX_gro_in)
except:
    e = sys.exc_info()[0]
    error.write('\n(%s) -- %s' %(name, e))
    pass

GMX_gro_out = os.path.join('Outputs/GromacsToGromacs/system', 'system_GMX_OUT.gro')
GMX_top_out = os.path.join('Outputs/GromacsToGromacs/system', 'system_GMX_OUT.top')
try:
    Driver.write(GMX_top_out, GMX_gro_out)
except:
    e = sys.exc_info()[0]
    error.write('\n(%s) -- %s' %(name, e))
    pass
GMX_cms_out = os.path.join('Outputs/GromacsToDesmond/system', 'system_GMX_OUT.cms')
try:
    Driver.write(GMX_cms_out)
except:
    e = sys.exc_info()[0]
    error.write('\n(%s) -- %s' %(name, e))
    pass

#system2_GMX located in system2
GMX2_gro_in = os.path.join('Inputs/Gromacs/system2','system2_GMX.gro')
if not os.path.isfile(GMX2_gro_in):
    #raise Exception ("File not found: {0}!".format(GMX2_gro_in))
    error.write('\n(%s) -- File not found: {0}!' %(cms_in))
    pass
GMX2_top_in = os.path.join('Inputs/Gromacs/system2','system2_GMX.top')
if not os.path.isfile(GMX2_top_in):
    #raise Exception ("File not found: {0}!".format(GMX2_top_in))
    error.write('\n(%s) -- File not found: {0}!' %(cms_in))
    pass
try:
    Driver.load(GMX2_top_in, GMX2_gro_in)
except:
    e = sys.exc_info()[0]
    error.write('\n(%s) -- %s' %(name, e))
    pass

GMX2_gro_out = os.path.join('Outputs/GromacsToGromacs/system2', 'system2_GMX_OUT.gro')
GMX2_top_out = os.path.join('Outputs/GromacsToGromacs/system2', 'system2_GMX_OUT.top')
try:
    Driver.write(GMX2_top_out, GMX2_gro_out)
except:
    e = sys.exc_info()[0]
    error.write('\n(%s) -- %s' %(name, e))
    pass

GMX2_cms_out = os.path.join('Outputs/GromacsToDesmond/system2', 'system2_GMX_OUT.cms')
try:
    Driver.write(GMX2_cms_out)
except:
    e = sys.exc_info()[0]
    error.write('\n(%s) -- %s' %(name, e))
    pass
