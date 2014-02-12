import os
import pdb
import numpy as np
import intermol.Driver as Driver
from desmond_energies import desmond_energies

def desmond_to_desmond(name='2PPN', cms=None, despath=None,
        energy=True, clean=True):
    """Test gromacs to gromacs conversion
    """

    if cms == None:
        cms = os.path.join(name,'desmond.cms')
    cms_in = os.path.join('Inputs/Desmond/', cms)    

    if not os.path.isfile(cms_in):
	    raise Exception("File not found: {0}!".format(cms_in))
    
    cms_out = os.path.join('Outputs/DesmondToDesmond/', name, 'desmond-converted.cms')

    #calc input energies
    if energy:
        elines_in = desmond_energies(name,cms_in,'in',despath)

    Driver.initSystem(name)
    Driver.load(cms_in)
    Driver.write(cms_out)

    pdb.set_trace()
    # calc output energies
    if energy:
        elines_out = desmond_energies(name,cms_out,'DtoD',despath)    

    e_in  = np.zeros([len(elines_in)],float)
    e_out = np.zeros([len(elines_out)],float)
    for i, e in enumerate(elines_in): 
        eins  = elines_in[i].split() 
        eouts = elines_out[i].split()
        e_in[i] = float(eins[1])
        e_out[i] = float(eouts[1])
        enames.append()

    diff = e_in - e_out
    rms = np.sqrt(np.mean(diff ** 2))
    return rms, enames, e_in, e_out

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-n', type='str', dest='name', default='2PPN',
            help="Name of system")
    parser.add_option('-p', type='str', dest='cms', default=None,
            help="Struture/Topology .cms file")
    parser.add_option('-d', type='str', dest='despath', default='/opt/schrodinger2013/',
            help="path for DESMOND binary")
    parser.add_option('-x', type='int', dest='clean', default=1,
            help="Clean backup files produced by DESMOND")
    parser.add_option('-e', type='int', dest='energy', default=1,
            help="Evaluate energies")

    (options, args) = parser.parse_args()
    name = options.name
    cms = options.cms
    despath = options.despath
    energy = options.energy
    clean = options.clean

    rms, types, e_in, e_out = desmond_to_desmond(name,cms, despath, energy, clean)

    print "======================================================================="
    print "Summary statistics"
    #types = ['Bond', 'Angle', 'Proper Dih.', 'Ryckaert-Bell.', 'LJ-14', 'Coulomb-14',
    #        'LJ (SR)', 'Disper. corr.', 'Coulomb (SR)', 'Coul. recip.', 'Potential',
    #        'Kinetic En.', 'Total Energy', 'Temperature']

    print "%20s %12s %12s %12s" % ("Type", "Input", "Output", "Diff")
    for name, i, o in zip(types, e_in, e_out):
        print "%20s %15.8f %15.8f %15.8f" % (name, i, o, i-o)

    print " "
    print "RMS signed error: %10.5f" % (rms)
    
