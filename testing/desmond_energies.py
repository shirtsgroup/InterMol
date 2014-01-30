import os
import pdb
import shutil

def desmond_energies(name, cms=None, in_out = 'DtoD', despath='/opt/schrodinger2013'):
    """

    despath = path to DESMOND binaries

    """
    cfg = 'Inputs/Desmond/onepoint.cfg'
    if name == None:
        name = 'onepoint'
    if in_out == 'in':
        base = 'Inputs/Desmond'
        if cms == None:
            cms = os.path.join(base, name, 'desmond.cms')
    elif in_out == 'DtoD':
        base = 'Outputs/DesmondToDesmond'
        if cms == None:
            base = os.path.join(base, name, 'desmond.cms')
    else:
        raise Exception("Unknown flag: {0}".format(in_out))

    os.environ['SCHRODINGER'] = despath
    # run desmond
    os.system('$SCHRODINGER/desmond' + " -WAIT -P 1 -in {cms} -JOBNAME {name} -c {cfg}".format(name = name, cms=cms,cfg=cfg))

    energrp =  name + '.enegrp.dat'

    # extract energy output and parse initial energies
    with open(energrp) as f:
        all_lines = f.readlines()

    # move output files to the correct directory
    files = ['checkpt.cpt',name + '-out.cfg', name + '-out.cms', name + '.ene',
             name + '-out.idx', name + '.log', name + '_simbox.dat', energrp,'trj']
    for file in files:
        file_out = os.path.join(base,name,file)
        if os.path.isfile(file_out):
            os.remove(file_out)
        elif os.path.isdir(file_out):
            shutil.rmtree(file_out)
        os.rename(file,file_out)
    # plus the extra .cms file that gets copied to the current directory
    # will need to tweak this line to be properly general.
    os.remove('desmond.cms')

    zerolines = []
    for line in all_lines:
        vals = line.split()
        # return all the lines for the first timepoint that dont start with a space.
        if (vals[1] == '(0.000000)' and vals[0][0] != ' '):
            zerolines.append(line)

    return zerolines
