from mmtools.relativetools.checkrelativestates import *

ligand_name_list = [('jnk.aff-%s'%i) for i in (13, 16, 53)]

# Gromacs Parameters
grompp = '/usr/local/gromacs4.5.1/bin/grompp_d'
mdrun = '/usr/local/gromacs4.5.1/bin/mdrun_d'
g_energy = '/usr/local/gromacs4.5.1/bin/g_energy_d'
mdp = '/home/mnz3v/mmtools/relativetools/Examples/JNK3/fep/jnk.aff-13/solvent/minimize.mdp'
gros = list()
originaltops = list()
relativetops = list()
ligand_basepath = '/home/mnz3v/mmtools/relativetools/Examples/JNK3/ligands-parameterized'

# Gromacs input files
for i in ligand_name_list:
    gros.append('/home/mnz3v/mmtools/relativetools/Examples/JNK3/fep/%s/solvent/system_GMX.gro'%(i))
    originaltops.append('/home/mnz3v/mmtools/relativetools/Examples/JNK3/fep/%s/solvent/system_GMX.top'%(i))

for i in ligand_name_list:
    for j in ligand_name_list:
        if j > i:
            relativetops.append('/home/mnz3v/mmtools/relativetools/Examples/JNK3/fep/%s/solvent/system_GMX.top'%(i))
            relativetops.append('/home/mnz3v/mmtools/relativetools/Examples/JNK3/fep/%s/solvent/system_GMX.top'%(j))


# Pairwise matches by atom index
system = checkrelativestates(grompp, mdrun, g_energy, ligand_name_list, mdp, gros, originaltops, relativetops, ligand_basepath = ligand_basepath)

