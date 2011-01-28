from mmtools.relativetools.checkrelativestates import *

ligand_name_list = [('LG%s'%i) for i in (2, 5, 6, 8, 9)]

# Gromacs Parameters
grompp = '/usr/local/gromacs4.5.1/bin/grompp_d'
mdrun = '/usr/local/gromacs4.5.1/bin/mdrun_d'
g_energy = '/usr/local/gromacs4.5.1/bin/g_energy_d'
mdp = '/home/mnz3v/mmtools/relativetools/Examples/FKBP/fep/LG2/solvent/minimize.mdp'
gros = list()
originaltops = list()
relativetops = list()
ligand_basepath = '/home/mnz3v/mmtools/relativetools/Examples/FKBP/ligands-parameterized'

# Gromacs input files
for i in ligand_name_list:
    gros.append('/home/mnz3v/mmtools/relativetools/Examples/FKBP/fep/%s/solvent/system_GMX.gro'%(i))
    originaltops.append('/home/mnz3v/mmtools/relativetools/Examples/FKBP/fep/%s/solvent/system_GMX.top'%(i))

for i in ligand_name_list:
    for j in ligand_name_list:
        if j > i:
            relativetops.append('/home/mnz3v/mmtools/relativetools/Examples/FKBP/fep/%s/solvent/system_GMX.top'%(i))
            relativetops.append('/home/mnz3v/mmtools/relativetools/Examples/FKBP/fep/%s/solvent/system_GMX.top'%(j))


# Pairwise matches by atom index
system = checkrelativestates(grompp, mdrun, g_energy, ligand_name_list, mdp, gros, originaltops, relativetops, ligand_basepath = ligand_basepath)

