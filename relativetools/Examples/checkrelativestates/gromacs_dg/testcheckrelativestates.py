from mmtools.relativetools.checkrelativestates import *

ligand_name_list = [('LG%s'%i) for i in (2, 3, 5, 6, 8, 9)]

# Gromacs Parameters
grompp = '/usr/local/gromacs_dg/bin/grompp_d'
mdrun = '/usr/local/gromacs_dg/bin/mdrun_d'
g_energy = '/usr/local/gromacs_dg/bin/g_energy_d'
mdp = '/home/mnz3v/mmtools/relativetools/Examples/checkrelativestates/gromacs_dg/minimize.mdp'
gros = list()
originaltops = list()
originalitps = list()
relativetops = list()
relativeitps = list()
matches = list()

# Gromacs input files
for i in ligand_name_list:
    gros.append('/home/mnz3v/mmtools/relativetools/Examples/checkrelativestates/gromacs_dg/gros/%s-0-eq.gro'%(i))
    originaltops.append('/home/mnz3v/mmtools/relativetools/Examples/checkrelativestates/gromacs_dg/originaltops/%s.top'%(i))
    originalitps.append('/home/mnz3v/mmtools/relativetools/Examples/checkrelativestates/gromacs_dg/originals/%s.itp'%(i))

for i in ligand_name_list:
    for j in ligand_name_list:
        if j > i:
            relativetops.append('/home/mnz3v/mmtools/relativetools/Examples/checkrelativestates/gromacs_dg/tops/%s-%s/%s.top'%(i, j, i))
            relativetops.append('/home/mnz3v/mmtools/relativetools/Examples/checkrelativestates/gromacs_dg/tops/%s-%s/%s.top'%(i, j, j))
            relativeitps.append('/home/mnz3v/mmtools/relativetools/Examples/checkrelativestates/gromacs_dg/relatives/%s-%s/%s.mut.itp'%(i, j, i))
            relativeitps.append('/home/mnz3v/mmtools/relativetools/Examples/checkrelativestates/gromacs_dg/relatives/%s-%s/%s.mut.itp'%(i, j, j))
            matches.append('/home/mnz3v/mmtools/relativetools/Examples/checkrelativestates/gromacs_dg/matches/%s-%s/%s-%s.mci'%(i, j, i, j))


# Pairwise matches by atom index
system = checkrelativestates(grompp, mdrun, g_energy, ligand_name_list, mdp, gros, originaltops, relativetops, originalitps = originalitps, relativeitps = relativeitps, matches = matches)

