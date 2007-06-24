#!/usr/bin/python

from tinkertools import *

#Generate starting coordinates for AA sidechain analogs
names = [ 'methane', 'butane', 'methanol', 'ethanol', 'toluene', 'p-cresol', 'methyl sulfide', 'acetamide', 'ammonia', 'methylamine', 'ethyl amine', 'dimethylamine', 'trimethylamine', 'pyrrolidine', 'N-methylacetamide', 'N,N-dimethylformamide', 'acetic acid', 'formaldehyde', 'acetaldehyde', 'hydrogen sulfide', 'dimethyl sulfide', 'dimethyl disulfide', 'methyl ethyl sulfide', 'benzene', 'ethylbenzene', '1H-imidazole']
#Methyl sulfide and methanethiol are the same thing
#Also track tinker names if they are different
tinker_names={'methylamine':'methyl amine', 'ethylamine':'ethyl amine', 'dimethylamine': 'dimethyl amine', 'trimethylamine':'trimethyl amine', 'N-methylacetamide':'NMeAcetamide', 'N,N-dimethylformamide':'DiMeFormamide', 'dimethyl sulfide':'DiMe Sulfide', 'dimethyl disulfide':'DiMe Disulfide', 'methyl ethyl sulfide':'MeEt Sulfide', '1H-imidazole':'imidazole'} 
#Later add water
#Skip propyl amine for now as it's missing the nitrogen parameter


#Later: Add any other small molecules with known solvation free energies, etc
# - water
# - propanol, isopropanol (no hydration free energies)
# - gases? 
# - ammonia
# - whatever else in there has known solvation free energies (i.e. good series)
# - amine series (but watch out for propyl amine, which is missing the nitrogen)


#First generate mol2 files so I can download them and use them in matching the atom names
print "Generating mol2 files..."
for name in names:
  #Make output name be same as input name, but without any commas or spaces
  outname = name.replace(' ', '')
  outname = outname.replace(',','')
  name_to_mol2(name, 'molecules/%s.mol2' % outname)

print "Done generating mol2 files; press any key to continue and generate xyz files with user input..."
raw_input()
for name in names:
  outname = name.replace(' ', '')
  outname = outname.replace(',','')    

  #Check for different amoeba name
  if tinker_names.has_key(name): amoebaname = tinker_names[name]
  else: amoebaname = None
  mol2_to_xyz('molecules/%s.mol2' % outname, 'amoeba.prm', 'molecules/%s.xyz' % outname, amoebaname = amoebaname)
