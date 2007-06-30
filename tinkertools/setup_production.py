#!/usr/bin/python

import os
from tinkertools import *

calcnum = 3 #Number of calculations per molecule
key = './TEMPLATE.key'
keyvac = './TEMPLATE_vac.key'

#Full set
#names = [ 'methane', 'butane', 'methanol', 'ethanol', 'toluene', 'p-cresol', 'methyl sulfide', 'acetamide', 'ammonia', 'methylamine', 'dimethylamine', 'trimethylamine', 'pyrrolidine', 'N-methylacetamide', 'N,N-dimethylformamide', 'acetic acid', 'formaldehyde', 'acetaldehyde', 'hydrogen sulfide', 'dimethyl sulfide', 'dimethyl disulfide', 'methyl ethyl sulfide', 'benzene', 'ethylbenzene', '1H-imidazole', 'ethylamine']
# Without amines
names = [ 'methane', 'butane', 'methanol', 'ethanol', 'toluene', 'p-cresol', 'methyl sulfide', 'acetamide', 'ammonia', 'N-methylacetamide', 'N,N-dimethylformamide', 'acetic acid', 'formaldehyde', 'acetaldehyde', 'hydrogen sulfide', 'dimethyl sulfide', 'dimethyl disulfide', 'methyl ethyl sulfide', 'benzene', 'ethylbenzene', '1H-imidazole']

#Convert names to names we actually use (remove commas, spaces)
tmp =[]
for name in names:
  name = name.replace(' ','')
  name = name.replace(',','')
  tmp.append(name)
names = tmp

#Get ready to set up calculations
ErrorLog=[] #Storage for errors (i.e. if equilibration didn't finish)

for name in names:
  print name
  for num in range(1,calcnum):
    print "  ", num
    equildir = os.path.join('equilibrate', name+str(num))
    targetdir = os.path.join('production', name+str(num))
    if not os.path.exists(targetdir): os.mkdir(targetdir)

    #Check that equilibration finished; if it did, do setup
    EquilError = "Error; unfinished equilibration."
    try:
      file = open( os.path.join(equildir, 'equil.log'), 'r')
      text=file.readlines()
      #Check equilibration finished
      if not text[-1].find('Induced Dipole File')>-1:
        ErrorLog.append('Unfinished calculation for %s; skipping...\n' % (name+str(num)))
        raise EquilError

      #If we are here, it finished

      #Obtain nonwater atom types for molecule by reading xyz file
      mol_atomtypes = get_nonwater_types( os.path.join(equildir, 'mol.xyz') )


      #Actually set up calculation
      setup_calc( './TEMPLATE.key', './TEMPLATE_vac.key', targetdir, mol_atomtypes, './amoeba.prm', os.path.join(equildir, 'mol.xyz'), os.path.join(equildir, 'mol.dyn') )


    except EquilError:
      print "Skipping %s..." % (name+str(num))

file = open('ErrorLog.txt','w')
file.writelines(ErrorLog)
file.close()

