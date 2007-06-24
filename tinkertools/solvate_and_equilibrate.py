#!/usr/bin/python

from tinkertools import *
import os

#Round 1
#names = [ 'phenol', 'methane', 'butane', 'methanol', 'ethanol', 'toluene', 'p-cresol', 'acetamide', 'methyl sulfide', 'methylamine', 'ammonia', 'dimethylamine', 'trimethylamine', 'pyrrolidine', 'N-methylacetamide', 'N,N-dimethylformamide', 'acetic acid', 'formaldehyde', 'acetaldehyde', 'hydrogen sulfide', 'dimethyl sulfide', 'dimethyl disulfide', 'methyl ethyl sulfide', 'benzene', 'ethylbenzene', '1H-imidazole'] 
#Round 2:
#names = ['ammonia', 'methylamine', 'dimethylamine', 'trimethylamine', 'pyrrolidine']
names = ['ethylamine']

toolsdir='/dmobley/svn/mmtools/tinkertools'

#number of copies to separately equilibrate
ncopies=3

runtext=['#!/usr/bin/python\n']
runtext.append('from tinkertools import *\n')
for name in names:
   print "Setting up %s..." % name
   ct = names.index(name)
   #Convert name to filename (remove spaces and commas)
   name = name.replace(' ','')
   name = name.replace(',','')
   
   #Solvate
   #solvate_molecule( os.path.join('molecules', name+'.xyz'), os.path.join(toolsdir, 'amoeba.prm'), os.path.join('molecules_solvated', name+'.xyz'))
   solvate_molecule( os.path.join('molecules_new', name+'.xyz'), os.path.join(toolsdir, 'amoeba.prm'), os.path.join('molecules_solvated', name+'.xyz'))
   #Set up equilibration
   for copy in range(ncopies):
     equildir = os.path.join('equilibrate',name+str(copy))
     if not os.path.exists(equildir): os.mkdir(equildir)
   
     #runtext=['#!/usr/bin/python\n']
     #runtext.append('from tinkertools import *\n')
     #Don't bother using fast water, as it's only 2x faster for constant pressure
     runtext.append('equilibrate( os.path.join("molecules_solvated", "%(name)s.xyz"), os.path.join("%(toolsdir)s", "amoeba.prm"),  "%(equildir)s", os.path.join( "%(toolsdir)s", "TEMPLATE_equilib.key"))\n' % vars())
     #file=open('run%s.py' % ct, 'w')
     #file.writelines(runtext)
     #file.close()
 
   #Just make one run file for now
   file=open('run.sh','w')
   file.writelines(runtext)
   file.close()

   #Submit to queue
   #os.system('queuesubmit.py -r run%s.sh -j %s' % (ct, 'l'+name))
