#!/usr/bin/python

def add_lonepair( inxyz, outxyz, prmfile):
  """Edit an xyz file containing only a small molecule with an amine and requiring an AMOEBA "lone pair" to be added. Add a lone pair bonded to the correct nitrogen and at the right distance from that nitrogen, and allow the BAT parameters to take care of the rest in minimization. Has hardwired into it a list of candidate nitrogens (by atom class) and will look through the input xyz file for these to figure out which atom the lone pair should be bonded to."""
  
  #candidate atom classes for being connected to a lone pair
  candidate_classes = [32, 36, 36, 62]

  #Parse parameter file to obtain types of these classes
  file = open(prmfile,'r')
  candidate_types = []
  for line in file.readlines():
     if line.find('atom')==0:
        tmp=line.split()
        if candidate_classes.count( int(tmp[2]) )>0:
          candidate_types.append(int(tmp[1]))
  file.close()   


  #Open input file
  file = open(inxyz, 'r')
  text = file.readlines()
  file.close()

  #Obtain number of atoms
  numatoms = int(text[0].split()[0])
  #Figure out number of lone pair atom we'll add 
  LPnum = numatoms+1

  outtext=[]
  #Add first line
  outtext.append(text[0].replace(str(numatoms), str(LPnum)))
  
  found = False
  ct = 1 #Skip first line
  while not found:
    line = text[ct]
    tmp = line.split()
    type = int(tmp[5])
    if candidate_types.count(type) > 0 :
      found = True
      anum = tmp[0]
      coords = [float(tmp[2]), float(tmp[3]), float(tmp[4])]
    else: ct+=1

  #Edit line where it was found to add an additional connection
  text[ct] = text[ct].replace('\n', str(LPnum)+'\n')
   
  #Add these to output
  for ct in range(1,numatoms+1):
     outtext.append(text[ct])
  
  #Add new atom
  newcoords = coords
  newcoords[0]-=0.300 #Bond length is 0.3
  addtext = "%5s   LP %12.6f%12.6f%12.6f  46    %s\n" % (LPnum, coords[0], coords[1], coords[2], anum)
  outtext.append(addtext)
  file=open(outxyz,'w')
  file.writelines(outtext)
  file.close()


#Now add lonepairs to the molecules with lone pairs:
names = ['ammonia', 'methylamine', 'dimethylamine', 'trimethylamine', 'pyrrolidine', 'ethylamine']
for name in names:
  add_lonepair('molecules/%s.xyz' % name, 'molecules_new/%s.xyz' % name, '/dmobley/svn/mmtools/tinkertools/amoeba.prm')
