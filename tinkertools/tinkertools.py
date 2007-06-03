#!/usr/bin/python

#==================================================================================
#Tinker tools, written by D. Mobley, 2007.
# 
#A variety of functions here; most are used by the main function, setup_calcs, as illustrated in the attached driver.py script.
#
#Currently, the attached driver.py script uses these tools to set up a sample calculation as follows
# - Begin from an xyz and dyn file that are assumed to be pre-equilibrated to the correct density, etc. in amoeba and already have the correct atom types (perhaps coming from pre-equilibration in another force field, or with a different water model, etc).
# - Figure out non-water atom types present by reading the xyz file.
# - Take a TEMPLATE key file associated with the repository, and the amoeba parameter file, and edit the TEMPLATE file to generate a key file for this molecule by copying over amoeba parameters for any non-water atoms; water parameters are assumed to already be in the TEMPLATE.key file. Do the same with a vacuum TEMPLATE file to generate a key file suitable for vacuum calculations. These are done by the function 'generate_parameter_file'
# - Set up the calculation as follows:
#   - Make output directory if it does not exist
#   - Generate directory structure for storing one copy each of charging, vdW and charging/vacuum calculations
#   - Set up charging in water calculations using the appropriate template key file from above and 
#     generating a set of key files with all of the electrostatics scaled back by a scale factor 
#     lambda; also generates run script for each lambda value that should also do reprocessing
#   - Set up charging in vacuum calculations using the same procedure but the vacuum template key file
#   - Set up vdW calculations beginning from the lambda=1.0 charging key file (that is, with 
#     electrostatics already turned off). Generate LJ parameters at each lambda using the combination
#     rules and our modified scheme for soft-core scaling using modification of gamma and delta for the       Halgren potentials. Create key files and run files.
#   - For each of these calculations above, put a different random number seed in the run key file at each lambda value.

# Later, added additional tools to help with earlier setup:
# - name_to_mol2 to generate an initial mol2 file with conformer for a small molecule based on the IUPAC name.
# - mol2_to_xyz to generate an xyz file from a mol2 file and amoeba parameters, after prompting user for help identifying appropriate amoeba parameters for each atom.
# - solvate_molecule to solvate a molecule in desired size box of water with roughly right density

# Repository also contains several parameter sets
# - amoeba.prm, with a bugfix in the name of the oxygen on phenol (insignificant)
# - amoeba_fastwater.prm, where polarization/multipoles on water are replaced with fixed charges derived by J. Chodera to minimize uncertainty in transforming polarizable into fixed-charge water

#LIMITATIONS:
# - There is currently no separate constant pressure equilibration (or separate equilibration at all) at each lambda value. Equilibration would be done by just throwing out some of the production data. Leaving out separate constant pressure equilibration at each lambda, though, will introduce significant errors as the ligand sizes become significant relative to the simulation box size (because then the water densities won't be right). Ultimately we will need to add a separate constant rpessure equilibration step at each lambda.

#CHANGELOG:
# - 6/1/2007: DLM adding additional tools for doing several additional setup tasks:
#   (1): Generate initial conformation for specified small molecule using OpeneEye tools (added name_to_mol2 function to do this).
#   (2): Create an xyz file from an amoeba parameter file and mol2 file (with user input for appropriate amoeba atom types)
#   (3): Generate solvated box of suitable size using output from (2) (function: solvate_molecule)
#   (4): Run pre-equilibration of molecule using fixed-charge (nonpolarizable) version of AMOEBA water

#DEPENDENCIES:
# - OpenEye: OEChem, Omega, Lexichem (oeiupac)
#==================================================================================


#===================================================================================
#IMPORTS
#===================================================================================
from numarray import *
import re
import random
import tempfile
from openeye.oechem import *
from openeye.oeomega import *
from openeye.oeiupac import *
import re
import os
import commands

#Seed random number generator
random.seed()
#=================================================================================
#Tools for generating initial configuration of systems
#=================================================================================

def generate_conf(mol,outmol2,name=None):
   """Takes an input OEMol molecule; uses openeye tools to generate an output mol2 file with coordinates for the atoms. Optionally, put the name specified in the name argumet into the name field of the resulting mol2 file."""
   #Initialize omega
   omega=OEOmega()
   #Only get one conformer out -- all we care about now is coordinates, not conformers
   omega.SetMaxConfs(1)
   #Don't include input in output
   omega.SetIncludeInput(False)
   
   # Run Omega on the molecule to generate a set of reasonable conformations.
   omega(mol)
   # Write the molecule to its own mol2 file.
   if name: mol.SetTitle(name)
   output_molecule_stream = oemolostream()
   output_molecule_stream.open(outmol2)
   OEWriteMolecule(output_molecule_stream, mol)
   output_molecule_stream.close()


def name_to_mol2(name, outmol2):
   """Use OEChem to convert the name of a molecule into an initial sybyl mol2 file with conformation generated by Omega. Write resulting molecule to output mol2 file."""
   mol = OEMol() #Create molecule
   status = OEParseIUPACName( mol, name )

   #Check aromaticity
   OEAssignAromaticFlags(mol)

   #Add hydrogens
   OEAddExplicitHydrogens(mol)

   #Generate conf with Omega and write
   generate_conf(mol, outmol2, name=name)

def mol2_to_xyz(mol2file, amoebaprm, outxyz, amoebaname = None):
   """Take an input mol2 file containing coordinates for some molecule and an amoeba parameter file; pick amoeba atom types out of key file; prompt user for input to identify amoeba atom types for appropriate atoms. Write out a tinker xyz file containing info from mol2 file. Requires that the molecule name as in the mol2 file match the name used for the molecule parameters in the tinker atom types section.
Limitations:
- Currently ignores case in molecular naming
Optional argument: Take amoebaname of molecule. Otherwise looks for amoeba parameters using the molecule name as read from the mol2 file. Specify this when the amoeba name for the molecule is different from that in the mol2 file."""

   #Create molecule
   mol = OEMol()

   #Load molecule from mol2 file
   input_molecule_stream=oemolistream()
   input_molecule_stream.open(mol2file)
   OEReadMolecule( input_molecule_stream, mol)
    
   #Molecule name
   if not amoebaname:
     name = mol.GetTitle()
   else: name = amoebaname
   print "\nMolecule: %s" % name
 
   #Open AMOEBA parameter file and store portion of atoms section dealing with this molecule; ignore case when looking through names
   amoebaatoms=[]
   file=open(amoebaprm,'r')
   for line in file.readlines():
     #May want to improve this; it is not particularly robust, for acetamide, for example, it also matches n-methylacetamide...
     if line.upper().find( name.upper() )>-1 and line.find('atom')==0:
       amoebaatoms.append(line)
   file.close()

   #Turn this into a dictionary by AMOEBA atom type key
   amoeba_key_to_type={}
   #re to parse
   atoms = re.compile(r'atom\s+(?P<type>\d+)\s+\d+\s+\w+\s+"(?P<name>.*)".*')
   for line in amoebaatoms:
     m = atoms.match(line)
     if m:
       #Get atom key/name
       key = m.group('name')
       #Remove molecule name and following space
       key = key[ len(name)+1: ]
       amoeba_key_to_type[key] = int( m.group('type') )
     else: #Lines should match; if they don't, raise an exception
       print "Error: No match found in this line:", line
       ParseError="Atom line from parameter file does not match expected pattern."
       raise ParseError

   #Store amoeba types for atoms
   amoeba_types=[]
   for atom in mol.GetAtoms():
     print "\nPlease choose AMOEBA atom type for atom %s with type %s, from the following options:" % (atom.GetType(), atom.GetName())
     keys = amoeba_key_to_type.keys()
     keys.sort()
     for key in keys:
       print "  %s: %s" % ( amoeba_key_to_type[key], key)
     type=int( raw_input() )
     if not amoeba_key_to_type.values().count(type)==1:
       print "Error: Please enter valid numeric type."
       type = int( raw_input() )
     if not amoeba_key_to_type.values().count(type)==1: raise TypeError
     amoeba_types.append(type)
   

        
   numatoms = mol.NumAtoms()

   #Storage for xyz file
   xyztext=[]
   #Add header to xyz file
   xyztext.append('  %s  %s - for AMOEBA, auto-generated from mol2\n' % (numatoms, name))

   #Build a dictionary storing which other atoms each atom is bonded to
   connected_atoms={}
   for bond in mol.GetBonds():
     start = bond.GetBgn().GetIdx()
     end = bond.GetEnd().GetIdx()
     if not connected_atoms.has_key(start): connected_atoms[start]=[]
     connected_atoms[start].append(end)
     if not connected_atoms.has_key(end): connected_atoms[end]=[]
     connected_atoms[end].append(start)

   #Loop over atoms in our molecule; print number, type and coordinates of each to xyz file, as well as bonds
   coords = mol.GetCoords()
   for atom in mol.GetAtoms():
      idx = atom.GetIdx()
      tinkernum = idx+1
      #Compute connected atoms
      connections = ''
      for num in connected_atoms[ idx ]:
        connections += str(num+1)+'    ' 
      #Coordinates
      coordstext = ''
      for c in coords[ idx ]:
         coordstext += "%12.6f" % around(c,6)
 
      atype = amoeba_types[ idx ]
      aname = atom.GetName()
      numtxt = "%5s" % tinkernum
      #Store line
      xyztext.append( '%(numtxt)s   %(aname)s %(coordstext)s   %(atype)s   %(connections)s\n' % vars())

   file=open(outxyz,'w')
   file.writelines(xyztext)
   file.close()
   
  
def solvate_molecule(xyzfile, prmfile, outxyzfile, boxsize = 24):
   """Take an input xyz file; generate a box of solvent of specified size and correct density from tiling of water boxes; soak molecule into box; write new xyz file. Uses gromacs tools genbox and trjconv for generating box. Optional argument: box size, in A. Default: 24. Uses spc water as the water with which to tile.
INPUT:
- xyzfile: Input solute coordinates
- prmfile: Amoeba parameters
- outxyzfile: Output xyz file name
- Optional boxsize, in angstroms.
"""

   boxsize=boxsize/10. #In NM for gromacs

   curdir=os.getcwd()
   #Make working directory
   tempdir = tempfile.mkdtemp() 

   #Copy amoeba file there from
   prmfile = os.path.join(curdir,prmfile)

   #Cd there
   os.chdir(tempdir)
 
   #Generate temporray gro file with box size in it
   grotext=['tmp gro file\n', ' 0\n', '  %.4f  %.4f  %.4f\n' % (boxsize, boxsize, boxsize)]
   file=open('box.gro','w')
   file.writelines(grotext)
   file.close()

   #Use genbox to solvate
   commands.getoutput('genbox -cp box.gro -cs spc216.gro -o solv.gro')
   
   #Convert to pdb
   commands.getoutput('echo "0" | trjconv -f solv.gro -s solv.gro -o solv.pdb')

   #Convert to xyz
   file=open('in','w')
   file.write('solv.pdb\n%(prmfile)s\n' % vars())
   file.close()
   commands.getoutput("pdbxyz.x < in")
   #Now it is solv.xyz

   #Edit and change water types to use amoeba names/numbering
   out=[]
   file=open('solv.xyz','r')
   for line in file.readlines():
     if line.find('OW')>-1: line = line.replace(' 0 ', ' 22 ')
     elif line.find('HW')>-1: 
        line=line.replace('HW1','HW ')
        line=line.replace('HW2','HW ')
        line = line.replace(' 0 ', ' 23 ')
     out.append(line)
   file.close()
   file=open('solv_names.xyz','w')
   file.writelines(out)
   file.close()

   #Copy solute here
   os.system('cp %s solute.xyz' % os.path.join(curdir,xyzfile))

   #Soak solute into box
   file=open('in','w')
   file.write('solute.xyz\n%(prmfile)s\n17\nsolv_names.xyz\n' % vars())
   file.close()
   commands.getoutput('xyzedit.x < in')

   #Move output  back to working directory
   os.system('mv solute.xyz_2 %s' % (os.path.join(curdir,outxyzfile)) )
   os.chdir(curdir)

   #Clean up
   os.system('rm -r %(tempdir)s' % vars())

#==================================================================================
#Tools for setting up parameter files
#==================================================================================

def scale_electrostatics(scalefactor,atomlist,infile,outfile):
  """Function to take an atom list and scale back the electrostatic parameters by "scalefactor" for all of those atoms as listed in the amoeba parameter file "infile" and write it to a new "outfile". Makes no attempt to preserve existing spacing in parameter file, as the spacing seems not to matter.
Input:
- scalefactor: The factor to scale the charges by
- atomlist: List of atom types for which to scale back the charges by this factor
- infile: Input parameter file
- outfile: Output modified parameter file
"""

  #Read infile
  file=open(infile,'r')
  text=file.readlines()
  file.close()
  
  #Regexes for parsing sections
  multipolehdr=re.compile('(?P<begin>multipole\s+\t*\s+)(?P<a1>\d+)\s+(?P<a2>\d+)\s+(?P<a3>\d+)\s+(?P<chg>[+-]*\d+\.\d+).*')
  polarizeline=re.compile('(?P<begin>polarize\s+\t*)(?P<atom>\d+)\s+(?P<polzn>[+-]*\d+\.\d+)(?P<rest>.*)')  

  outarray=[] 
  ct=0
  while ct<len(text):
    line=text[ct]
    #Look for multipole and polarization
    m=multipolehdr.match(line)
    p=polarizeline.match(line)
    if m:
      #Look to see if any of the atoms are in our list
      fnd=False
      for atom in atomlist:
        if ( int(m.group('a1'))==atom or int(m.group('a2'))==atom or int(m.group('a3'))==atom):
            fnd=True
      if fnd:
        #If so, then modify the charges in this block
        outarray.append(m.group('begin')+'     '+m.group('a1')+'     '+m.group('a2')+'     '+m.group('a3')+'              '+str(float(m.group('chg'))*scalefactor)+'\n')
        #Other lines in block
        for i in range(1,5):
          tmp=text[ct+i].split()
          scaled=''
          for elem in tmp:
            scaled+=str(float(elem)*scalefactor)+'    '
          outarray.append('                                       '+scaled+'\n')
      else:
        outarray.append(line)
        for i in range(1,5): outarray.append(text[ct+i])
      ct+=5
    #If polarization match
    elif p:
      fnd=False
      for atom in atomlist:
        if int(p.group('atom'))==atom:
           fnd=True
      if fnd:
        #If matches any of our atoms of interest, reduce the polarization accordingly
        outarray.append(p.group('begin')+p.group('atom')+'               '+str(float(p.group('polzn'))*scalefactor)+p.group('rest')+'\n')
      else: outarray.append(line)
      ct+=1
    else:
       outarray.append(line)
       ct+=1

  #Save output file
  file=open(outfile,'w')
  file.writelines(outarray)
  file.close()


def read_LJ_params(file,atomlist):
   """Reads a tinker parameter file to obtain the sigma and epsilon parameters for the atom TYPES (not classes!) specified in atomlist; returns them in a dictionary by atom type number. Also returns, in case desired, a dictionary mapping the atomlist (types) to classes as a second return argument. 
INPUT:
- file: A tinker parameter file
- atomlist: A list of atom types to obtain sigma/epsilon parameters for
Returns:
- vdwdict[atom]['r','eps']: Dictionary by atom type containing r and epsilon values for atoms
- typedict[atom]: What is the class for the type 'atom'?

"""
   file=open(file,'r')
   text=file.readlines()
   file.close()

   #Sort list to make sure atoms are consecutive to avoid having to loop over file multiple times
   atomlist.sort()
   vdwdict={}
   typedict={}

   #Parse file looking up atom classes lines
   ct=0
   #Atom classes occur earlier in the file, so look these up first
   for atom in atomlist:
     fnd=False
     while not fnd:
       line=text[ct]
       if line.find('atom')==0:
         tmp=line.split()
         if tmp[1]==str(atom):
           fnd=True
           #Store atom class under atom type
           typedict[atom]=tmp[2]
         else: ct+=1
       else:ct+=1

   #Now we have the type dictionary, so go on to look for vdw information
   #Now sort atom list by class
   def class_compare(x,y):
     if int(typedict[x])>int(typedict[y]): return 1
     elif int(typedict[x])<int(typedict[y]): return -1
     else: return 0
   atomlist.sort(class_compare)
     
   #Look for atoms in order
   for atom in atomlist:
     #Look for LJ parameters by atom class, store under atom type
     fnd=False
     while not fnd:
        line=text[ct]
        if line.find('vdw')==0:
          tmp=line.split()
          if tmp[1]==typedict[atom]:
            fnd=True
            #store params
            vdwdict[atom]={}
            vdwdict[atom]['r']=tmp[2]
            vdwdict[atom]['eps']=tmp[3]
          else: ct+=1
        else: ct+=1
   return vdwdict,typedict
        

def get_nonwater_types(xyzfile, waterlist=[22,23]):
   """Read a tinker XYZ file and return a list of types of non-water atoms. 
INPUT:
- xyz file name
- optional: List of types of water atoms. If not provided, types are assumed to be 22 and 23 (as in amoeba.prm and in the TEMPLATE.key file used by this).
RETURNS:
- atomlist: List of types of atoms.
Works by parsing through the file and checking all the atom types; adding any to the list that are not the water atom types. 
"""

   file=open(xyzfile,'r')
   text=file.readlines()
   file.close()
   text=text[1:] #Strip header
   #Storage
   mol_atomtypes=[]

   for line in text:
     tmp=line.split()
     typenum=int(tmp[5])
     water=False
     #Check to make sure it isn't water
     if waterlist.count(typenum)>0: water=True
     if not water:
        if mol_atomtypes.count(typenum)==0: mol_atomtypes.append(typenum) 
   
   return mol_atomtypes

def generate_parameter_file( template, amoebaprm, nonwatertypes, outfile ):
   """Opens a specified template.key file and amoeba parameter file; generates a new key file which it writes to an output file. This new key file additionally contains the parameters for the nonwater atoms as copied from the amoeba parameter file.
INPUT:
- template: Path name of template key file
- amoebaprm: Path name of amoeba parameter file
- nonwatertypes: List of nonwater types with parameters to copy (note: Does not copy excess dihedrals/angles, etc., i.e. angles involving only one of these atoms).
- outfile: Path name of output file.
OUTPUT:
- output file written to outfile
"""
   
   #Get template text
   file=open(template,'r')
   text=file.readlines()
   file.close()

   file=open(amoebaprm, 'r')
   amotext=file.readlines()
   file.close()

   #Storage for output
   addtext=[]
   
   #Obtain vdw classes
   (dum, type_to_class) = read_LJ_params(amoebaprm, nonwatertypes)
   vdwclasses = type_to_class.values()


   #Count line numbers
   ct=0
   #Loop over lines; copy relevant entries over from the following sections: 'atom', 'vdw', 'bond', 'angle', 'strbnd', 'opbend', 'torsion', 'pitors', 'multipole', 'polarize'. Currently just uses a bunch of if statements to do this.
   while ct < len(amotext):
     line=amotext[ct]
     tmp=line.split()

     #Skip blank lines
     if len(tmp) > 0:

       #'atom' entries
       if tmp[0] == 'atom':
         type=int(tmp[1])
         #Store any atom entries matching atom types in our small molecule to the output
         if nonwatertypes.count(type)>0:
           addtext.append(line)

       #'vdw' entries
       elif tmp[0]=='vdw':
        vdwtype=tmp[1]
        if vdwclasses.count(vdwtype)>0:
          addtext.append(line) 

       #'bond' entries
       elif tmp[0]=='bond':
        vdwtypes=[tmp[1],tmp[2]]
        counts=0
        for type in vdwtypes: 
          if vdwclasses.count(type)>0: counts+=1
        if counts>1:
          addtext.append(line)     

       #'angle' entries
       elif tmp[0] == 'angle':
         vdwtypes=[tmp[1],tmp[2],tmp[3]]
         counts=0
         for type in vdwtypes:
          if vdwclasses.count(type)>0: counts+=1
         if counts>2:
           addtext.append(line)

       #'strbnd' entries
       elif tmp[0] == 'strbnd':
         vdwtype=tmp[1]
         if vdwclasses.count(vdwtype)>0:
           addtext.append(line)
       
       #'opbend' entries
       elif tmp[0] == 'opbend':
         vdwtypes=[tmp[1],tmp[2]]
         counts=0
         for type in vdwtypes:
           if vdwclasses.count(type)>0: counts+=1
           if counts>1:
             addtext.append(line)

       #'torsion' entries
       elif tmp[0] == 'torsion':
         vdwtypes=[tmp[1],tmp[2],tmp[3],tmp[4]]
         counts=0
         for type in vdwtypes:
           if vdwclasses.count(type)>0: counts+=1
         if counts>3:
          addtext.append(line)

       #'pitors' entries
       elif tmp[0] == 'pitors':
         vdwtypes=[tmp[1],tmp[2]]
         counts=0
         for type in vdwtypes:
          if vdwclasses.count(type)>0: counts+=1
         if counts>0:
          addtext.append(line)

       #'multipole' entries, by atom type:
       elif tmp[0]=='multipole':
         atypes=[int(tmp[1]),int(tmp[2]),int(tmp[3])]
         counts=0
         for atype in atypes:
           if nonwatertypes.count(atype)>0: counts+=1
         if counts>2:
           for idx in range(0,5):
             addtext.append(amotext[ct+idx])
         ct+=4 #Increment counter since we took multiple lines

       #'polarize' entries, by atom type:
       elif tmp[0]=='polarize':
         atype=int(tmp[1])
         if nonwatertypes.count(atype)>0:
           addtext.append(line)

     #increment line num
     ct+=1

   file=open(outfile,'w')
   file.writelines(text)
   file.writelines(addtext)
   file.close()

#==================================================================================
#Setting up a calculation
#==================================================================================

def equilibrate( xyzfile, prmfile, workdir, template, boxsize=24, nsteps=50000):
  """Do minimization and equilibration for specified xyz and parameter files in specified working directory. Adds random seed to key file.
INPUT:
- xyzfile for tinker
- prmfile with corresponding atom types
- workdir: directory in which to work (will put xyz and prm files there)
- template: A template key file to which will be added key words relating to pressure regulation, box size, etc.
OPTIONAL INPUT:
- boxsize, as used in setting up xyz file. Used in determining a-axis setting (which is set using box size plus 5%). Default: 24A.
- nsteps: Number of timesteps of constant pressure equilibration to run (units 1 fs). Default: 50000.
  """
  
  if not os.path.exists(workdir): os.mkdir(workdir) 
  os.system('cp %s %s' % (xyzfile, os.path.join(workdir,'mol.xyz') ))
  os.system('cp %s %s' % (prmfile, os.path.join(workdir,'prmfile') ))
   
  file=open(template,'r')
  text=file.readlines()
  file.close()
  
  #Add parameter file
  text.insert(0,'parameters\tprmfile\n')
  #Compute and add a-axis value
  size=boxsize*1.05
  text.insert(5,'a-axis\t %.3f\n' % size)
  
  #Add random seed to key file
  seed = random.randint(0, 5000)
  text.insert(7, 'RANDOMSEED    %s\n' % seed)


  #Store key file
  file=open( os.path.join(workdir,'mol.key'), 'w')
  file.writelines(text)
  file.close()

  #Switch to working directory
  os.chdir(workdir)
  #Make run script
  runtext=['minimize.x mol 1.0 > mol_min.log\n']
  runtext.append('mv mol.xyz_2 mol.xyz\n')
  runtext.append('dynamic.x mol %(nsteps)s 1.0 1.0 4 298 1.0 > equil.log\n' % vars() )
  file=open('run_equilib.sh','w')
  file.writelines(runtext)
  file.close()
  #Run
  os.system('chmod 755 run_equilib.sh')
  os.system('./run_equilib.sh')

  

def stripatoms(atomlist,filein, fileout):
  """Strip all atom types listed in atomlist from the xyz file 'filein' and write a new xyz file 'fileout' that contains only the remaining atoms."""
  file=open(filein,'r')
  text=file.readlines()
  file.close()
  
  outtext=[]
  outtext.append(text[0])
  text=text[1:]
  for line in text:
     atom=line.split()[5]
     save=True
     for at in atomlist:
       if atom==str(at):
         save=False
     if save:
       outtext.append(line)

  #Fix number of atoms in first line
  oldnum=outtext[0].split()[0]
  newnum=len(outtext)-1
  outtext[0]=outtext[0].replace(str(oldnum),str(newnum))

  file=open(fileout,'w')
  file.writelines(outtext)


def get_reprocessing_text( name, reprlambdas ):
  """Return text to shuffle key files around and reprocess existing trajectories for a particular calculation with key files from other lambda values.
INPUT:
- Name of key/trajectory files to reprocess
- reprlambdas: List of lambda values at which to reprocess. Will look for key files for these in ../(lambda)/name(lamda).key for each lambda in this list.
OUTPUT:
- Returns a text string containing commands to do this reprocessing (i.e. with a shell script)."""

  outtext=''
  for lmb in reprlambdas:
    reprkey = os.path.join(  '../', str(lmb), name+str(lmb)+'.key' )
    outtext+='cp %s %s.key\n' % (reprkey, name) 
    outtext+= "printf '%s.arc\\nE\\n' > tmp\n" % name
    outtext+= "analyze.x < tmp > repr%s.log\n" % lmb
     
  return outtext

def setup_calc( keyfile, keyfile_vac, targetdir, mol_atomtypes, amoebafile, xyzfile, dynfile, simlen = 100000, cg_lambda = [0.0,0.106,0.226,0.368,0.553,0.8,1.0], vdw_lambda = [0.0,0.1,0.15,0.25,0.3,0.35,0.375,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0], wateratoms=[22,23]):
   """Setup a set of amoeba calculations in a target directory; will generate three subdirectories ('nochg_wat','nochg_vac','novdw_wat') within that directory, each containing all of the appropriate run input for calculations at a series of lambda values.
INPUT:
 - keyfile: Path of a key file containing all relevant water and small molecule parameters. Note that key file should contain an a-axis setting, but this can be large as this should be read from the dyn file.
 - keyfile_vac: Equivalent key file but for vacuum calculations (i.e. no PME)
 - targetdir: Directory into which to put output/directories
 - mol_atomtypes: List of atom types in small molecule
 - amoebafile: Path to amoeba prm file (used for getting some LJ parameters for applying combination rules, etc)
 - xyzfile: xyz file to use for calculations
 - dynfile: dyn file to use for calculations
OPTIONAL INPUT:
  - simulation length
  - cg_lambda: Array of lambda values for charging
  - vdw_lambda: Array of lambda values for LJ calcs
  - wateratoms: Atom types of water atoms; default 22 and 23
OUTPUT: 
  - Set of directories/scripts for running calculations, with appropriate files. Note that all output xyz/key/dyn files will be named "mol.xxx" for consistency.
OTHER NOTES:
- Calculations are currently done with only minimization and then constant volume production. Some amount of time from production should be discarded as equiliibration when analyzing. Also, this approach will fail for large molecules, where the system will need to be equilibrated separately at constant pressure for each lambda value in order to maintain correct densities.
"""
   #Settings for soft core potentials
   alpha=0.88
   lpow=1.0


   #Check target dir exists; if not, make it
   if not os.path.isdir(targetdir): os.mkdir(targetdir)
   
   #Obtain some vdw parameters we need
   (vdwparams, type_to_class) = read_LJ_params(amoebafile, mol_atomtypes)
   (watvdwparams, dum) = read_LJ_params(amoebafile, wateratoms)

   #Charging calculations
   cgdir=os.path.join(targetdir, 'nochg_wat')
   if not os.path.exists(cgdir): os.mkdir(cgdir)

   print 'Setting up charging in water...'
   for lmb in cg_lambda:
     opath=os.path.join(cgdir, str(lmb))
     if not os.path.exists(opath): os.mkdir(opath)
     #Copy xyz and dyn files there
     os.system('cp %s %s' % (xyzfile, os.path.join(opath, 'mol.xyz')))     
     os.system('cp %s %s' % (dynfile, os.path.join(opath, 'mol.dyn')))     

     #Factor to scale polarization/charges by
     scalefactor = (1. - lmb)

     #Generate output key file there with scaled electrostatics
     scale_electrostatics( scalefactor, mol_atomtypes, keyfile, os.path.join(opath, 'mol.key') ) 

     #Add random seed to key file
     file=open( os.path.join( opath, 'mol.key' ) , 'a')
     seed = random.randint(0, 5000)
     file.write('RANDOMSEED    %s\n' % seed)
     file.close()

     #Generate a copy of that file using the lambda value in the filename, as we'll need to do file shuffling later
     os.system('cp %s %s' % ( os.path.join(opath, 'mol.key'), os.path.join(opath, 'mol%s.key' % lmb) ) ) 

     #GENERATE RUN SCRIPT HERE; ALSO NEEDS TO HAVE REPROCESSING STUFF IN IT
     runtext='minimize.x mol 1.0 > mol_min.log\n'
     runtext+='mv mol.xyz_2 mol.xyz\n'
     #NOTE HERE WE ARE NOT DOING CONSTANT PRESSURE EQUILIBRATION FIRST
     runtext+='dynamic.x mol %s 1 0.1 2 298 1  > mol.log\n' % simlen
     #reprocessing -- get text for reprocessing at neighboring lambda values and self.
     idx = cg_lambda.index(lmb)
     nbrs=[lmb]
     if idx+1 < len(cg_lambda): nbrs.append(cg_lambda[idx+1])
     if idx-1 >=0 : nbrs.append(cg_lambda[idx-1])
     reprtext = get_reprocessing_text( 'mol', nbrs ) 

     #Write run script
     file=open( os.path.join(opath, 'run.sh'), 'w')
     file.write(runtext)
     file.write(reprtext)
     file.close()

   #Charging in vacuum
   print 'Setting up charging in vacuum...'
   cgdir=os.path.join(targetdir, 'nochg_vac')
   if not os.path.exists(cgdir): os.mkdir(cgdir)

   for lmb in cg_lambda:
     opath=os.path.join(cgdir, str(lmb))
     if not os.path.exists(opath): os.mkdir(opath)
     #Copy xyz file there after stripping it of waters
     stripatoms( wateratoms, xyzfile, os.path.join(opath, 'mol.xyz'))
      
     #Factor to scale polarization/charges by
     scalefactor = (1. - lmb)

     #Generate output key file there with scaled electrostatics
     scale_electrostatics( scalefactor, mol_atomtypes, keyfile_vac, os.path.join(opath, 'mol.key') )

     #Add random seed to key file
     file=open( os.path.join( opath, 'mol.key' ) , 'a')
     seed = random.randint(0, 5000)
     file.write('RANDOMSEED    %s\n' % seed)
     file.close()

     #Generate a copy of that file using the lambda value in the filename, as we'll need to do file shuffling later
     os.system('cp %s %s' % ( os.path.join(opath, 'mol.key'), os.path.join(opath, 'mol%s.key' % lmb) ) )

     #GENERATE RUN SCRIPT HERE; ALSO NEEDS TO HAVE REPROCESSING STUFF IN IT
     runtext='minimize.x mol 1.0 > mol_min.log\n'
     runtext+='mv mol.xyz_2 mol.xyz\n'
     runtext+='dynamic.x mol %s 1 0.1 2 298 1  > mol.log\n' % simlen
     #reprocessing -- get text for reprocessing at neighboring lambda values and self.
     nbrs=[lmb]
     if idx+1 < len(cg_lambda): nbrs.append(cg_lambda[idx+1])
     if idx-1 >=0 : nbrs.append(cg_lambda[idx-1])
     reprtext = get_reprocessing_text( 'mol', nbrs )

     #Write run script
     file=open( os.path.join(opath, 'run.sh'), 'w')
     file.write(runtext)
     file.write(reprtext)
     file.close()

   #Set up LJ calculation
   print 'Setting up vdW calculation...'
   vdwdir=os.path.join(targetdir, 'novdw_wat')
   if not os.path.exists(vdwdir): os.mkdir(vdwdir)
   for lmb in vdw_lambda:
     opath=os.path.join(vdwdir, str(lmb))
     if not os.path.exists(opath): os.mkdir(opath)
     #Copy xyz and dyn files there
     os.system('cp %s %s' % (xyzfile, os.path.join(opath, 'mol.xyz')))
     os.system('cp %s %s' % (dynfile, os.path.join(opath, 'mol.dyn')))
   
     #Compute gamma and delta for halgren potentials
     ghal = 0.12+alpha*lmb**lpow
     dhal = 0.07+alpha*lmb**lpow

     #STorage for parameters
     tmptext=[]

     #Use combination rules to compute modified r and epsilon for each pair of interactions between
     #water and other
     for w in wateratoms:
              #Get r and eps for water atom type
      watr=float(watvdwparams[w]['r'])
      wateps=float(watvdwparams[w]['eps'])
      for at in mol_atomtypes:
         #Retrieve vdw class of this atom type
         aclass=type_to_class[at]
         #Retrieve vdw parameters, which are indexed by atom TYPE.
         atomr=float(vdwparams[at]['r'])
         atomeps=float(vdwparams[at]['eps'])

         #Use combination rules to compute combined r and epsilon
         comb_eps=4*(atomeps*wateps)/(sqrt(atomeps)+sqrt(wateps))**2
         comb_r=(atomr**3+watr**3)/(atomr**2+watr**2)
         #Compute scaled epsilon
         eps_new=comb_eps*(1.-lmb)**4


         #Store data to output text array to insert in key file later
         tmptext.append('vdwpr        %s    %s       %.4f     %.5f\n' % (w,aclass,comb_r,eps_new))
         tmptext.append('halpr        %s    %s       %.4f     %.4f\n' % (w,aclass,ghal,dhal))

     #Now edit the key file we set up for charging to add these parameters at the end; also use new random seed
     cgfile=os.path.join(targetdir, 'nochg_wat', '1.0', 'mol1.0.key') 
     file=open(cgfile,'r')
     outtext=[]
     for line in file.readlines():
       if line.find('RANDOMSEED')>-1:
          seed=random.randint(0,5000)
          line='RANDOMSEED    %s\n' % seed
       outtext.append(line)
     for line in tmptext:
       outtext.append(line)
   
     #Write
     file=open( os.path.join(opath, 'mol.key'), 'w')
     file.writelines(outtext)
     file.close()

     #Generate a copy of that file using the lambda value in the filename, as we'll need to do file shuffling later
     os.system('cp %s %s' % ( os.path.join(opath, 'mol.key'), os.path.join(opath, 'mol%s.key' % lmb) ) )
     #GENERATE RUN SCRIPT HERE; ALSO NEEDS TO HAVE REPROCESSING STUFF IN IT
     runtext='minimize.x mol 1.0 > mol_min.log\n'
     runtext+='mv mol.xyz_2 mol.xyz\n'
     #NOTE HERE WE ARE NOT DOING CONSTANT PRESSURE EQUILIBRATION FIRST
     runtext+='dynamic.x mol %s 1 0.1 2 298 1  > mol.log\n' % simlen
     #reprocessing -- get text for reprocessing at neighboring lambda values and self.
     idx = vdw_lambda.index(lmb)
     nbrs=[lmb]
     if idx+1 < len(vdw_lambda): nbrs.append(vdw_lambda[idx+1])
     if idx-1 >=0 : nbrs.append(vdw_lambda[idx-1])
     reprtext = get_reprocessing_text( 'mol', nbrs )
     #Write run script
     file=open( os.path.join(opath, 'run.sh'), 'w')
     file.write(runtext)
     file.write(reprtext)
     file.close()


#==================================================================================
#Tools for analysis
#==================================================================================

def read_potential_energy(file):
  """Read total potential energy from Tinker output (stdout, stored to a log file) and return it as a (numarray) array of floats with length equal to the number of frames."""
  infile=open(file,'r')
  intext=infile.readlines()
  infile.close()
  poten=[]
  for line in intext:
    if line.find('Total Potential Energy :')>-1:
       poten.append(float(line.split()[4]))
  #Convert to numarray
  return array(poten)
  
