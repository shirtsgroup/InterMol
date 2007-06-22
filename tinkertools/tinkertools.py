#!/usr/bin/python

#==================================================================================
# Tinker tools, written by D. Mobley, UCSF, 2007.
# Updates contributed by John D. Chodera, Stanford, 2007.
# 
# A variety of functions here; most are used by the main function, setup_calcs, as illustrated in the attached driver.py script.
#
# Currently, the attached driver.py script uses these tools to set up a sample calculation as follows
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
#
# Later, added additional tools to help with earlier setup:
# - name_to_mol2 to generate an initial mol2 file with conformer for a small molecule based on the IUPAC name.
# - mol2_to_xyz to generate an xyz file from a mol2 file and amoeba parameters, after prompting user for help identifying appropriate amoeba parameters for each atom.
# - solvate_molecule to solvate a molecule in desired size box of water with roughly right density
#
# Repository also contains several parameter sets
# - amoeba.prm, with a bugfix in the name of the oxygen on phenol (insignificant)
# - amoeba_fastwater.prm, where polarization/multipoles on water are replaced with fixed charges derived by J. Chodera to minimize uncertainty in transforming polarizable into fixed-charge water
#
# LIMITATIONS:
# - There is currently no separate constant pressure equilibration (or separate equilibration at all) at each lambda value. Equilibration would be done by just throwing out some of the production data. Leaving out separate constant pressure equilibration at each lambda, though, will introduce significant errors as the ligand sizes become significant relative to the simulation box size (because then the water densities won't be right). Ultimately we will need to add a separate constant rpessure equilibration step at each lambda.
#
# CHANGELOG:
# - 2007-06-21: JDC Heavily modified code to conform to Pande lab Python style guide.
# - 2007-06-01: DLM adding additional tools for doing several additional setup tasks:
#   (1): Generate initial conformation for specified small molecule using OpeneEye tools (added name_to_mol2 function to do this).
#   (2): Create an xyz file from an amoeba parameter file and mol2 file (with user input for appropriate amoeba atom types)
#   (3): Generate solvated box of suitable size using output from (2) (function: solvate_molecule)
#   (4): Run pre-equilibration of molecule using fixed-charge (nonpolarizable) version of AMOEBA water
#
# DEPENDENCIES:
# - OpenEye: OEChem, Omega, Lexichem (oeiupac)
# - mmtools: utilities.Units, utilities.Constants
#
# TODO
# - (JDC) What else uses generate_conf?  Can we streamline name_to_mol2?
# - Change "from X import *" to "import X" to reduce namespace clutter.
#==================================================================================


#===================================================================================
# IMPORTS
#===================================================================================
from numarray import *
import re
import random
import tempfile
from openeye.oechem import *
from openeye.oeomega import *
from openeye.oeiupac import *
import os
import commands
from math import *
import mmtools.utilities.Units as Units
import mmtools.utilities.Constants as Constants

#===================================================================================
# INITIALIZER
#===================================================================================

# Initialize random number generator.
random.seed()

#=================================================================================
# Tools for generating initial configuration of systems
#=================================================================================

def generate_conf(molecule, output_mol2_filename):
   """Generate a mol2 file from OEMol molecule, building conformation with Omega.

   ARGUMENTS
     molecule (OEMol) - the input molecule (coordinates need not be defined)
     output_mol2_filename (string) - filename of mol2 file to be written

   NOTES
     OpenEye's Omega is used to build a single conformation for the given molecule.

   EXAMPLES
     
   """

   # Generate a single new conformation using Omega.
   omega = OEOmega() # Initialize Omega.   
   omega.SetMaxConfs(1) # Generate a single conformation.
   omega.SetIncludeInput(False) # Omit the input conformation from the output.
   omega(molecule) # Run Omega to generate conformations for molecule.

   # Write the molecule to a mol2 file.
   output_molecule_stream = oemolostream()
   output_molecule_stream.open(output_mol2_filename)
   OEWriteMolecule(output_molecule_stream, molecule)
   output_molecule_stream.close()

   return

def name_to_mol2(IUPAC_name, output_mol2_filename):
   """Generate a mol2 file of a small molecule from its IUPAC name.

   ARGUMENTS
     IUPAC_name (string) - IUPAC name of molecule to generate
     output_mol2_filename (string) - filename of mol2 file to be written

   NOTES
     OpenEye LexiChem's OEParseIUPACName is used to generate the molecle, and Omega is used to generate a single conformation.

   EXAMPLES
     # Generate a mol2 file for phenol.
     name_to_mol2('phenol', 'phenol.mol2')
     
   """

   # Create an OEMol molecule from IUPAC name.
   molecule = OEMol() # create a molecule
   status = OEParseIUPACName(molecule, IUPAC_name) # populate the molecule from the IUPAC name
   OEAssignAromaticFlags(molecule) # check aromaticity.
   OEAddExplicitHydrogens(molecule) # add hydrogens
   molecule.SetTitle(IUPAC_name) # Set molecule title to IUPAC name.

   # Generate conformation with Omega and write to mol2 file.
   generate_conf(molecule, output_mol2_filename)

   return

def mol2_to_xyz(mol2file, amoebaprm, outxyz, amoebaname = None):
   """   
   Take an input mol2 file containing coordinates for some molecule and an amoeba parameter file;
   pick amoeba atom types out of key file;
   prompt user for input to identify amoeba atom types for appropriate atoms.
   Write out a tinker xyz file containing info from mol2 file.
   If 'amoebaname' is not specified, looks for amoeba parameters using the molecule name as read from the mol2 file. Specify this when the amoeba name for the molecule is different from that in the mol2 file.

   REQUIRED ARGUMENTS
     mol2file (string) - name of mol2 file containing coordinates
     amoebaprm (string) - name of AMOEBA parameter file containing atom parameters
     outxyz (string) - name of Tinker .xyz file to be written

   OPTIONAL ARGUMENTS
     amoebaname (string) - Name given the molecule in AMOEBA parameter file amoebaprm.  

   NOTES
     Requires that the molecule name as in the mol2 file match the name used for the molecule parameters in the tinker atom types section.

   LIMITATIONS
     Currently ignores case in molecular naming

   """

   # Create molecule.
   mol = OEMol()

   # Load molecule from mol2 file.
   input_molecule_stream=oemolistream()
   input_molecule_stream.open(mol2file)
   OEReadMolecule( input_molecule_stream, mol)
    
   # Determine AMOEBA molecule name (for amoebaprm file) from mol2 file title or, if specified, 'amoebaname' optional argument.
   if not amoebaname:
     name = mol.GetTitle()
   else: name = amoebaname
   print "\nMolecule: %s" % name
 
   # Open AMOEBA parameter file and store portion of atoms section dealing with this molecule; ignore case when looking through names.
   amoebaatoms = []
   file = open(amoebaprm,'r')
   for line in file.readlines():
      # May want to improve this; it is not particularly robust, for acetamide, for example, it also matches n-methylacetamide...
      if line.upper().find( name.upper() )>-1 and line.find('atom')==0:
         amoebaatoms.append(line)
   file.close()

   # Turn this into a dictionary by AMOEBA atom type key
   amoeba_key_to_type = {}
   # re to parse
   atoms = re.compile(r'atom\s+(?P<type>\d+)\s+\d+\s+\w+\s+"(?P<name>.*)".*')
   for line in amoebaatoms:
      m = atoms.match(line)
      if m:
         # Get atom key/name
         key = m.group('name')
         # Remove molecule name and following space
         key = key[ len(name)+1: ]
         amoeba_key_to_type[key] = int( m.group('type') )
      else: # Lines should match; if they don't, raise an exception
         print "Error: No match found in this line:", line
         raise "Atom line from parameter file does not match expected pattern."

   # Store amoeba types for atoms.
   amoeba_types = []
   for atom in mol.GetAtoms():
      print "\nPlease choose AMOEBA atom type for atom name %s with type %s, from the following options:" % (atom.GetName(), atom.GetType())
      keys = amoeba_key_to_type.keys()
      keys.sort()
      for key in keys:
         print "  %s: %s" % ( amoeba_key_to_type[key], key)
      type = int( raw_input() )
      if not amoeba_key_to_type.values().count(type)==1:
         print "Error: Please enter valid numeric type."
         type = int( raw_input() )
      if not amoeba_key_to_type.values().count(type)==1: raise TypeError
      amoeba_types.append(type)
      print "Using type: %s" % type
   
   # Determine number of atoms in molecule.
   numatoms = mol.NumAtoms()

   # Storage for xyz file.
   xyztext = []
   # Add header to xyz file
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

   # Loop over atoms in our molecule; print number, type and coordinates of each to xyz file, as well as bonds.
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

      # Store line
      xyztext.append( '%(numtxt)s   %(aname)s %(coordstext)s   %(atype)s   %(connections)s\n' % vars())

   # Write .xyz file.
   file = open(outxyz,'w')
   file.writelines(xyztext)
   file.close()

   return
  
def solvate_molecule(xyzfile, prmfile, outxyzfile, boxsize = 24.0 * Units.A):
   """
   Take an input xyz file;
   generate a box of solvent of specified size and correct density from tiling of water boxes;
   soak molecule into box;
   write new xyz file.

   ARGUMENTS
     xyzfile (string) - input solute coordinates
     prmfile (string) - Amoeba parameters
     outxyzfile (string) - Output xyz file name

   OPTIONAL ARGUMENTS
     boxsize (Units: length) - box edge length (default: 24.0 * Units.A)
     
   NOTES
     Uses gromacs tools genbox and trjconv for generating box.
     Uses spc water as the water with which to tile.
   """

   # Get absolute path to AMOEBA parameter file.
   prmfile = os.path.join(curdir,prmfile)

   # Store current directory.
   curdir = os.getcwd()

   # Make temporary working directory.
   tempdir = tempfile.mkdtemp() 

   # Change to temporary working directory.
   os.chdir(tempdir)
 
   # Generate temporay gro file with box size (in nm).
   grotext=['tmp gro file\n', ' 0\n', '  %.4f  %.4f  %.4f\n' % (boxsize / Units.nm, boxsize / Units.nm, boxsize / Units.nm)]
   file = open('box.gro','w')
   file.writelines(grotext)
   file.close()

   # Use genbox to create solvent in box.
   commands.getoutput('genbox -cp box.gro -cs spc216.gro -o solv.gro')
   
   # Convert solvated .gro file to PDB.
   commands.getoutput('echo "0" | trjconv -f solv.gro -s solv.gro -o solv.pdb')

   # Convert solvent PDB file to to Tinker .xyz using Tinker pdbxyz.
   file = open('in','w')
   file.write('solv.pdb\n%(prmfile)s\n' % vars())
   file.close()
   commands.getoutput("pdbxyz.x < in")
   # Now it is solv.xyz

   # Edit and change water types to use amoeba names/numbering.
   out = []
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

   # Copy solute here.
   os.system('cp %s solute.xyz' % os.path.join(curdir,xyzfile))

   # Soak solute into box using Tinker xyzedit.
   file = open('in','w')
   file.write('solute.xyz\n%(prmfile)s\n17\nsolv_names.xyz\n' % vars())
   file.close()
   commands.getoutput('xyzedit.x < in')

   # Move output back to working directory
   os.system('mv solute.xyz_2 %s' % (os.path.join(curdir,outxyzfile)) )

   # Restore working directory.
   os.chdir(curdir)

   # Clean up temporary directory.
   os.system('rm -r %(tempdir)s' % vars())

   return

#==================================================================================
# Tools for setting up parameter files
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

   #If there are redundant atoms in the list, only look them up once (make a new, non-redundant list)
   tmplist = []
   for elem in atomlist:
     if tmplist.count(elem)==0: tmplist.append(elem)
   atomlist = tmplist

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
        

def get_nonwater_types(xyzfile, waterlist=[22,23], include_redundant = False):
   """Read a tinker XYZ file and return a list of types of non-water atoms. 
INPUT:
- xyz file name
- optional: List of types of water atoms. If not provided, types are assumed to be 22 and 23 (as in amoeba.prm and in the TEMPLATE.key file used by this).
- optional: include_redundant: Provide a full list of atom types, including those which are used more than once. Default: False.
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
        if not include_redundant:
          if mol_atomtypes.count(typenum)==0: mol_atomtypes.append(typenum) 
        else:
          mol_atomtypes.append(typenum)   

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
# Setting up a calculation
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
  startdir = os.getcwd()  

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
  #Compute how often to output; get two frames out
  outfreq = 0.5*nsteps/1000
  runtext.append('dynamic.x mol %(nsteps)s 1.0 %(outfreq)s 4 298 1.0 > equil.log\n' % vars() )
  file=open('run_equilib.sh','w')
  file.writelines(runtext)
  file.close()
  #Run
  os.system('chmod 755 run_equilib.sh')
  os.system('./run_equilib.sh')
  #Back to original dir
  os.chdir(startdir)
  

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
    reprprm = os.path.join(  '../', str(lmb), name+str(lmb)+'.prm' )
    outtext+='cp %s %s.key\n' % (reprkey, name) 
    outtext+='cp %s %s%s.prm\n' % (reprprm, name, lmb) 
    outtext+= "printf '%s.arc\\nE\\n' > tmp\n" % name
    outtext+= "analyze.x < tmp > repr%s.log\n" % lmb
     
  return outtext

def setup_calc( keyfile, vackeyfile, targetdir, mol_atomtypes, amoebafile, xyzfile, dynfile, simlen = 100000, equillen=10000, cg_lambda = [0.0,0.106,0.226,0.368,0.553,0.8,1.0], vdw_lambda = [0.0,0.1,0.15,0.25,0.3,0.35,0.375,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0], wateratoms=[22,23]):
   """Setup a set of amoeba calculations in a target directory; will generate three subdirectories ('nochg_wat','nochg_vac','novdw_wat') within that directory, each containing all of the appropriate run input for calculations at a series of lambda values.
INPUT:
 - keyfile: Key file with settings for simulations
 - vackeyfile: Key file with settings for vacuum simulations
 - targetdir: Directory into which to put output/directories
 - mol_atomtypes: List of atom types in small molecule
 - amoebafile: Path to amoeba prm file (used for getting some LJ parameters for applying combination rules, etc); uses this file to generate modified amoeba prm files for the individual simulations, with modified interactions for the atoms in mol_atomtypes
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

     #Generate output parameter file there with scaled electrostatics
     scale_electrostatics( scalefactor, mol_atomtypes, amoebafile, os.path.join(opath, 'mol%s.prm' % lmb) ) 

     #Add random seed to key file, and parameters
     os.system('cp %s %s' % (keyfile, os.path.join( opath, 'mol.key' ) ) )
     file=open( os.path.join( opath, 'mol.key' ) , 'a')
     seed = random.randint(0, 5000)
     file.write('RANDOMSEED    %s\n' % seed)
     file.write('PARAMETERS mol%s.prm\n' % lmb)
     file.close()

     #Generate a copy of that file using the lambda value in the filename, as we'll need to do file shuffling later
     os.system('cp %s %s' % ( os.path.join(opath, 'mol.key'), os.path.join(opath, 'mol%s.key' % lmb) ) ) 

     #GENERATE RUN SCRIPT HERE; ALSO NEEDS TO HAVE REPROCESSING STUFF IN IT
     runtext='minimize.x mol 1.0 > mol_min.log\n'
     runtext+='mv mol.xyz_2 mol.xyz\n'
     #DO CONSTANT PRESSURE EQUILIBRATION
     runtext+='dynamic.x mol %s 1 1000 4 298 1  > equil.log\n' % equillen
     #DO PRODUCTION
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
     scale_electrostatics( scalefactor, mol_atomtypes, amoebafile, os.path.join(opath, 'mol%s.prm' % lmb) )

     #Add random seed to key file, and parameters
     os.system('cp %s %s' % (vackeyfile, os.path.join( opath, 'mol.key') ) )
     file=open( os.path.join( opath, 'mol.key' ) , 'a')
     seed = random.randint(0, 5000)
     file.write('RANDOMSEED    %s\n' % seed)
     file.write('PARAMETERS    mol%s.prm\n' % lmb)
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

     #Now edit the key file we set up for charging to add these parameters at the end
     cgfile=os.path.join(targetdir, 'nochg_wat', '1.0', 'mol1.0.prm') 
     file=open(cgfile,'r')
     outtext=file.readlines()
     for line in tmptext:
       outtext.append(line)
     #Write
     file=open( os.path.join(opath, 'mol%s.prm' % lmb), 'w')
     file.writelines(outtext)
     file.close()

     #Add random seed and parameter file name to key file
     os.system('cp %s %s' % (keyfile, os.path.join(opath, 'mol.key')))
     file = open( os.path.join(opath, 'mol.key'), 'r')
     text=file.readlines()
     text.insert(0, 'parameters mol%s.prm\n' % lmb)
     seed=random.randint(0,5000)
     text.insert(len(text)-1, 'RANDOMSEED   %s\n' % seed)
     file.close()
     file = open( os.path.join(opath, 'mol.key'), 'w')
     file.writelines(text)
     file.close()
  
     #Generate copy of key file tagged by lambda
     os.system('cp %s %s' % (os.path.join(opath, 'mol.key'), os.path.join(opath, 'mol%s.key' % lmb) ) ) 

     #GENERATE RUN SCRIPT HERE; ALSO NEEDS TO HAVE REPROCESSING STUFF IN IT
     runtext='minimize.x mol 1.0 > mol_min.log\n'
     runtext+='mv mol.xyz_2 mol.xyz\n'
     #DO CONSTANT PRESSURE EQUILIBRATION
     runtext+='dynamic.x mol %s 1 1000 4 298 1  > equil.log\n' % equillen
     #DO PRODUCTION
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
# Tools for analysis of Tinker trajectories.
#==================================================================================

def read_potential_energy(file):
  """Read the total potential energies of all frames reported in a Tinker .log file, returnning them as a numarray array of floats.

  ARGUMENTS
    file (string) - name of the Tinker .log file

  RETURNS
    potential_energies (numarray 1D float array) - array of potential energies, with length equal to number of frames

  """

  # Read lines from Tinker .log file.
  infile = open(file,'r')
  lines = infile.readlines()
  infile.close()

  # Construct Python list of potential energies.
  potential_energies = []
  for line in lines:
    if line.find('Total Potential Energy :')>-1:
       potential_energies.append(float(line.split()[4]))

  # Return numarray version of potential energy array.
  return array(potential_energies)
 

def get_volume(dynfile):
   """Get the box volume from a Tinker .dyn file.

   ARGUMENTS
     dynfile (string) - filename of Tinker .dyn file.

   RETURNS
     volume (Units: length^3) - box volume

   """

   print "dynfile: %s" % dynfile
   # Read the contents of the .dyn file.
   file = open(dynfile, 'r')
   lines = file.readlines()
   file.close()

   # Locate the line containing the box dimensions.
   line_index = 0
   for line in lines:
      if line.find('Periodic Box Dimensions') != -1:
         break
      line_index += 1
   # Get the line containing box dimensions.
   line = lines[line_index + 1]

   # Swap to scientific notation.
   line = line.replace('D','e')

   # Extract box dimensions.
   box_dimensions = line.split()

   # Compute volume.
   volume = 1.0
   for box_dimension in box_dimensions:
      volume *= float(box_dimension) * Units.A

   # Return volume.
   return volume

def get_number_of_solvent_molecules(xyzfile, solvent_atom_types = [ 22, 23, 23 ]):
   """Determine the number of solvent molecules in a Tinker .xyz file.

   ARGUMENTS
     xyzfile (string) - filename of the Tinker .xyz file

   OPTIONAL ARGUMENTS
     solvent_atom_types (list of integers) - atom types in solvent molecules (default: [22, 23, 23])

   RETURNS
     nsolvent (integer) - the number of solventmolecules in the Tinker .xyz file

   """

   # Read the .xyz file into memory.
   file = open(xyzfile,'r')
   lines = file.readlines()
   file.close()

   # Strip header line.
   lines.pop(0)

   # Count number of water molecules.
   nsolvent_atoms = 0   
   for line in lines:
     tmp = line.split()
     typenum = int(tmp[5])

     # If it is water, count it
     if solvent_atom_types.count(typenum)>0:
        nsolvent_atoms += 1

   # Determine number of molecules from number of atoms.
   nsolvent = nsolvent_atoms / len(solvent_atom_types)

   # Return number of solvent molecules in the .xyz file.
   return nsolvent

def LR_correction(basedir, dynfile = None, xyzfile = None, amoebaprm = None, cutoff = 9.0 * Units.A):
   """Compute the analytical long-range dispersion correction to the free energy of solvation for a solute for the Halgren 14-7 buffered potential from Tinker output.

   REQUIRED ARGUMENTS
     basedir (string) - path to base directory in which default filenames are found

   OPTIONAL ARGUMENTS
     dynfile (string) - path to Tinker .dyn file (from which box volume is extracted)
     xyzfile (string) - path to Tinker .xyz file (used to determine the number of waters and solute atom types and number)
     amoebaprm (string) - AMOEBA parameter file
     cutoff (Units: length) - nonbonded cutoff used in simulations

   RETURNS
     correction (Units: energy) - energy correction
     
   NOTES
     A sharp cutoff is assumed, and factors of near-unity are assumed to be unity (a safe assumption with default gamma and delta).

   """

   # Construct default names of files if not specified.
   if not xyzfile: xyzfile = os.path.join(basedir, 'mol.xyz')
   if not dynfile: dynfile = os.path.join(basedir, 'mol.dyn')
   if not amoebaprm: amoebaprm = os.path.join(basedir, 'mol0.0.prm')

   # TODO: Read gamma and delta corrections.
   
   # Determine the number of atoms in the simulation box from .xyz file.
   line = commands.getoutput('head -1 %(xyzfile)s' % vars())
   numatoms = int( line.split()[0] )

   # Determine box volume from .dyn file.
   boxvol = get_volume(dynfile)
   print "boxvol = %f A^3" % (boxvol / (Units.A**3))

   # Obtain list of all non-water atom types.
   solute_atomtypes = get_nonwater_types(xyzfile, include_redundant = True)

   # Determine number of water molecules.
   numwaters = get_number_of_solvent_molecules(xyzfile)
   print "numwaters = %d" % numwaters
   print "free volume / water = %f A^3" % ((boxvol / numwaters) / Units.A**3)
   print "density = %f g / cm^3" % ((numwaters / boxvol) * (18.0 * Units.amu) / (Units.g / Units.cm**3))

   # Read vdw parameters for solute and water types.
   (vdwparams, typedict) = read_LJ_params(amoebaprm, solute_atomtypes+[22,23])

   # List of water atoms, including redundancies.
   wateratoms = [22, 23, 23]

   # Define U(r) function for Halgren
   def U(r, eps, Rstar, gamma = 0.12, delta = 0.07):
      rho = r / Rstar
      return eps * ((1.+delta)/(rho+delta))**7 * ((1+gamma)/(rho**7 + gamma)-2)         
   # Define integrand I(r) = 4 pi r^2 rho U(r)
   def integrand(r, eps, Rstar, numwaters, boxvol):
      surface_area = 4. * pi * r**2
      number_density = numwaters / boxvol
      return surface_area * number_density * U(r, eps, Rstar)

   # Loop over solute atoms, then solvent atoms, and compute and accumulate eps_ij*Rij**7
   from scipy.integrate import quad
   from scipy.integrate import inf   

   correction = 0.0
   for solute_atom in solute_atomtypes:
     # Grab r and epsilon for solute atom
     Rstar_i = float(vdwparams[solute_atom]['r']) * Units.A
     eps_i = float(vdwparams[solute_atom]['eps']) * Units.kcal / Units.mol

     # Loop over water atoms.
     for water_atom in wateratoms:
       # Grab r and epsilon for water atoms.
       Rstar_j = float(vdwparams[water_atom]['r']) * Units.A
       eps_j = float(vdwparams[water_atom]['eps']) * Units.kcal / Units.mol
       
       # Use Halgren combination rules to compute effective pair r and epsilon.
       eps = 4. * (eps_i*eps_j) / (sqrt(eps_i) + sqrt(eps_j))**2
       Rstar = (Rstar_i**3 + Rstar_j**3) / (Rstar_i**2 + Rstar_j**2)

       # Numerically integrate
       density = numwaters / boxvol # water number density
       integral = quad(lambda x: integrand(x, eps, Rstar, numwaters, boxvol), cutoff, 10 * cutoff, epsabs = 1.0e-4 * Units.kcal/Units.mol)
       contribution = integral[0]
       
#       print "%d: r = %f A" % (solute_atom, Rstar / Units.A)
#       print "%d: eps = %f kcal/mol" % (solute_atom, eps / (Units.kcal/Units.mol))
#       print "%d:     %f kcal/mol" % (solute_atom, prefactor * eps * (Rstar**7) / (Units.kcal/Units.mol))
#       print "%d:     %f kcal/mol" % (solute_atom, contribution / (Units.kcal/Units.mol))
#       print ""

       # Accumulate contribution.
       #correction += -2. * pi * (numwaters/boxvol) * cutoff**(-4) * eps * (Rstar**7)
       correction += contribution

   return correction

def pV_correction(directory, pressure = 1.0 * Units.atm):
   """Compute the pressure-volume contribution the the hydration free energy.

   ARGUMENTS
     directory (string) - directory containing all lambda value directories

   OPTIONAL ARGUMENTS
     pressure (Units: pressure) - pressure (default: 1.0 * Units.atm)

   RETURNS
     pV_work (Units: energy) - pV work done by system for going from lambda = 0 to lambda = 1

   """

   # Determine directories for lambda = 0 and lambda = 1.
   lambda_0_dir = os.path.join(directory, '0.0')
   lambda_1_dir = os.path.join(directory, '1.0')

   # Determine .dyn file locations
   dynfile_0 = os.path.join(lambda_0_dir, 'mol.dyn')
   dynfile_1 = os.path.join(lambda_1_dir, 'mol.dyn')

   # Determine box volumes for each lambda value
   V_0 = get_volume(dynfile_0)
   V_1 = get_volume(dynfile_1)

   print "%e %e" % (V_0, V_1)
   # Compute work done by system in going from lambda = 0 to lambda = 1.
   work = pressure * (V_1 - V_0)

   return work

