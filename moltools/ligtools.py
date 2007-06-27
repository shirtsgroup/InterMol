import tempfile
import commands
from openeye.oechem import *
from openeye.oeiupac import *
import os

"""Ligtools:
- get_ligand: Takes a pdb file, extracts specified ligand and tries to protonate and assign bond types. Optionally writes out a output file (of type specified by the filename, i.e. mol2 or pdb) of it; also returns it as an OE mol.
- add_ligand_to_gro: Add a ligand at the end of an existing gro file.
- lig_mol2_to_gromacs: Convert a ligand mol2 file to gromacs top and gro files using antechamber and amb2gmx.pl. 

Requirements: 
- Amber and Antechamber installations in your path
- MMTOOLSPATH environment variable set to the location of mmtools (for amb2gmx.pl).

By D. Mobley, 6/28/2007."""

def get_ligand(pdbfile, resnum = None, resname = None, outfile = None):
   """Open specified pdb file (from provided path), and output the HETATM entries and associated connect entries for the specified ligand. Ligand can be specified by residue number, residue name, or both. Return an OE molecule containing the ligand. If outfile is provided, writes the molecule to the output file also."""

   file = open(pdbfile,'r')
   lines = file.readlines()
   file.close()

   ligatoms = []
   
   ofilename = tempfile.mktemp(suffix = '.pdb')
   ofile = open(ofilename, 'w')
   for line in lines:
      fieldtype = line[0:6].split()[0]
      if fieldtype == 'HETATM':
        tresname = line[17:20].split()[0]
        tresnum = int(line[22:26].split()[0])
        tatom = int(line[6:11].split()[0])
        if resnum and resname:
          if tresname == resname and tresnum == resnum:
            ofile.write(line)
            ligatoms.append(tatom)
        elif resnum:
          if tresnum == resnum: 
            ofile.write(line)
            ligatoms.append(tatom)
        elif resname:
          if tresname == resname: 
            ofile.write(line)
            ligatoms.append(tatom)
      #Also grab associated CONECT entries
      if fieldtype == 'CONECT':
         anums = line.replace('CONECT','').split()
         anums_i = [int(atom) for atom in anums] 
         found = False
         for atom in anums_i:
           if ligatoms.count(atom)>0 and not found:
              ofile.write(line)
              found=True

   ofile.close()

   molecule = OEMol()
   #Read molecule from pdb
   ifs = oemolistream(ofilename)
   OEReadMolecule(ifs, molecule)
   ifs.close()
   #Assign aromaticity/bonds, do some naming   
   OEAssignAromaticFlags(molecule) # check aromaticity.
   OEAddExplicitHydrogens(molecule) # add hydrogens   
   name = OECreateIUPACName(molecule) #Attempt to figure out IUPAC name for this; parts that cannot be recognized will be named 'BLAH'
   molecule.SetTitle(name) #Set title to IUPAC name

   if outfile:
      ostream = oemolostream()
      ostream.open(outfile)
      OEWriteMolecule(ostream, molecule)
      ostream.close()     

   #Don't forget to delete temp file
   commands.getoutput('rm %s' % ofilename)

   return molecule

def am1bcc_charge_mol2(ligmol2, outmol2, cleanup = True, judgetypes = None, nc = 0):
   """Use ANTECHAMBER to calculate AM1-BCC charges for specified mol2 and write output to specified mol2 file. Intermediate files are written to a temp dir which will be deleted unless cleanup is set to False. Judge bond/atom type options are default, unless a numerical argument specifying the -j type for antehcamber is provided under keyword judgetypes. Net charge assumed to be zero unless nc is specified. Output mol2 will have AMBER atom types."""
   workdir = tempfile.mkdtemp()
   tmpin = tempfile.mktemp( suffix = '.mol2', dir = workdir)
   commands.getoutput('cp %s %s' % (ligmol2, tmpin))
   dir = os.getcwd()
   os.chdir(workdir)
   tmpout = tempfile.mktemp( suffix ='.mol2', dir = workdir)

   if not judgetypes:
     print commands.getoutput('antechamber -i %s -fi mol2 -o %s -fo mol2 -c bcc -nc %s' % (tmpin, tmpout, nc))
   else:
     print commands.getoutput('antechamber -i %s -fi mol2 -o %s -fo mol2 -c bcc -nc %s -j %s' % (tmpin, tmpout, nc, judgetypes))
   commands.getoutput('cp %s %s' % (tmpout, os.path.join(dir, outmol2) ) ) 

   if cleanup:
     commands.getoutput('rm -r %s' % workdir)
   else:
     print "Work done in %s..." % workdir
   os.chdir(dir)

def ligmol2_to_gromacs(ligmol2, outname, cleanup = True, showWarnings = True):
   """Take a specified ligand mol2 file containing partial charges and AMBER atom types; convert it to gromacs topology and coordinate files. Uses antechamber (which must be in path), then parmchk and tleap which also must be in path. From this, generates prmtop and crd files. Then uses the amb2gmx.pl script to convert to GROMACS (hence requiring a full AMBER installation)."""
   workdir = tempfile.mkdtemp()
   commands.getoutput('cp %s %s' % (ligmol2, os.path.join(workdir, 'tmp.mol2')))
   dir = os.getcwd()
   os.chdir(workdir)
   
   #Check parameters
   print commands.getoutput('parmchk -i tmp.mol2 -f mol2 -o tmp.frcmod')
  
   leapscript = """source leaprc.gaff 
mods = loadAmberParams tmp.frcmod
ligand=loadMol2 tmp.mol2
desc ligand
check ligand
saveAmberParm ligand tmp.prmtop tmp.crd
quit"""
   file=open('leap.in', 'w')
   file.write(leapscript)
   file.close()

   tleapout = commands.getoutput('tleap -f leap.in')
   tleapout = tleapout.split('\n')
   if showWarnings:
     for line in tleapout:
       tmp = line.upper()
       if tmp.find('WARNING')>-1: print line
       if tmp.find('ERROR')>-1: print line

   amb2gmx = os.path.join(os.getenv('MMTOOLSPATH'), 'converters', 'amb2gmx.pl')
   commands.getoutput('%s --prmtop tmp.prmtop --crd tmp.crd --outname tmp' % amb2gmx)

   #Copy files back
   commands.getoutput('cp tmp.gro %s' % os.path.join(dir, outname+'.gro')) 
   commands.getoutput('cp tmp.top %s' % os.path.join(dir, outname+'.top')) 
   os.chdir(dir)

   #Cleanup
   if cleanup:
     commands.getoutput('rm -r %s' % workdir)
   else:
     print "Work done in %s..." % workdir

def add_ligand_to_gro(targetgro, liggro, outgro, resname = 'TMP'):
   """Take an existing gro file and a ligand gro file, and add the ligand at the end of the existing gro file. Makes no effort to adjust box size. Assumes ligand is a single residue; derives residue number from the last residue in the target gro file. Optionally specify ligand residue name, else TMP is used."""
   ligfile = open(liggro, 'r')
   liglines = ligfile.readlines()
   ligfile.close()
   targetfile = open(targetgro,'r')
   targetlines = targetfile.readlines()
   targetfile.close()

   ligatoms = int(liglines[1].split()[0])
   targetatoms = int(targetlines[1].split()[0])
   #Figure out existing number of residues
   lastres = targetlines[-2].split()[0]
   i = len(lastres)
   while not lastres[0:i].isdigit():
     i-=1
   lastresnum = int(lastres[0:i])
   #Residue number for ligand
   ligresnum = lastresnum+1
   
   #Compute new number of atoms
   newatomnum = ligatoms+targetatoms
   
   #Start creating output
   outtext = [ targetlines[0] ]
   outtext.append(' %s\n' % newatomnum)
   for line in targetlines[2:-1]:
      outtext.append(line)

   #Now append ligand
   resnumname='%4s%-4s' % (ligresnum, resname)
   for line in liglines[2:-1]:
      anum = int( line[15:20].split()[0] )
      newanum = targetatoms+anum
      line = ' '+resnumname+line[9:15]+('%5s' % newanum)+line[20:]
      outtext.append(line)
   #Add back on the box
   outtext.append(targetlines[-1])

   file=open(outgro, 'w')
   file.writelines(outtext)
   file.close()

def top_to_itp(topfile, outputitp, moleculetype = None):
   """Edit a topology file (usually a ligand topology file) to make it an itp file suitable for including in a topology. If optional argument moleculetype is specified, the name of the molecule in the moleculetype section will be replaced by the specified value. Assumes (for now) the topology file only contains one molecule definition."""
   file = open(topfile,'r')
   text = file.readlines()
   file.close()
  
   outtext=[]

   def read_through_section(text, linenum):
     linenum+=1
     line = text[linenum]
     while (line.find('[')==-1):
       linenum+=1
       try:
        line=text[linenum]
       except IndexError:
        return linenum-1
     return linenum

   linenum = 0
   while linenum < len(text):
     line = text[linenum]
     if line.find('defaults')>-1 and line.find(';')!=0:
       #Read through to next section
       linenum = read_through_section(text, linenum)
     elif line.find('system')>-1 and line.find(';')!=0:
       linenum = read_through_section(text, linenum)
     elif line.find('molecules')>-1 and line.find(';')!=0:
       linenum = read_through_section(text, linenum)
     elif moleculetype and line.find('moleculetype')>-1:
       outtext.append(line) 
       #Read through comments and append
       linenum+=1
       line = text[linenum]
       while line.find(';')==0:
          outtext.append(line)
          linenum+=1
          line = text[linenum]
       #Replace molecule name once we read through comments
       line = line.replace('solute', moleculetype)
       outtext.append(line)
     else:
       outtext.append(line)
       linenum +=1     
