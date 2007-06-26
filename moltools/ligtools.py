import tempfile
import commands
from openeye.oechem import *
from openeye.oeiupac import *

"""Ligtools:
- get_ligand: Takes a pdb file, extracts specified ligand and tries to protonate and assign bond types. Optionally writes out a output file (of type specified by the filename, i.e. mol2 or pdb) of it; also returns it as an OE mol.

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

