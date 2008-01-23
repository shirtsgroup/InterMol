import tempfile
import commands
from openeye.oechem import *
from openeye.oeomega import *
from openeye.oeiupac import *
from openeye.oeshape import *
from openeye.oeproton import *
from openeye.oeiupac import *
import os

"""Ligtools:
- get_ligand: Takes a pdb file, extracts specified ligand and tries to protonate and assign bond types. Optionally writes out a output file (of type specified by the filename, i.e. mol2 or pdb) of it; also returns it as an OE mol.
- add_ligand_to_gro: Add a ligand at the end of an existing gro file.
- ligmol2_to_gromacs: Convert a ligand mol2 file to gromacs top and gro files using antechamber and amb2gmx.pl. 
- generate_conf_from_file: Generates (and optionally writes to file) a ligand conformation (or more than one) for a molecule file of arbitrary (OE readable) type; returns it.
- fit_mol_to_refmol: Fit a OE molecule (multi-conformer) to a reference molecule (single conformer, i.e. a ligand structure from a pdb file); write out an output file of the best N matches, where N is specified.
- EnumerateProtonation to enumerate possible protonation states.
- name_to_mol2: Generates mol2 file from IUPAC name using lexichem, Omega

Requirements: 
- Amber and Antechamber installations in your path
- MMTOOLSPATH environment variable set to the location of mmtools (for amb2gmx.pl).

By D. Mobley, 6/28/2007.
Revised by J. D. Chodera, 1/20/2008.
"""

#=============================================================================================
def extractMoleculeFromPDB(pdbfile, resnum = None, resname = None, chain = None, outfile = None):
   """Extract a ligand specified in the HETATM records of a PDB file.

   ARGUMENTS
     pdbfile (String) - the name of the PDB file from which the ligand is to be extracted

   RETURNS
     ligand (OEMol) - the ligand extracted from the PDB file
     
   OPTIONAL ARGUMENTS
     resnum - limit HETATM extraction to this residue, if specified (default: None)
     resname - limit HETATM extraction to this residue name, if specified (default: None)
     chain - limit HETATM extraction to this chain, if specified (default: None)
     outfile - if specified, the molecule is written to this output file (default: None)

   NOTES
     The molecule will be 'normalized' by protonating it and naming it according to its IUPAC name.
     Parts that are not recognized will be termed 'BLAH'.

   """

   # Read contents of source PDB file.
   file = open(pdbfile,'r')
   lines = file.readlines()
   file.close()

   ligatoms = [] # list of ligand atoms

   # Write temporary PDB file containing only HETATM and CONECT records from desired ligand.
   ofilename = tempfile.mktemp(suffix = '.pdb')
   ofile = open(ofilename, 'w')
   for line in lines:
      fieldtype = line[0:6].split()[0]
      if fieldtype == 'HETATM':
        tresname = line[17:20].split()[0]
        tresnum = int(line[22:26].split()[0])
        tatom = int(line[6:11].split()[0])
        tchain = line[21]
        match = True
        if chain:
          if tchain != chain:
            match = False
        if resnum:
          if tresnum != resnum:
            match = False
        if resname:
          if tresname != resname:
            match = False

        if match:
          ofile.write(line)
          ligatoms.append(tatom)
        
      # Also grab associated CONECT entries.
      if fieldtype == 'CONECT':
         anums = line.replace('CONECT','').split()
         anums_i = [int(atom) for atom in anums] 
         found = False
         for atom in anums_i:
           if ligatoms.count(atom)>0 and not found:
              ofile.write(line)
              found=True

   ofile.close()

   # Create a molecule from the temporary PDB file.
   molecule = OEMol()
   ifs = oemolistream(ofilename)
   OEReadMolecule(ifs, molecule)
   ifs.close()
   
   # Assign aromaticity/bonds, do naming.
   OEAssignAromaticFlags(molecule) # check aromaticity
   OEAddExplicitHydrogens(molecule) # add hydrogens   
   name = OECreateIUPACName(molecule) # attempt to determine IUPAC name
   molecule.SetTitle(name) # Set title to IUPAC name

   # Write molecule to file, if desired.
   if outfile:
      ostream = oemolostream()
      ostream.open(outfile)
      OEWriteMolecule(ostream, molecule)
      ostream.close()     

   # Delete the temporary PDB file.
   # TODO: Modify this to use os.path commands.
   commands.getoutput('rm %s' % ofilename)

   # Return the molecule.
   return molecule

# For backward-compatibility.
def get_ligand(pdbfile, resnum = None, resname = None, outfile = None, chain = None):
   return extractMoleculeFromPDB(pdbfile, resnum, resname, chain, outfile)
#=============================================================================================
def normalizeMolecule(molecule):
   """Normalize the molecule by checking aromaticity, adding explicit hydrogens, and renaming by IUPAC name.

   ARGUMENTS
     molecule (OEMol) - the molecule to be normalized.
   """
   
   # Assign aromaticity using Tripos model.
   # OEAssignAromaticFlags(molecule, OEAroModelTripos)
   OEAssignAromaticFlags(molecule, OEAroModelOpenEye)   

   # Assign Tripos atom types.
   # OETriposAtomTypes(molecule)

   # Assign Tripos atom names.
   # OETriposAtomTypeNames(molecule)

   # Add hydrogens.
   OEAddExplicitHydrogens(molecule)

   # Set title to IUPAC name.
   name = OECreateIUPACName(molecule)
   molecule.SetTitle(name)

   return molecule
#=============================================================================================
def readMolecule(filename, normalize = False):
   """Read in a molecule from a file (such as .mol2).

   ARGUMENTS
     filename (string) - the name of the file (such as .mol2)

   OPTIONAL ARGUMENTS
     normalize (boolean) - if True, molecule is normalized (renamed, aromaticity, protonated) after reading (default: False)

   RETURNS
     molecule (OEMol) - OEMol representation of molecule
   """

   # Open input stream.
   istream = oemolistream()
   istream.open(filename)

   # Create molecule.
   molecule = OEMol()   

   # Read the molecule.
   OEReadMolecule(istream, molecule)

   # Close the stream.
   istream.close()

   # Normalize if desired.
   if normalize: normalizeMolecule(molecule)

   return molecule   
#=============================================================================================
def writeMolecule(molecule, filename):
   """Write a molecule to a file (such as .mol2).

   ARGUMENTS
     molecule (OEMol) - the molecule to be written
     filename (string) - the file to write the molecule to (type autodetected from filename)

   RETURNS
     None

   NOTES
     If a .mol2 file is written, the substructure name is replaced with 'MOL'.
   """

   # Open output stream.
   ostream = oemolostream()
   ostream.open(filename)

   # Write molecule
   OEWriteMolecule(ostream, molecule)

   # Close the stream.
   ostream.close()

   # Replace substructure name if mol2 file.
   suffix = os.path.splitext(filename)[-1]
   if (suffix == '.mol2'):
      substructure_name = 'MOL'
      modifySubstructureName(filename, substructure_name)

   return
#=============================================================================================
def assignPartialCharges(molecule, charge_model = 'am1bcc'):
   """Assign partial charges to a molecule using OEChem oeproton.

   ARGUMENTS
     molecule (OEMol) - molecule for which charges are to be assigned

   OPTIONAL ARGUMENTS
     charge_model (string) - partial charge model (default: 'am1bcc') ['am1bcc']

   RETURNS
     charged_molecule (OEMol) - the charged molecule with GAFF atom types
   """

   # Check input pameters.
   supported_charge_models  = ['am1bcc']
   if not (charge_model in supported_charge_models):
      raise "Charge model %(charge_model)s not in supported set of %(supported_charge_models)s" % vars()

   charged_molecule = OEMol(molecule)   
   if charge_model == 'am1bcc':
      # assign partial charges using oeproton facility      
      OEAssignPartialCharges(charged_molecule, OECharges_AM1BCC)

   # Return the charged molecule
   return charged_molecule
#=============================================================================================
def assignPartialChargesWithAntechamber(molecule, charge_model = 'bcc', judgetypes = None, cleanup = True, verbose = False):
   """Assign partial charges to a molecule.

   ARGUMENTS
     molecule (OEMol) - molecule for which charges are to be computed

   OPTIONAL ARGUMENTS
     charge_model (string) - antechamber partial charge model (default: 'bcc')
     judgetypes (integer) - if specified, this is provided as a -j argument to antechamber (default: None)
     cleanup (boolean) - clean up temporary files (default: True)
     verbose (boolean) - if True, verbose output of subprograms is displayed

   RETURNS
     charged_molecule (OEMol) - the charged molecule with GAFF atom types

   REQUIREMENTS
     antechamber (on PATH)
   """

   # Create temporary working directory and move there.
   old_directory = os.getcwd()
   working_directory = tempfile.mkdtemp()   
   os.chdir(working_directory)

   # Write input mol2 file to temporary directory.
   uncharged_molecule_filename = tempfile.mktemp(suffix = '.mol2', dir = working_directory)
   print "Writing uncharged molecule to %(uncharged_molecule_filename)s" % vars()
   writeMolecule(molecule, uncharged_molecule_filename)

   # Create filename for output mol2 file.
   charged_molecule_filename = tempfile.mktemp(suffix ='.mol2', dir = working_directory)

   # Determine net charge of ligand from formal charges.
   formal_charge = formalCharge(molecule)

   # Run antechamber to assign GAFF atom types and charge ligand.
   command = 'antechamber -i %(uncharged_molecule_filename)s -fi mol2 -o %(charged_molecule_filename)s -fo mol2 -c %(charge_model)s -nc %(formal_charge)d' % vars()
   if judgetypes: command += ' -j %(judgetypes)d' % vars()
   if verbose: print command
   output = commands.getoutput(command)
   if verbose: print output

   # Read new mol2 file.
   print "Reading charged molecule from %(charged_molecule_filename)s" % vars()   
   charged_molecule = readMolecule(charged_molecule_filename)

   # Clean up temporary working directory.
   if cleanup:
      commands.getoutput('rm -r %s' % working_directory)
   else:
      print "Work done in %s..." % working_directory

   # Restore old working directory.
   os.chdir(old_directory)

   # Return the charged molecule
   return charged_molecule

# For backwards-compatibility.
def am1bcc_charge_mol2(ligmol2, outmol2, cleanup = True, judgetypes = None, nc = 0):
   """Deprecated. Use assignPartialCharges()."""
   assignPartialCharges(ligmol2, outmol2, 'am1bcc', cleanup = cleanup, judgetypes = judgetypes, nc = nc)
   return

#=============================================================================================
def ligmol2_to_gromacs(ligmol2, outname, cleanup = True, showWarnings = True, verbose = False):
   """Convert mol2 file with partial charges and AMBER atomtypes to gromacs topology/coordinate files.

   ARGUMENTS
     ligmol2 (string) - mol2 filename of ligand with partial charges and AMBER atomtypes
     outname (string) - prefix of filenames for gromacs topology/coordinate files to be produced

   OPTIONAL ARGUMENTS
     cleanup (boolean) - clean up temporary directories (default: True)
     showWarnings (boolean) - show any warnings during conversion process (default: True)
     verbose (boolean) - show complete output of tools (default: False)

   REQUIREMENTS
     antechamber (must be in PATH)
     amb2gmx.pl conversion script (must be in MMTOOLSPATH)
     AMBER installation (in PATH)

   """

   # Create temporary working directory and copy ligand mol2 file there.
   workdir = tempfile.mkdtemp()
   if verbose: print workdir
   commands.getoutput('cp %s %s' % (ligmol2, os.path.join(workdir, 'tmp.mol2')))
   dir = os.getcwd()
   os.chdir(workdir)
   
   # Generate frcmod file for additional GAFF parameters.
   print commands.getoutput('parmchk -i tmp.mol2 -f mol2 -o tmp.frcmod')

   # Create AMBER topology/coordinate files using LEaP.
   leapscript = """\
source leaprc.gaff 
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
   if verbose: print tleapout
   tleapout = tleapout.split('\n')
   # Shop any warnings.
   if showWarnings:
     print "Any warnings follow:"
     for line in tleapout:
       tmp = line.upper()
       if tmp.find('WARNING')>-1: print line
       if tmp.find('ERROR')>-1: print line

   # Use amb2gmx.pl to convert from AMBER to gromacs topology/coordinates.
   amb2gmx = os.path.join(os.getenv('MMTOOLSPATH'), 'converters', 'amb2gmx.pl')
   amb2gmx_output = commands.getoutput('%s --prmtop tmp.prmtop --crd tmp.crd --outname tmp' % amb2gmx)
   if verbose: print amb2gmx_output

   # Copy gromacs topology/coordinates to desired output files.
   commands.getoutput('cp tmp.gro %s' % os.path.join(dir, outname+'.gro')) 
   commands.getoutput('cp tmp.top %s' % os.path.join(dir, outname+'.top')) 
   os.chdir(dir)

   # Clean up temporary files.
   if cleanup:
     commands.getoutput('rm -r %s' % workdir)
   else:
     print "Work done in %s..." % workdir

   return
#=============================================================================================
def formalCharge(molecule):
   """Report the net formal charge of a molecule.

   ARGUMENTS
     molecule (OEMol) - the molecule whose formal charge is to be determined

   RETURN VALUES
     formal_charge (integer) - the net formal charge of the molecule
   """

   # Create a copy of the molecule.
   molecule_copy = OEMol(molecule)

   # Assign formal charges.
   OEFormalPartialCharges(molecule_copy)

   # Compute net formal charge.
   formal_charge = int(round(OENetCharge(molecule_copy)))

   # return formal charge
   return formal_charge
#=============================================================================================
def parameterizeForAmber(molecule, topology_filename, coordinate_filename, charge_model = False, judgetypes = None, cleanup = True, show_warnings = True, verbose = False, resname = None):
   """Parameterize small molecule with GAFF and write AMBER coordinate/topology files.

   ARGUMENTS
     molecule (OEMol) - molecule to parameterize (only the first configuration will be used if multiple are present)
     topology_filename (string) - name of AMBER topology file to be written
     coordinate_filename (string) - name of AMBER coordinate file to be written

   OPTIONAL ARGUMENTS
     charge_model (string) - if not False, antechamber is used to assign charges (default: False) -- if set to 'bcc', for example, AM1-BCC charges will be used
     judgetypes (integer) - if provided, argument passed to -j of antechamber to judge types (default: None)
     cleanup (boolean) - clean up temporary files (default: True)
     show_warnings (boolean) - show warnings during parameterization (default: True)
     verbose (boolean) - show all output from subprocesses (default: False)
     resname (string) - if set, residue name to use for parameterized molecule (default: None)
     
   REQUIREMENTS
     antechamber (must be in PATH)
     amb2gmx.pl conversion script (must be in MMTOOLSPATH)
     AMBER installation (in PATH)
   
   """

   # Create temporary working directory and copy ligand mol2 file there.
   working_directory = tempfile.mkdtemp()
   old_directory = os.getcwd()
   os.chdir(working_directory)
   if verbose: print "Working directory is %(working_directory)s"

   # Write molecule to mol2 file.
   tripos_mol2_filename = os.path.join(working_directory, 'tripos.mol2')
   writeMolecule(molecule, tripos_mol2_filename)

   if resname:
      # Set substructure name (which will become residue name) if desired.
      modifySubstructureName(tripos_mol2_filename, resname)
                 
   # Run antechamber to assign GAFF atom types.
   gaff_mol2_filename = os.path.join(working_directory, 'gaff.mol2')   
   command = 'antechamber -i %(tripos_mol2_filename)s -fi mol2 -o %(gaff_mol2_filename)s -fo mol2' % vars()
   if judgetypes: command += ' -j %(judgetypes)d' % vars()
   if charge_model:
      formal_charge = formalCharge(molecule)
      command += ' -c %(charge_model)s -nc %(formal_charge)d' % vars()   
   if verbose: print command
   output = commands.getoutput(command)
   if verbose: print output
   # TODO: SHOW WARNINGS, CHECK FOR ERRORS
   
   # Generate frcmod file for additional GAFF parameters.
   frcmod_filename = os.path.join(working_directory, 'gaff.frcmod')
   print commands.getoutput('parmchk -i %(gaff_mol2_filename)s -f mol2 -o %(frcmod_filename)s' % vars())

   # Create AMBER topology/coordinate files using LEaP.
   leapscript = """\
source leaprc.gaff 
mods = loadAmberParams %(frcmod_filename)s
molecule = loadMol2 %(gaff_mol2_filename)s
desc molecule
check molecule
saveAmberParm molecule amber.prmtop amber.crd
quit""" % vars()
   leapin_filename = os.path.join(working_directory, 'leap.in')
   outfile = open(leapin_filename, 'w')
   outfile.write(leapscript)
   outfile.close()

   tleapout = commands.getoutput('tleap -f %(leapin_filename)s' % vars())
   if verbose: print tleapout
   tleapout = tleapout.split('\n')
   # Shop any warnings.
   if show_warnings:
      print "Any LEaP warnings follow:"    # TODO: Only print this if there are any warnings to be printed.
      for line in tleapout:
         tmp = line.upper()
         if tmp.find('WARNING')>-1: print line
         if tmp.find('ERROR')>-1: print line

   # Restore old directory.
   os.chdir(old_directory)   

   # Copy gromacs topology/coordinates to desired output files.
   commands.getoutput('cp %s %s' % (os.path.join(working_directory, 'amber.crd'), coordinate_filename))
   commands.getoutput('cp %s %s' % (os.path.join(working_directory, 'amber.prmtop'), topology_filename))

   # Clean up temporary files.
   os.chdir(old_directory)
   if cleanup:
      commands.getoutput('rm -r %s' % working_directory)
   else:
      print "Work done in %s..." % working_directory

   return
#=============================================================================================
def parameterizeForGromacs(molecule, topology_filename, coordinate_filename, charge_model = False, cleanup = True, show_warnings = True, verbose = False, resname = None):
   """Parameterize small molecule with GAFF and write gromacs coordinate/topology files.

   ARGUMENTS
     ligmol2 (string) - mol2 filename of ligand with partial charges and AMBER atomtypes
     outname (string) - prefix of filenames for gromacs topology/coordinate files to be produced

   OPTIONAL ARGUMENTS
     charge_model (string) - if not False, antechamber is used to assign charges (default: False) -- if set to 'bcc', for example, AM1-BCC charges will be used
     cleanup (boolean) - clean up temporary directories (default: True)
     show_warnings (boolean) - show any warnings during conversion process (default: True)
     verbose (boolean) - show complete output of tools (default: False)
     resname (string) - if set, residue name to use for parameterized molecule (default: None)

   REQUIREMENTS
     antechamber (must be in PATH)
     amb2gmx.pl conversion script (must be in MMTOOLSPATH)
     AMBER installation (in PATH)

   """

   # Create temporary directory.
   working_directory = tempfile.mkdtemp()
   old_directory = os.getcwd()
   os.chdir(working_directory)
   if verbose: print "Working directory is %(working_directory)s"

   # Create AMBER coordinate/topology files.
   amber_topology_filename = os.path.join(working_directory, 'amber.prmtop')
   amber_coordinate_filename = os.path.join(working_directory, 'amber.crd')
   parameterizeForAmber(molecule, amber_topology_filename, amber_coordinate_filename, charge_model=charge_model, cleanup=cleanup, show_warnings=show_warnings, verbose=verbose, resname=resname)
   
   # Use amb2gmx.pl to convert from AMBER to gromacs topology/coordinates.
   amb2gmx = os.path.join(os.getenv('MMTOOLSPATH'), 'converters', 'amb2gmx.pl')
   command = '%(amb2gmx)s --prmtop %(amber_topology_filename)s --crd %(amber_coordinate_filename)s --outname gromacs' % vars()
   if verbose: print command
   amb2gmx_output = commands.getoutput(command)
   if verbose: print amb2gmx_output

   # Restore old directory.
   os.chdir(old_directory)   

   # Copy gromacs topology/coordinates to desired output files.
   commands.getoutput('cp %s %s' % (os.path.join(working_directory, 'gromacs.gro'), coordinate_filename))
   commands.getoutput('cp %s %s' % (os.path.join(working_directory, 'gromacs.top'), topology_filename))

   # Clean up temporary files.
   if cleanup:
      commands.getoutput('rm -r %s' % working_directory)
   else:
      print "Work done in %s..." % working_directory

   return
#=============================================================================================
def add_ligand_to_gro(targetgro, liggro, outgro, resname = 'TMP'):
   """Append ligand coordinates to the end of an existing gromacs .gro file.

   ARGUMENTS
     targetgro (string) - gromacs .gro file to which ligand coordinates are to be appended (not modified)
     liggro (string) - gromacs .gro file containing ligand coordinates (not modified)
     outgro (string) - filename of gromacs .gro file to contain merged target and ligand coordinates (created)

   OPTIONAL ARGUMENTS
     resname (string) - name for ligand residue (default: 'TMP')

   NOTES
     No effort is made to adjust box size.
     Ligand must be a single residue.
     New residue number is derived from last residue in target .gro file.
     
   """

   # Read ligand coordinates.
   ligfile = open(liggro, 'r')
   liglines = ligfile.readlines()
   ligfile.close()

   # Read target file coordinates.
   targetfile = open(targetgro,'r')
   targetlines = targetfile.readlines()
   targetfile.close()

   # Determine number of atoms in each file.
   ligatoms = int(liglines[1].split()[0])
   targetatoms = int(targetlines[1].split()[0])

   # Get number of last residue in target file.
   lastres = targetlines[-2].split()[0]
   i = len(lastres)
   while not lastres[0:i].isdigit():
     i-=1
   lastresnum = int(lastres[0:i])

   # Compute new residue number of ligand.
   ligresnum = lastresnum+1
   
   # Compute new number of atoms.
   newatomnum = ligatoms+targetatoms
   
   # Create new gromacs .gro file in memory.
   outtext = [ targetlines[0] ]
   outtext.append(' %s\n' % newatomnum)
   for line in targetlines[2:-1]:
      outtext.append(line)

   # Append the ligand coordinate lines, renumbering atom and residue numbers.
   resnumname='%4s%-4s' % (ligresnum, resname)
   for line in liglines[2:-1]:
      anum = int( line[15:20].split()[0] )
      newanum = targetatoms+anum
      line = ' '+resnumname+line[9:15]+('%5s' % newanum)+line[20:]
      outtext.append(line)

   # Add box line from target .gro file to end of file.
   outtext.append(targetlines[-1])

   # Write modified .gro file.
   file = open(outgro, 'w')
   file.writelines(outtext)
   file.close()

   return

#=============================================================================================
def top_to_itp(topfile, outputitp, moleculetype = None):
   """Transform a gromacs .top topology file into an .itp file suitable for inclusion.

   ARGUMENTS
     topfile (string) - the name of the gromacs .top topology file from which to draw topology
     outputitp (string) - the name of the .itp file to be created

   OPTIONAL ARGUMENTS
     moleculetype (string) - if specified, the molecule name in the [ moleculetype ] section will be renamed to this (default: None)

   NOTES
     Assumes (for now) the topology file only contains one molecule definition."""

   # Read the topology file.
   file = open(topfile,'r')
   text = file.readlines()
   file.close()
   
   # Create output .itp text.
   outtext=[]
   
   # Define a function to skip to next section.
   def read_through_section(text, linenum):
      """Skip to next section in text.
      
      ARGUMENTS
         text - list of lines of file
         linenum - line number to start reading at
         
      RETURNS
         linenum - line number at which new section was found
         """
      linenum+=1
      line = text[linenum]
      while (line.find('[')==-1):
         linenum+=1
         try:
            line=text[linenum]
         except IndexError:
            return linenum-1
      return linenum

   # Start reading at beginning of file.
   linenum = 0
   while linenum < len(text):
      line = text[linenum]
      if line.find('defaults')>-1 and line.find(';')!=0:
         # Read through to next section
         linenum = read_through_section(text, linenum)
      elif line.find('system')>-1 and line.find(';')!=0:
         linenum = read_through_section(text, linenum)
      elif line.find('molecules')>-1 and line.find(';')!=0:
         linenum = read_through_section(text, linenum)
         # This will be the end of file, so make sure we don't let a last line straggle in here.
         break;
      elif moleculetype and line.find('moleculetype')>-1:
         outtext.append(line) 
         # Read through comments and append.
         linenum+=1
         line = text[linenum]
         while line.find(';')==0:
            outtext.append(line)
            linenum+=1
            line = text[linenum]
         # Replace molecule name once we read through comments
         line = line.replace('solute', moleculetype)
         outtext.append(line)
         linenum+=1
      else:
         outtext.append(line)
         linenum +=1    

   # Write output .itp file.
   outfile = open(outputitp,'w')
   outfile.writelines(outtext)
   outfile.close()

   return
#=============================================================================================
def modifySubstructureName(mol2file, name):
   """Replace the substructure name (subst_name) in a mol2 file.

   ARGUMENTS
     mol2file (string) - name of the mol2 file to modify
     name (string) - new substructure name

   NOTES
     This is useful becuase the OpenEye tools leave this name set to <0>.
   """

   # Read mol2 file.
   file = open(mol2file, 'r')
   text = file.readlines()
   file.close()
   
   # Find the atom records.
   atomsec = []
   ct = 0
   while text[ct].find('<TRIPOS>ATOM')==-1:
     ct+=1
   ct+=1
   atomstart = ct
   while text[ct].find('<TRIPOS>BOND')==-1:
     ct+=1
   atomend = ct

   atomsec = text[atomstart:atomend]
   outtext=text[0:atomstart]
   repltext = atomsec[0].split()[7] # mol2 file uses space delimited, not fixed-width

   # Replace substructure name.
   for line in atomsec:
     # If we blindly search and replace, we'll tend to clobber stuff, as the subst_name might be "1" or something lame like that that will occur all over. 
     # If it only occurs once, just replace it.
     if line.count(repltext)==1:
       outtext.append( line.replace(repltext, name) )
     else:
       # Otherwise grab the string left and right of the subst_name and sandwich the new subst_name in between. This can probably be done easier in Python 2.5 with partition, but 2.4 is still used someplaces.
       # Loop through the line and tag locations of every non-space entry
       blockstart=[]
       ct=0
       c=' '
       for ct in range(len(line)):
         lastc = c
         c = line[ct]
         if lastc.isspace() and not c.isspace():
           blockstart.append(ct)
       line = line[0:blockstart[7]] + line[blockstart[7]:].replace(repltext, name, 1)
       outtext.append(line)
       
   # Append rest of file.
   for line in text[atomend:]:
     outtext.append(line)
     
   # Write out modified mol2 file, overwriting old one.
   file = open(mol2file,'w')
   file.writelines(outtext)
   file.close()

   return

#=============================================================================================
def generate_conf_from_file(infile, GenerateOutfile = False, outfile = None, maxconfs = 1, threshold=None, TorsionLib = None):
   """Use OE Omega to generate a conformation for the specified input file; return an OE mol containing the conformation(s). Optionally (if GenerateOutfile = True) write output to specified outfile as well. Input and output formats come from file names. 
Required arguments:
- Infile, with input file
Optional arguments:
- GenerateOutfile: Boolean specifying whether to write output file
- outfile: Name of outfile (format specified by suffix)
- maxconfs: Number of conformations to save. Default: 1.
- threshold: RMSD threshold for retaining conformers; lower thresholds retain more conformers. Default of None uses Omega's default threshold. Otherwise this should be a float.
- TorsionLib: Optionally specify a path to an Omega torsion library, which will be applied to perform an *additional* drive of torsions specified in that library after applying the default library. Default: None. 
"""
   #Open input file
   input_molecule_stream=oemolistream()
   input_molecule_stream.open(infile)

   #Initialize omega
   omega=OEOmega()

   #Adjst to desired number of conformerst
   if maxconfs:
     #Note that with my Omega version, although this "works", it doesn't actually control the maximum number of conformers if larger than 120; things top out at 120. Weird.
     omega.SetMaxConfs(maxconfs)
   else:
     omega.SetMaxConfs(1)
   #Don't include input in output
   omega.SetIncludeInput(False)

   #Adjust RMS threshold
   if threshold:
     omega.SetRMSThreshold(threshold) 
 
   #Create molecule
   molecule = OEMol()   

   OEReadMolecule(input_molecule_stream, molecule)
   name = OECreateIUPACName(molecule) #Attempt to figure out IUPAC name for this; parts that cannot be recognized will be named 'BLAH'
   molecule.SetTitle(name) #Set title to IUPAC name

   # Run Omega on the molecule to generate a set of reasonable conformations.
   omega(molecule)

   #If desired, do an additional torsion drive
   if TorsionLib:
     omega.SetTorsionLibrary(TorsionLib)
     omega(molecule)

   # Write the molecule to its own mol2 file.
   if GenerateOutfile:
    output_molecule_stream = oemolostream()
    output_molecule_stream.open(outfile)
    OEWriteMolecule(output_molecule_stream, molecule)
    output_molecule_stream.close()
   
   return molecule
#=============================================================================================
def fit_mol_to_refmol(refmol, fitmol_conformers, outfile, maxconfs = None):
   """Fit a multi conformer OE molecule (second argument) to a reference OE molecule using the OE Shape toolkit; write the resulting matches to specified outfile. Optionally specify the maximum number of conformers, 'maxconfs', to only get the best maxconfs matches written to that output file. Scores will be printed to stdout. Loosely based on OE Shape tookit documentation."""

   outfs = oemolostream(outfile)
   #Set up storage for overlay
   best = OEBestOverlay()
   #Set reference molecule
   best.SetRefMol(refmol)

   print "Ref. Title:", refmol.GetTitle()
   print "Fit Title:", fitmol_conformers.GetTitle()
   print "Num confs:", fitmol_conformers.NumConfs()

   resCount = 0

   #Each conformer-conformer pair generates multiple scores since there are multiple possible overlays; we only want the best. Load the best score for each conformer-conformer pair into an iterator and loop over it.
   scoreiter = OEBestOverlayScoreIter()
   OESortOverlayScores(scoreiter, best.Overlay(fitmol_conformers), OEHighestTanimoto())
   tanimotos = []
   for score in scoreiter:
      #Get the particular conformation of this match; transform it to overlay onto the reference structure
      outmol = OEGraphMol(fitmol_conformers.GetConf(OEHasConfIdx(score.fitconfidx)))
      score.Transform(outmol)

      #Write output
      OEWriteMolecule(outfs, outmol)

      print "FitConfIdx: %-4d" %score.fitconfidx,
      print "RefConfIdx: %-4d" % score.refconfidx,
      print "Tanimoto: %.2f" % score.tanimoto
      tanimotos.append(score.tanimoto)
      resCount +=1

      if resCount == maxconfs: break 

   return tanimotos   
#=============================================================================================
def fit_file_to_refmol(refmol, fitmol_file, outfile, maxconfs = None):
   """Fit molecules/conformers from a file (second argument) to a OE molecule (first argument). Like fit_mol_to_refmol, but will handle cases where the file contains, for example, multiple protonation states of the same molecule as well as multiple conformers. Writes the resulting matches to specified outfile. Optionally specify the maximum number of conformers, 'maxconfs', to only get the best maxconfs matches written to that output file. Scores will be printed to stdout. Loosely based on OE Shape tookit documentation. Attempts to recombine conformers from input file into multi-conformer molecules, then generate up to maxconfs conformers for each molecule."""

   infs = oemolistream()
   infs.open(fitmol_file)
   #Attempt to use OE tools to re-join multiple conformations into the same molecule when recombining, but keep molecules with different numbers of atoms/connectivities separate.
   infs.SetConfTest(OEAbsoluteConfTest())

   #Set up output
   outfs = oemolostream(outfile)
   #Set up storage for overlay
   best = OEBestOverlay()
   #Set reference molecule
   best.SetRefMol(refmol)


   print "Ref. Title:", refmol.GetTitle()
   #Now attempt to loop over input, reading a *molecule* at a time (molecules will now be multi-conformer, hopefully) and scoring separately for each molecule
   for fitmol_conformers in infs.GetOEMols():
     resCount = 0 #Track count of results for each molecule so we only keep this many conformers for each molecule

     print "Fit Title:", fitmol_conformers.GetTitle()
     print "Num confs:", fitmol_conformers.NumConfs()

     #Each conformer-conformer pair generates multiple scores since there are multiple possible overlays; we only want the best. Load the best score for each conformer-conformer pair into an iterator and loop over it.
     scoreiter = OEBestOverlayScoreIter()
     OESortOverlayScores(scoreiter, best.Overlay(fitmol_conformers), OEHighestTanimoto())
     for score in scoreiter:
        #Get the particular conformation of this match; transform it to overlay onto the reference structure
        outmol = OEGraphMol(fitmol_conformers.GetConf(OEHasConfIdx(score.fitconfidx)))
        score.Transform(outmol)

        #Write output
        OEWriteMolecule(outfs, outmol)

        print "FitConfIdx: %-4d" %score.fitconfidx,
        print "RefConfIdx: %-4d" % score.refconfidx,
        print "Tanimoto: %.2f" % score.tanimoto
        resCount +=1

        if resCount == maxconfs: break 

   return
#=============================================================================================
def split_mol2_poses(infile, outpath, outname):
  """Split a mol2 file into a separate file for each conformation.

  ARGUMENTS
    infile (string) - name of input file
    outpath (string) - path to output directory, which must already exist
    outname (string) - prefixed of output file names: will create output files [outname]0.mol2, [outname]1.mol2, ...

  OPTIONAL ARGUMENTS

  RETURNS

  """

  # Open input file
  file=open(infile,'r')

  outmol=[]
  ct=0
  for line in file.readlines():
    if line.find('MOLECULE')>-1:
       #Write output file every time we get to a new name MOLECULE field, but only if there is stuff here (i.e. we aren't at the first occurrence
       if len(outmol)>3:
         ofile=open(os.path.join(outpath,outname+str(ct)+'.mol2'),'w')
         ofile.writelines(outmol)
         ofile.close()
         outmol=[]
         ct+=1
    outmol.append(line)
  #Do last one too
  file.close()
  ofile=open(os.path.join(outpath,outname+str(ct)+'.mol2'),'w')
  ofile.writelines(outmol)
  ofile.close()

  #Return resulting number of files
  return ct+1

def EnumerateProtonation( mol, outfile, maxstates = 500 ):
   """Take a OE molecule; loop over conformations and enumerate likely protonation states; write output to resulting outfile. Default maxstates (maximum number of protonation states per conformer) of 500."""
   
   usearomaticity = True
   onlycountstates = False
   verbose = False
   #Set up output file
   ofile = oemolostream()
   ofile.open(outfile)
   #Initialize protonation typer
   typer = OETyperMolFunction( ofile, usearomaticity, onlycountstates, maxstates)

   #Suppress hydrogens so we get proper enumeration
   #TESTING: For c41, if we comment these out, it only enumerates protonations on the tetrazole and this ends up being good; otherwise it reassigns the bonds on the thiazole and never has a double bonded nitrogen. 
   #But does this cause problems for the others?
   #OESuppressHydrogens(mol)
   #OEAssignImplicitHydrogens(mol)
   #OEAssignFormalCharges(mol)
   for conf in mol.GetConfs():
     OEEnumerateFormalCharges(conf, typer, verbose)
     typer.Reset()
   
   ofile.close()

   return

#=============================================================================================
def createMoleculeFromIUPAC(name):
   """Generate a small molecule from its IUPAC name.

   ARGUMENTS
     IUPAC_name (string) - IUPAC name of molecule to generate

   RETURNS
     molecule (OEMol) - the molecule

   NOTES
     OpenEye LexiChem's OEParseIUPACName is used to generate the molecle.
     The molecule is normalized by adding hydrogens.
     Omega is used to generate a single conformation.

   EXAMPLES
     # Generate a mol2 file for phenol.
     molecule = createMoleculeFromIUPAC('phenol')
     
   """

   # Create an OEMol molecule from IUPAC name.
   molecule = OEMol() # create a molecule
   status = OEParseIUPACName(molecule, name) # populate the molecule from the IUPAC name

   # Normalize the molecule.
   normalizeMolecule(molecule)

   # Generate a conformation with Omega
   omega = OEOmega()
   omega.SetIncludeInput(False) # don't include input
   omega.SetMaxConfs(1) # set maximum number of conformations to 1
   omega(molecule) # generate conformation      

   # Return the molecule.
   return molecule

# Backwards-compatibility.
def name_to_mol2(IUPAC_name, output_mol2_filename):
   """Deprecated.  Use createMoleculeFromIUPAC()."""
   createMoleculeFromIUPAC(name = IUPAC_name, maxconfs = 1, output_filename = output_mol2_filename)
#=============================================================================================
# TEST DRIVER
if __name__ == '__main__':
   
   from mmtools.moltools.ligtools import *

   # Create a phenol molecule.
   molecule = createMoleculeFromIUPAC('phenol')

   # Write mol2 file for the molecule.
   writeMolecule(molecule, 'phenol.mol2')

   # Charge the molecule with OpenEye tools.
   charged_molecule = assignPartialCharges(molecule, charge_model = 'am1bcc')

   # Write charged molecule.
   writeMolecule(charged_molecule, 'phenol-am1bcc-openeye.mol2')

   # COMMENTED OUT BECAUSE ATOM TYPES GET JACKED UP
   # Charge the molecule with Antechamber.   
   # charged_molecule = assignPartialChargesWithAntechamber(molecule, charge_model = 'bcc')
   # Write charged molecule.
   # writeMolecule(charged_molecule, 'phenol-am1bcc-antechamber.mol2')
   
   # Write GAFF parameters for AMBER
   parameterizeForAmber(charged_molecule, topology_filename = 'amber.prmtop', coordinate_filename = 'amber.crd', resname = 'PHE')

   # Write GAFF parameters for gromacs.
   parameterizeForGromacs(charged_molecule, topology_filename = 'gromacs.top', coordinate_filename = 'gromacs.gro', resname = 'PHE')

   # Write GAFF parameters for gromacs, using antechamber to generate AM1-BCC charges.
   parameterizeForGromacs(molecule, topology_filename = 'gromacs-antechamber.top', coordinate_filename = 'gromacs-antechamber.gro', charge_model = 'bcc', resname = 'PHE')
   
