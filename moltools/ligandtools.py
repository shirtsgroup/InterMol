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
# METHODS FOR READING, EXTRACTING, OR CREATING MOLECULES
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
# METHODS FOR INTERROGATING MOLECULES
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
# METHODS FOR MODIFYING MOLECULES
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

   WARNING
     This module is currently broken, as atom names get all jacked up during readMolecule() for these mol2 files.
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
   # TODO: Atom names get all jacked up here -- is there a way to fix this?

   # Clean up temporary working directory.
   if cleanup:
      commands.getoutput('rm -r %s' % working_directory)
   else:
      print "Work done in %s..." % working_directory

   # Restore old working directory.
   os.chdir(old_directory)

   # Return the charged molecule
   return charged_molecule
#=============================================================================================
# METHODS FOR WRITING OR EXPORTING MOLECULES
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
# METHODS FOR MANIPULATING GROMACS TOPOLOGY AND COORDINATE FILES
#=============================================================================================
def perturbGromacsTopology(topology_filename, molecule, perturb_torsions = True, perturb_vdw = True, perturb_charges = True):
   """Modify a gromacs topology file to add perturbed-state parameters.

   ARGUMENTS
     topology_file (string) - the name of the topology file to modify
     molecule (OEMol) - molecule corresponding to contents of topology file

   OPTIONAL ARGUMENTS
     perturb_torsions (boolean) - if True, torsions whose central bond is not in an aromatic ring will be turned off in B state (default: True)
     perturb_vdw (boolean) - if True, van der Waals interactions will be turned off in B state (default: True)
     perturb_charges (boolean) - if True, charges will be turned off in B state (default: True)

   NOTES
     This code currently only handles the special format gromacs topology files produced by amb2gmx.pl -- there are allowed variations in format that are not treated here.
   """

   def stripcomments(line):
      """Return line with whitespace and comments stripped.
      """
      # strip comments
      index = line.find(';')
      if index > -1:
         line = line[0:index]
      # strip whitespace
      line = line.strip()
      # return stripped line
      return line         

   def extract_section(lines, section):
      """Identify lines associate with a section.

      ARGUMENTS
        lines (list of strings) - the lines in the file
        section (string) - the section name to locate

      RETURNS
        indices (list of integers) - line indices within lines belonging to section
      """

      indices = list()

      nlines = len(lines)
      for start_index in range(nlines):
         # get line
         line = stripcomments(lines[start_index])
         # split into elements
         elements = line.split()
         # see if keyword is matched
         if (len(elements) == 3):
            if (elements[0]=='[') and (elements[1]==section) and (elements[2]==']'):
               # increment counter to start of section data and abort search
               start_index += 1
               break

      # throw an exception if section not found
      if (start_index == nlines):
         raise "Section %(section)s not found." % vars()

      # Locate end of section.
      for end_index in range(start_index, nlines):
         # get line
         line = stripcomments(lines[end_index])
         # split into elements
         elements = line.split()
         # see if keyword is matched
         if (len(elements) == 3):
            if (elements[0]=='['):
               break
      
      # compute indices of lines in section
      indices = range(start_index, end_index)

      # return these indices
      return indices

   # Read the contents of the topology file.
   infile = open(topology_filename, 'r')
   lines = infile.readlines()
   infile.close()

   # Parse atomtypes.
   atomtypes = list() # storage for atom types
   indices = extract_section(lines, 'atomtypes')
   for index in indices:
      # extract the line
      line = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip line if number of elements is less than expected
      if (nelements < 7): continue
      # parse elements
      atomtype = dict()
      atomtype['name'] = elements[0]
      atomtype['bond_type'] = elements[1]
      atomtype['mass'] = float(elements[2])
      atomtype['charge'] = float(elements[3])
      atomtype['ptype'] = elements[4]
      atomtype['sigma'] = float(elements[5])
      atomtype['epsilon'] = float(elements[6])
      # append
      atomtypes.append(atomtype)

   # augment the 'atomtypes' list with perturbed atom types,    
   if perturb_vdw:
      indices = extract_section(lines, 'atomtypes')         
      for atomtype in atomtypes:
         # make perturbed record
         perturbed_atomtype = atomtype
         perturbed_atomtype['name'] += '_pert'
         perturbed_atomtype['epsilon'] = 0.0
         # form the line
         line = "%(name)-10s%(bond_type)6s      0.0000  0.0000  A %(sigma)13.5e%(epsilon)13.5e ; perturbed\n" % perturbed_atomtype
         # insert the new line
         lines.insert(indices[-1], line)
         indices.append(indices[-1]+1)

   # Process [ atoms ] section
   atoms = list()
   atom_indices = dict()
   indices = extract_section(lines, 'atoms')
   for index in indices:
      # extract the line
      line = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 8): continue
      # parse line
      atom = dict()
      atom['nr'] = int(elements[0])
      atom['type'] = elements[1]
      atom['resnr'] = int(elements[2])
      atom['residue'] = elements[3]
      atom['atom'] = elements[4]
      atom['cgnr'] = int(elements[5])
      atom['charge'] = float(elements[6])
      atom['mass'] = float(elements[7])
      
      # set perturbation type
      atom['typeB'] = atom['type']
      if perturb_vdw: atom['typeB'] += '_pert' # perturbed vdw type
      atom['chargeB'] = atom['charge']
      if perturb_charges: atom['chargeB'] = 0.0 # perturbed charges
      
      # construct a new line
      line = "%(nr)6d %(type)10s %(resnr)6d %(residue)6s %(atom)6s %(cgnr)6d %(charge)10.5f %(mass)10f %(typeB)6s %(chargeB)10.5f ; perturbed\n" % atom
      
      # replace the line
      lines[index] = line

      # store atoms
      atoms.append(atom)

      # store lookup of atom atom names -> atom numbers
      atom_indices[atom['atom']] = atom['nr']
      
   # Set rotatable bond torsions in B state to zero, if desired.
   if perturb_torsions:
      # Determine list of rotatable bonds to perturb.
      rotatable_bonds = list()
      for bond in molecule.GetBonds():
         # DEBUG
         print "%d %d" % (bond.GetBgnIdx(), bond.GetEndIdx())
         # This test is used because bond.IsRotor() doesn't seem to work correctly (e.g. phenol).
         if (not bond.IsAromatic()) and (bond.GetOrder() == 1) and (bond.GetBgn().GetDegree() > 1) and (bond.GetEnd().GetDegree() > 1):
            i = atom_indices[bond.GetBgn().GetName()]
            j = atom_indices[bond.GetEnd().GetName()]
            if (i < j):
               rotatable_bonds.append( (i,j) )
            else:
               rotatable_bonds.append( (j,i) )
               
      print "%d rotatable bond(s) found." % len(rotatable_bonds)
      
      # Search for [ dihedrals ] section.
      indices = extract_section(lines, 'dihedrals') # extract non-blank, non-comment lines
      for index in indices:
         # extract the line
         line = stripcomments(lines[index])         
         # parse the line
         elements = line.split()
         nelements = len(elements)
         # skip unrecognized lines
         if (nelements < 11): continue         
         # extract dihedral atom indices
         i = int(elements[0])
         j = int(elements[1])
         k = int(elements[2])
         l = int(elements[3])
         # if not in our list of rotatable bonds, skip it
         if not (j,k) in rotatable_bonds: continue         
         # function number
         function = int(elements[4])
         if function != 3: raise "Only [dihedrals] function = 3 is supported."
         # C0 - C5 parameters
         C = list()
         for element in elements[5:11]:
            C.append(float(element))

         # append perturbed parameters to end of line.
         line = "    %-4s %-4s %-4s %-4s %3d%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f" % (i, j, k, l, function, C[0], C[1], C[2], C[3], C[4], C[5])
         line += " %12.5f%12.5f%12.5f%12.5f%12.5f%12.5f ; perturbed" % (0.0, 0.0, 0.0, 0.0, 0.0, 0.0) + "\n"

         # replace the line
         lines[index] = line


   # Replace topology file.
   outfile = open(topology_filename, 'w')
   for line in lines:
      outfile.write(line)
   outfile.close()

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
def expandConformations(molecule, maxconfs = None, threshold = None, include_original = False):   
   """Enumerate conformations of the molecule with OpenEye's Omega.

   ARGUMENTS
   molecule (OEMol) - molecule to enumerate conformations for

   OPTIONAL ARGUMENTS
     include_original (boolean) - if True, original conformation is included (default: False)
     maxconfs (integer) - if set to an integer, limits the maximum number of conformations to generated -- maximum of 120 (default: None)
     threshold (real) - threshold in RMSD (in Angstroms) for retaining conformers -- lower thresholds retain more conformers (default: None)
     torsionlib (string) - if a path to an Omega torsion library is given, this will be used instead (default: None)
     
   RETURN VALUES
     expanded_molecule - molecule with expanded conformations
     
   """
   # Initialize omega
   omega = OEOmega()

   # Set maximum number of conformers.
   if maxconfs: omega.SetMaxConfs(maxconfs)
     
   # Set whether given conformer is to be included.
   omega.SetIncludeInput(include_original)
   
   # Set RMSD threshold for retaining conformations.
   if threshold: omega.SetRMSThreshold(threshold) 
 
   # If desired, do a torsion drive.
   if torsionlib: omega.SetTorsionLibrary(torsionlib)

   # Create copy of molecule.
   expanded_molecule = OEMol(molecule)   

   # Enumerate conformations.
   omega(expanded_molecule)

   # return conformationally-expanded molecule
   return expanded_molecule
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

   # Modify gromacs topology file for alchemical free energy calculation.
   perturbGromacsTopology('gromacs.top', molecule)

