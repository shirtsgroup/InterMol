import tempfile
import commands
from math import *
from openeye.oechem import *
from openeye.oeomega import *
from openeye.oeiupac import *
from openeye.oeshape import *
from openeye.oequacpac import * #DLM added 2/25/09 for OETyperMolFunction; replacing oeproton
from openeye.oeiupac import *
from openeye.oeszybki import *
import os
import re
import shutil

"""Tools for working with small-molecule ligands in setting up free energy calculations with AMBER and gromacs.

REQUIREMENTS

  AmberTools installations (in PATH).  Be sure to download the latest version of Antechamber separately.
  MMTOOLSPATH environment variable set to the location of mmtools (for acpypi.py).

TODO
  Implement scheme described by Christopher Bayly for minimizing artifacts in AM1BCC parameterization by enumerating
  conformations with Omega and minimizing with MMFF with all charges set to their absolute values to minimize intramolecular
  contacts between fragments.
  Move methods that operate on gromacs topology and coordinate files to 'gromacstools'.

AUTHORS

  Originally by David L. Mobley, UCSF, 6/28/2007.
  Rewritten by John D. Chodera, Stanford, 1/20/2008.
  Update by Hideki Fujioka, Tulane, 6/04/2009.
  Update by D. Mobley, University of New Orleans, 6/05/2009
  Additions and fixes by D. Mobley, UNO, 7/2009
  Extended and generalized, D. Mobley, UNO, 8/2009

CHANGELOG
  6/4/2009: HF switched parameterizeForGromacs to use Acpypi from Google code, which is superior to previous amb2gmx.pl
  6/5/2009: DLM minor changes to test case at end; switched perturbGromacsTopology to handle acpypi (rather than amb2gmx) generated topologies -- mainly modifying handling of the comments so these are preserved rather than clobbered.
  6/30/2009: DLM modified parameterizeForGromacs to use shutil rather than commands.getoutput for copying, improving errorchecking; fixed a bug for some residue names that caused files to fail to be copied back.
  7/9/2009: DLM added fitMolToRefmol function to fit a target molecule onto a reference molecule
  7/29/2009: DLM modified assignPartialCharges to work on molecules coming from createMoleculeFromIUPAC; previously they wouldn't as the IUPAC name function doesn't assign atom names and the assignPartialCharges function requires atom names.
  7/31/2009: Removed redundant (older) copy of add_ligand_to_gro
  8/3/2009: DLM: Minor edits to documentation, optional arguments.
  8/4/2009: DLM edited perturbGromacsTopology to add optional argument, vdw_decoupling, that will modify pairs and nonbond_params sections to maintain intramolecular vdw interactions for a molecule which is being deleted. Also made it optional to provide perturbGromacsTopology with a molecule, since this is only used when dihedrals are perturbed (so it is now only required in that case).
  8/11/2009: DLM edited extractMoleculeFromPDB to add option of specifying an altloc typefor hetatm extraction, for example for cases where there are two ligands with residue name AB1 modeled at partial occupancy, distinguished only by altloc flags "A" and "B"
  8/17/2009: DLM edited add_ligand_to_gro to add option to add ligand elsewhere in a gro file, aside from at the very end.
  10/20/2009: DLM fixed a bug in add_ligand_to_topology wherein ligands with two dihedrals sections would not have the contents of one of the sections added to the resulting topology file.
  10/20/2009: DLM incorporating minor changes from Gabe Rocklin into add_ligand_to_gro to fix problems with ligand numbering when combining with protein under some circumstances
  10/27/2009: DLM fixing crash introduced by last fix that occurred when ligands were added at end of gro file (resnum was not defined). 
  10/28/2009: DLM fixing bug in perturbGromacsTopology wherein A & B state charges were nonzero for transformations involving turning off vdw interactions; made some other minor modifications there to make it easier to avoid this problem.
"""

#=============================================================================================
# METHODS FOR READING, EXTRACTING, OR CREATING MOLECULES
#=============================================================================================
def extractMoleculeFromPDB(pdbfile, resnum = None, resname = None, chain = None, outfile = None, altloc = None):
   """Extract a ligand specified in the HETATM records of a PDB file.

   ARGUMENTS
     pdbfile (String) - the name of the PDB file from which the ligand is to be extracted

   OPTIONAL ARGUMENTS
     resnum - limit HETATM extraction to this residue, if specified (default: None)
     resname - limit HETATM extraction to this residue name, if specified (default: None)
     chain - limit HETATM extraction to this chain, if specified (default: None)
     outfile - if specified, the molecule is written to this output file (default: None)
     altloc - if specified, limit HETATM extraction to this altloc identifier (default: None)

   RETURNS
     ligand (OEMol) - the ligand extracted from the PDB file
     
   NOTES
     The molecule will be 'normalized' by protonating it and naming it according to its IUPAC name.
     Parts that are not recognized will be termed 'BLAH'.

   LIMITATIONS/WARNINGS:
     Note that in this approach, if the ligand does not have hydrogen atoms, bond types are assigned based on the bond angles in the PDB structure. Ligands are often poorly modeled in PDB structures, so this will often result in poor bond assignments, for example sp2 carbons may be treated as sp3 carbons, or vise versa. If using this tool, you should ALWAYS check that you are getting the molecule you expect in some other way, such as by inspecting a mol2 file.
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
        taltloc = line[16:17]
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
        if altloc:
            if taltloc != altloc:
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
   os.remove(ofilename)

   # Return the molecule.
   return molecule
#=============================================================================================
def createMoleculeFromIUPAC(name, verbose = False, charge = None):
   """Generate a small molecule from its IUPAC name.

   ARGUMENTS
     IUPAC_name (string) - IUPAC name of molecule to generate

   OPTIONAL ARGUMENTS
     verbose (boolean) - if True, subprocess output is shown (default: False)
     charge (int) - if specified, a form of this molecule with the desired charge state will be produced (default: None)

   RETURNS
     molecule (OEMol) - the molecule

   NOTES
     OpenEye LexiChem's OEParseIUPACName is used to generate the molecle.
     The molecule is normalized by adding hydrogens.
     Omega is used to generate a single conformation.
     Also note that atom names will be blank coming from this molecule. They are assigned when the molecule is written, or one can assign using OETriposAtomNames for example.

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
   #omega.SetVerbose(verbose)
   #DLM 2/27/09: Seems to be obsolete in current OEOmega
   
   omega.SetIncludeInput(False) # don't include input
   omega.SetMaxConfs(1) # set maximum number of conformations to 1
   omega(molecule) # generate conformation      

   if (charge != None):
      # Enumerate protonation states and select desired state.
      protonation_states = enumerateStates(molecule, enumerate = "protonation", verbose = verbose)
      for molecule in protonation_states:
         if formalCharge(molecule) == charge:
            # Return the molecule if we've found one in the desired protonation state.
            return molecule
      if formalCharge(molecule) != charge:
         print "enumerateStates did not enumerate a molecule with desired formal charge."
         print "Options are:"
         for molecule in protonation_states:
            print "%s, formal charge %d" % (molecule.GetTitle(), formalCharge(molecule))
         raise "Could not find desired formal charge."
    
   # Return the molecule.
   return molecule
#=============================================================================================
def readMolecule(filename, normalize = False):
   """Read in a molecule from a file (such as .mol2).

   ARGUMENTS
     filename (string) - the name of the file containing a molecule, in a format that OpenEye autodetects (such as .mol2)

   OPTIONAL ARGUMENTS
     normalize (boolean) - if True, molecule is normalized (renamed, aromaticity, protonated) after reading (default: False)

   RETURNS
     molecule (OEMol) - OEMol representation of molecule

   EXAMPLES
     # read a mol2 file
     molecule = readMolecule('phenol.mol2')
     # works with any type of file that OpenEye autodetects
     molecule = readMolecule('phenol.sdf')
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

   EXAMPLE
     net_charge = formalCharge(molecule)
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

   EXAMPLES
     # read a partial molecule and normalize it
     molecule = readMolecule('molecule.sdf')
     normalizeMolecule(molecule)
   """
   
   # Find ring atoms and bonds
   # OEFindRingAtomsAndBonds(molecule) 
   
   # Assign aromaticity.
   OEAssignAromaticFlags(molecule, OEAroModelOpenEye)   

   # Add hydrogens.
   OEAddExplicitHydrogens(molecule)

   # Set title to IUPAC name.
   name = OECreateIUPACName(molecule)
   molecule.SetTitle(name)

   return molecule
#=============================================================================================
def expandConformations(molecule, maxconfs = None, threshold = None, include_original = False, torsionlib = None, verbose = False):   
   """Enumerate conformations of the molecule with OpenEye's Omega after normalizing molecule. 

   ARGUMENTS
   molecule (OEMol) - molecule to enumerate conformations for

   OPTIONAL ARGUMENTS
     include_original (boolean) - if True, original conformation is included (default: False)
     maxconfs (integer) - if set to an integer, limits the maximum number of conformations to generated -- maximum of 120 (default: None)
     threshold (real) - threshold in RMSD (in Angstroms) for retaining conformers -- lower thresholds retain more conformers (default: None)
     torsionlib (string) - if a path to an Omega torsion library is given, this will be used instead (default: None)
     verbose (boolean) - if True, omega will print extra information

   RETURN VALUES
     expanded_molecule - molecule with expanded conformations

   EXAMPLES
     # create a new molecule with Omega-expanded conformations
     expanded_molecule = expandConformations(molecule)

     
   """
   # Initialize omega
   omega = OEOmega()

   # Set verbosity.
   #omega.SetVerbose(verbose)
   #DLM 2/27/09: Seems to be obsolete in current OEOmega

   # Set maximum number of conformers.
   if maxconfs:
      omega.SetMaxConfs(maxconfs)
     
   # Set whether given conformer is to be included.
   omega.SetIncludeInput(include_original)
   
   # Set RMSD threshold for retaining conformations.
   if threshold:
      omega.SetRMSThreshold(threshold) 
 
   # If desired, do a torsion drive.
   if torsionlib:
      omega.SetTorsionLibrary(torsionlib)

   # Create copy of molecule.
   expanded_molecule = OEMol(molecule)   

   # Enumerate conformations.
   omega(expanded_molecule)


   # verbose output
   if verbose: print "%d conformation(s) produced." % expanded_molecule.NumConfs()

   # return conformationally-expanded molecule
   return expanded_molecule
#=============================================================================================
def assignPartialCharges(molecule, charge_model = 'am1bcc', multiconformer = False, minimize_contacts = False, verbose = False):
   """Assign partial charges to a molecule using OEChem oeproton.

   ARGUMENTS
     molecule (OEMol) - molecule for which charges are to be assigned

   OPTIONAL ARGUMENTS
     charge_model (string) - partial charge model, one of ['am1bcc'] (default: 'am1bcc')
     multiconformer (boolean) - if True, multiple conformations are enumerated and the resulting charges averaged (default: False)
     minimize_contacts (boolean) - if True, intramolecular contacts are eliminated by minimizing conformation with MMFF with all charges set to absolute values (default: False)
     verbose (boolean) - if True, information about the current calculation is printed

   RETURNS
     charged_molecule (OEMol) - the charged molecule with GAFF atom types

   NOTES
     multiconformer and minimize_contacts can be combined, but this can be slow

   EXAMPLES
     # create a molecule
     molecule = createMoleculeFromIUPAC('phenol')
     # assign am1bcc charges
     assignPartialCharges(molecule, charge_model = 'am1bcc')
   """

   #Check that molecule has atom names; if not we need to assign them
   assignNames = False
   for atom in molecule.GetAtoms():
       if atom.GetName()=='':
          assignNames = True #In this case we are missing an atom name and will need to assign
   if assignNames:
      if verbose: print "Assigning TRIPOS names to atoms"
      OETriposAtomNames(molecule)

   # Check input pameters.
   supported_charge_models  = ['am1bcc']
   if not (charge_model in supported_charge_models):
      raise "Charge model %(charge_model)s not in supported set of %(supported_charge_models)s" % vars()

   # Expand conformations if desired.   
   if multiconformer:
      expanded_molecule = expandConformations(molecule)
   else:
      expanded_molecule = OEMol(molecule)
   nconformers = expanded_molecule.NumConfs()
   if verbose: print 'assignPartialCharges: %(nconformers)d conformations will be used in charge determination.' % vars()
   
   # Set up storage for partial charges.
   partial_charges = dict()
   for atom in molecule.GetAtoms():
      name = atom.GetName()
      partial_charges[name] = 0.0

   # Assign partial charges for each conformation.
   conformer_index = 0
   for conformation in expanded_molecule.GetConfs():
      conformer_index += 1
      if verbose and multiconformer: print "assignPartialCharges: conformer %d / %d" % (conformer_index, expanded_molecule.NumConfs())

      # Assign partial charges to a copy of the molecule.
      if verbose: print "assignPartialCharges: determining partial charges..."
      charged_molecule = OEMol(conformation)   
      if charge_model == 'am1bcc':
         OEAssignPartialCharges(charged_molecule, OECharges_AM1BCC)         
      
      # Minimize with positive charges to splay out fragments, if desired.
      if minimize_contacts:
         if verbose: print "assignPartialCharges: Minimizing conformation with MMFF and absolute value charges..." % vars()         
         # Set partial charges to absolute value.
         for atom in charged_molecule.GetAtoms():
            atom.SetPartialCharge(abs(atom.GetPartialCharge()))
         # Minimize in Cartesian space to splay out substructures.
         szybki = OESzybki() # create an instance of OESzybki
         szybki.SetRunType(OERunType_CartesiansOpt) # set minimization         
         szybki.SetUseCurrentCharges(True) # use charges for minimization
         results = szybki(charged_molecule)
         # DEBUG
         writeMolecule(charged_molecule, 'minimized.mol2')
         for result in results: result.Print(oeout)
         # Recompute charges;
         if verbose: print "assignPartialCharges: redetermining partial charges..."         
         OEAssignPartialCharges(charged_molecule, OECharges_AM1BCC)         
         
      # Accumulate partial charges.
      for atom in charged_molecule.GetAtoms():
         name = atom.GetName()
         partial_charges[name] += atom.GetPartialCharge()
   # Compute and store average partial charges in a copy of the original molecule.
   charged_molecule = OEMol(molecule)
   for atom in charged_molecule.GetAtoms():
      name = atom.GetName()
      atom.SetPartialCharge(partial_charges[name] / nconformers)

   # Return the charged molecule
   return charged_molecule
#=============================================================================================
def assignPartialChargesWithAntechamber(molecule, charge_model = 'bcc', judgetypes = None, cleanup = True, verbose = False, netcharge = None):
   """Assign partial charges to a molecule.

   ARGUMENTS
     molecule (OEMol) - molecule for which charges are to be computed

   OPTIONAL ARGUMENTS
     charge_model (string) - antechamber partial charge model (default: 'bcc')
     judgetypes (integer) - if specified, this is provided as a -j argument to antechamber (default: None)
     cleanup (boolean) - clean up temporary files (default: True)
     verbose (boolean) - if True, verbose output of subprograms is displayed
     netcharge (integer) -- if given, give -nc (netcharge) option to antechamber in calculation of charges

   RETURNS
     charged_molecule (OEMol) - the charged molecule with GAFF atom types

   REQUIREMENTS
     antechamber (on PATH)

   WARNING
     This module is currently broken, as atom names get all jacked up during readMolecule() for these mol2 files.
     DLM 4/2/2009: I believe I have fixed the module by switching antechamber to use sybyl atom types. However please note that these will not work as input to tleap and must use GAFF/AMBER atom types there. 
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
   if netcharge:
       chargestr='-nc %d' % netcharge
   else:
       chargestr=''
   command = 'antechamber -i %(uncharged_molecule_filename)s -fi mol2 -o %(charged_molecule_filename)s -fo mol2 -c %(charge_model)s -nc %(formal_charge)d -at sybyl %(chargestr)' % vars()
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
#=============================================================================================
def enumerateStates(molecules, enumerate = "protonation", consider_aromaticity = True, maxstates = 200, verbose = True):
    """Enumerate protonation or tautomer states for a list of molecules.

    ARGUMENTS
      molecules (OEMol or list of OEMol) - molecules for which states are to be enumerated

    OPTIONAL ARGUMENTS
      enumerate - type of states to expand -- 'protonation' or 'tautomer' (default: 'protonation')
      verbose - if True, will print out debug output

    RETURNS
      states (list of OEMol) - molecules in different protonation or tautomeric states

    TODO
      Modify to use a single molecule or a list of molecules as input.
      Apply some regularization to molecule before enumerating states?
      Pick the most likely state?
      Add more optional arguments to control behavior.
    """

    # If 'molecules' is not a list, promote it to a list.
    if type(molecules) != type(list()):
       molecules = [molecules]

    # Check input arguments.
    if not ((enumerate == "protonation") or (enumerate == "tautomer")):
        raise "'enumerate' argument must be either 'protonation' or 'tautomer' -- instead got '%s'" % enumerate

    # Create an internal output stream to expand states into.
    ostream = oemolostream()
    ostream.openstring()
    ostream.SetFormat(OEFormat_SDF)
    
    # Default parameters.
    only_count_states = False # enumerate states, don't just count them

    # Enumerate states for each molecule in the input list.
    states_enumerated = 0
    for molecule in molecules:
        if (verbose): print "Enumerating states for molecule %s." % molecule.GetTitle()
        
        # Dump enumerated states to output stream (ostream).
        if (enumerate == "protonation"): 
            # Create a functor associated with the output stream.
            functor = OETyperMolFunction(ostream, consider_aromaticity, False, maxstates)
            # Enumerate protonation states.
            if (verbose): print "Enumerating protonation states..."
            states_enumerated += OEEnumerateFormalCharges(molecule, functor, verbose)        
        elif (enumerate == "tautomer"):
            # Create a functor associated with the output stream.
            functor = OETautomerMolFunction(ostream, consider_aromaticity, False, maxstates)
            # Enumerate tautomeric states.
            if (verbose): print "Enumerating tautomer states..."
            states_enumerated += OEEnumerateTautomers(molecule, functor, verbose)    
    print "Enumerated a total of %d states." % states_enumerated

    # Collect molecules from output stream into a list.
    states = list()
    if (states_enumerated > 0):    
        state = OEMol()
        istream = oemolistream()
        istream.openstring(ostream.GetString())
        istream.SetFormat(OEFormat_SDF)
        while OEReadMolecule(istream, state):
           states.append(OEMol(state)) # append a copy

    # Return the list of expanded states as a Python list of OEMol() molecules.
    return states


def fitMolToRefmol( fitmol, refmol, maxconfs = None, verbose = False):

    """Fit a multi-conformer target molecule to a reference molecule using OpenEye Shape tookit, and return an OE molecule with the top conformers of the resulting fit. Tanimoto scores also returned.

    ARGUMENTS
      fitmol (OEMol) -- the (multi-conformer) molecule to be fit.
      refmol (OEMol) -- the molecule to fit to

    OPTIONAL ARGUMENTS
      maxconfs -- Limit on number of conformations to return; default return all
      verbose -- Turn verbosity on/off

    RETURNS
      outmol (OEMol) -- output (fit) molecule resulting from fitmol
      scores

    NOTES
      Passing this a multi-conformer fitmol is recommended for any molecule with rotatable bonds as fitting only includes rotations and translations, so one of the provided conformers must already have right bond rotations."""

    #Set up storage for overlay
    best = OEBestOverlay()
    #Set reference molecule
    best.SetRefMol(refmol)

    if verbose:
        print "Reference title: ", refmol.GetTitle()
        print "Fit title: ", fitmol.GetTitle()
        print "Num confs: ", fitmol.NumConfs()

    resCount = 0
    #Each conformer-conformer pair generates multiple scores since there are multiple possible overlays; we only want the best. Load the best score for each conformer-conformer pair into an iterator and loop over it
    scoreiter = OEBestOverlayScoreIter()
    OESortOverlayScores(scoreiter, best.Overlay(fitmol), OEHighestTanimoto())
    tanimotos = [] #Storage for scores
    for score in scoreiter:
        #Get the particular conformation of this match and transform to overlay onto reference structure

        #tmpmol = OEGraphMol(fitmol.GetConf(OEHasConfIdx(score.fitconfidx)))
        tmpmol = OEMol(fitmol.GetConf(OEHasConfIdx(score.fitconfidx)))
        score.Transform(tmpmol)
        #Store to output molecule
        try: #If it already exists
            outmol.NewConf(tmpmol)
        except: #Otherwise
            outmol = tmpmol

        #Print some info
        if verbose:
            print "FitConfIdx: %-4d" % score.fitconfidx,
            print "RefConfIdx: %-4d" % score.refconfidx,
            print "Tanimoto: %.2f" % score.tanimoto
        #Store score
        tanimotos.append(score.tanimoto)
        resCount+=1

        if resCount == maxconfs: break

    return ( outmol, tanimotos )


#=============================================================================================
# METHODS FOR WRITING OR EXPORTING MOLECULES
#=============================================================================================
def writeMolecule(molecule, filename, substructure_name = 'MOL', preserve_atomtypes = False):
   """Write a molecule to a file in any format OpenEye autodetects from filename (such as .mol2).
   WARNING: The molecule will be standardized before writing by the high-level OEWriteMolecule function.
   OEWriteConstMolecule is used, to avoid changing the molecule you pass in.

   ARGUMENTS
     molecule (OEMol) - the molecule to be written
     filename (string) - the file to write the molecule to (type autodetected from filename)

   OPTIONAL ARGUMENTS
     substructure_name (String) - if a mol2 file is written, this is used for the substructure name (default: 'MOL')
     preserve_atomtypes (bool) - if True, a mol2 file will be written with atom types preserved

   RETURNS
     None

   NOTES
     Multiple conformers are written.

   EXAMPLES
     # create a molecule
     molecule = createMoleculeFromIUPAC('phenol')
     # write it as a mol2 file
     writeMolecule(molecule, 'phenol.mol2')
   """

   # Open output stream.
   ostream = oemolostream()
   ostream.open(filename)

   # Define internal function for writing multiple conformers to an output stream.
   def write_all_conformers(ostream, molecule):
      # write all conformers of each molecule
      for conformer in molecule.GetConfs():
         if preserve_atomtypes: OEWriteMol2File(ostream, conformer)
         else: OEWriteConstMolecule(ostream, conformer)
      return

   # If 'molecule' is actually a list of molecules, write them all.
   if type(molecule) == type(list()):
      for individual_molecule in molecule:
         write_all_conformers(ostream, individual_molecule)
   else:
      write_all_conformers(ostream, molecule)

   # Close the stream.
   ostream.close()

   # Replace substructure name if mol2 file.
   suffix = os.path.splitext(filename)[-1]
   if (suffix == '.mol2' and substructure_name != None):
      modifySubstructureName(filename, substructure_name)

   return
#=============================================================================================
def parameterizeForAmber(molecule, topology_filename, coordinate_filename, charge_model = False, judgetypes = None, cleanup = True, show_warnings = True, verbose = False, resname = None, netcharge = None):
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
     netcharge (integer) -- if set, pass this net charge to calculation in antechamber (with -nc (netcharge)), otherwise assumes zero.

   REQUIREMENTS
     acpypi.py conversion script (must be in MMTOOLSPATH)
     AmberTools installation (in PATH)

   EXAMPLES
     # create a molecule
     molecule = createMoleculeFromIUPAC('phenol')
     # parameterize it for AMBER, using antechamber to assign AM1-BCC charges
     parameterizeForAmber(molecule, topology_filename = 'phenol.prmtop', coordinate_filename = 'phenol.crd', charge_model = 'bcc')
   
   """

   # Create temporary working directory and copy ligand mol2 file there.
   working_directory = tempfile.mkdtemp()
   old_directory = os.getcwd()
   os.chdir(working_directory)
   if verbose: print "Working directory is %(working_directory)s" % vars()

   # Write molecule to mol2 file.
   tripos_mol2_filename = os.path.join(working_directory, 'tripos.mol2')
   writeMolecule(molecule, tripos_mol2_filename)

   if resname:
      # Set substructure name (which will become residue name) if desired.
      modifySubstructureName(tripos_mol2_filename, resname)
                 
   # Run antechamber to assign GAFF atom types.
   gaff_mol2_filename = os.path.join(working_directory, 'gaff.mol2')   
   if netcharge:
      chargstr = '-nc %d' % netcharge
   else: chargestr=''
   command = 'antechamber -i %(tripos_mol2_filename)s -fi mol2 -o %(gaff_mol2_filename)s -fo mol2 %(chargestr)s' % vars()
   if judgetypes: command += ' -j %(judgetypes)d' % vars()
   if charge_model:
      formal_charge = formalCharge(molecule)
      command += ' -c %(charge_model)s -nc %(formal_charge)d' % vars()   
   if verbose: print command
   output = commands.getoutput(command)
   if verbose or (output.find('Warning')>=0): print output   
   
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
def parameterizeForGromacs(molecule, topology_filename, coordinate_filename, charge_model = False, cleanup = True, show_warnings = True, verbose = False, resname = None, netcharge=None):
   """Parameterize small molecule with GAFF and write gromacs coordinate/topology files.

   ARGUMENTS
     molecule (OEMol) - OEMol molecule of ligand
     topology_filename (string) - name of output topology file
     coordinate_filename (string) - name of output coordinate file

   OPTIONAL ARGUMENTS
     charge_model (string) - if not False, antechamber is used to assign charges (default: False) -- if set to 'bcc', for example, AM1-BCC charges will be used
     cleanup (boolean) - clean up temporary directories (default: True)
     show_warnings (boolean) - show any warnings during conversion process (default: True)
     verbose (boolean) - show complete output of tools (default: False)
     resname (string) - if set, residue name to use for parameterized molecule (default: None)
     netcharge (integer) - if set, pass this net charge to calculation in antechamber (with -nc (netcharge)), otherwise assume zero.

   REQUIREMENTS
     antechamber (must be in PATH)
     acpypi.py conversion script (must be in MMTOOLSPATH)

   EXAMPLES
     # create a molecule
     molecule = createMoleculeFromIUPAC('phenol')
     # parameterize it for gromacs, using antechamber to assign AM1-BCC charges
     parameterizeForGromacs(molecule, topology_filename = 'phenol.top', coordinate_filename = 'phenol.gro', charge_model = 'bcc')

   """

   # Create temporary directory.
   working_directory = tempfile.mkdtemp()
   old_directory = os.getcwd()
   os.chdir(working_directory)
   if verbose: print "Working directory is %(working_directory)s"

   # Create AMBER coordinate/topology files.
   amber_topology_filename = os.path.join(working_directory, 'amber.prmtop')
   amber_coordinate_filename = os.path.join(working_directory, 'amber.crd')
   parameterizeForAmber(molecule, amber_topology_filename, amber_coordinate_filename, charge_model=charge_model, cleanup=cleanup, show_warnings=show_warnings, verbose=verbose, resname=resname, netcharge = netcharge)
   
   # Use acpypi to convert from AMBER to gromacs topology/coordinates.
   acpypi =  os.path.join(os.getenv('MMTOOLSPATH'), 'converters', 'acpypi.py')
   command = '%(acpypi)s -p %(amber_topology_filename)s -x %(amber_coordinate_filename)s -b OUT' % vars()
   if verbose: print command
   acpypi_output = commands.getoutput(command)
   if verbose: print acpypi_output

   # Restore old directory.
   os.chdir(old_directory)   

   # Copy gromacs topology/coordinates to desired output files.
   shutil.copy( os.path.join(working_directory, 'OUT_GMX.gro'), coordinate_filename )
   shutil.copy( os.path.join(working_directory, 'OUT_GMX.top'), topology_filename)
  
   # Clean up temporary files.
   if cleanup:
      commands.getoutput('rm -r %s' % working_directory)
   else:
      print "Work done in %s..." % working_directory

   return
#=============================================================================================
# METHODS FOR MANIPULATING GROMACS TOPOLOGY AND COORDINATE FILES
#=============================================================================================
def stripcomments(line):
   """Return (line, comments) with whitespace and comments stripped.
   """
   # strip comments
   index = line.find(';')
   comments =''
   if index > -1:
      comments = line[index:]
      line = line[0:index]
   # strip whitespace
   line = line.strip()
   comments = comments.strip()
   # return stripped line
   return line,comments         

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
      line,comments = stripcomments(lines[start_index])
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
      line,comments = stripcomments(lines[end_index])
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

def perturbGromacsTopology(topology_filename, molecule = None, perturb_torsions = True, perturb_vdw = True, perturb_charges = True, perturb_atom_indices = None, vdw_decoupling = False, decouple_atom_types = None):
   """Modify a gromacs topology file to add perturbed-state parameters.

   ARGUMENTS
     topology_file (string) - the name of the topology file to modify

   OPTIONAL ARGUMENTS
     perturb_torsions (boolean) - if True, torsions whose central bond is not in an aromatic ring will be turned off in B state (default: True); note that molecule must also be specified in this case.
     perturb_vdw (boolean) - if True, van der Waals interactions will be turned off in B state (default: True)
     perturb_charges (boolean) - if True, charges will be turned off in B state (default: True)
     perturb_atom_indices (list of ints) - if not None, only atoms with gromacs .top file indices in included range will be perturbed (default: None)
     molecule (OEMol) - molecule corresponding to contents of topology file -- must be the same one used to generate the topology file with parameterizeForGromacs(). NOTE: Required when using perturb_torsions.
     vdw_decoupling (boolean) -- if True (default False), juggle gromacs pairs list and explicit specification of interactions in order to maintain A state intramolecular interactions within the molecule being perturbed. Assumes combination rule 2 and fudgeLJ=0.5.
     decouple_atom_types (list of types): Used only with vdw_decoupling, to specify list of atom types to retain interactions between. If these are not provided, only interactions between perturbed atoms will be retained.

   NOTES
     This code currently only handles the special format gromacs topology files produced by acpypi -- there are allowed variations in format that are not treated here.
     Note that this code also only handles the first section of each kind found in a gromacs .top file.
     Note also that the vdw decoupling portion of this code will work properly only if the atom types being decoupled do not occur elsewhere in the system. 
     The vdw decoupling portion of this code also retains interactions only between *perturbed* atoms unless a separate list of atoms is specified.

   TODO
     Perhaps this method should be combined with 'parameterizeForGromacs' as an optional second step, ensuring that the correct 'molecule' is used.
     Generalize this code to allow it to operate more generally on a specified moleculetype.

   EXAMPLES
     # create a molecule
     molecule = createMoleculeFromIUPAC('phenol')
     # parameterize it for gromacs, using antechamber to assign AM1-BCC charges
     parameterizeForGromacs(molecule, topology_filename = 'phenol.top', coordinate_filename = 'phenol.gro', charge_model = 'bcc')
     # modify topology to prepare it for free energy calculations
     perturbGromacsTopology('phenol.top', molecule = molecule)



   """

   # Read the contents of the topology file.
   infile = open(topology_filename, 'r')
   lines = infile.readlines()
   infile.close()

   # Parse atomtypes to build a list of all atom types for later use.
   atomtypes = list() # storage for atom types
   sigma_by_type = dict() #Store sigmas and epsilons by atom type for later use
   epsilon_by_type = dict()
   indices = extract_section(lines, 'atomtypes')
   for index in indices:
      # extract the line
      line,comments = stripcomments(lines[index])
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
      epsilon_by_type[ atomtype['name'] ] = atomtype['epsilon']
      sigma_by_type[ atomtype['name'] ] = atomtype['sigma']
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
         line = "%(name)-10s%(bond_type)6s      0.0000  0.0000  A %(sigma)13.5e%(epsilon)13.5e ; perturbed\n" % perturbed_atomtype
         # insert the new line
         lines.insert(indices[-1], line)
         indices.append(indices[-1]+1)


   # Process [ atoms ] section
   atoms = list()
   atom_indices = dict()
   type_by_index = dict() #Store atom types by index for later use
   indices = extract_section(lines, 'atoms')
   perturbed_types = [] #Keep track of which atom types are perturbed
   for index in indices:
      # extract the line
      line,comments = stripcomments(lines[index])
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
      type_by_index[ atom['nr'] ] = atom['type']

      # skip if an atom range is specified, and this is not among the atoms we're looking for
      if perturb_atom_indices:
         if atom['nr'] not in perturb_atom_indices: continue
            
         
      # save types we are perturbing for later
      if perturbed_types.count( atom['type'])==0:
          perturbed_types.append( atom['type'] )
      # set perturbation type
      atom['typeB'] = atom['type']
      atom['chargeB'] = atom['charge']
      if perturb_vdw: 
          atom['typeB'] += '_pert' # perturbed vdw type
          #DLM edit 10/28/09 -- set A&B state charges to zero whenever we are turning off vdw
          atom['chargeB'] = 0.0
          atom['charge'] = 0.0
      if perturb_charges: 
          atom['chargeB'] = 0.0 # perturbed charges
      
      # construct a new line
      line = "%(nr)6d %(type)10s %(resnr)6d %(residue)6s %(atom)6s %(cgnr)6d %(charge)10.5f %(mass)10.6f %(typeB)10s %(chargeB)10.5f %(mass)10.6f" % atom + " %(comments)s perturbed\n" % vars()
      
      # replace the line
      lines[index] = line

      # store atoms
      atoms.append(atom)

      # store lookup of atom atom names -> atom numbers
      atom_indices[atom['atom']] = atom['nr']



   #If we are doing vdw decoupling, prepare a nonbond_params section explicitly specifying intramolecular interactions in order to retain them
   if vdw_decoupling:
        if not decouple_atom_types:
            decouple_atom_types = perturbed_types #Types to retain interactions between

        #Build section
        nonbond_sec = [ '[ nonbond_params ]\n' ]
        #Note that below could be rewritten more simply now that we are storing sigma and epsilon by each type
        for type1 in decouple_atom_types:
            idx1 = decouple_atom_types.index( type1)
            for type2 in decouple_atom_types:
                idx2 = decouple_atom_types.index(type2)
                if idx2 >= idx1: #Prevent adding same lines twice
                    sigma = 0.5*( sigma_by_type[type1] + sigma_by_type[type2] )
                    epsilon = sqrt( epsilon_by_type[type1]*epsilon_by_type[type2] )
                    #Add to parameters section
                    nonbond_sec.append( type1 + '   ' + type2 + '   1   ' + str( sigma) + '  ' + str(epsilon) +'\n' )
                    #Also add for perturbed atoms
                    nonbond_sec.append( type1+'_pert' + '   ' + type2+'_pert' + '   1   ' + str( sigma) + '  ' + str(epsilon) +'\n' )
                
        #Add nonbond_sec into lines in appropriate place
        molstart = lines.index('[ moleculetype ]\n')
        for line in nonbond_sec:
            lines.insert( molstart, line)
            molstart+=1
        lines.insert(molstart, '\n')


   #If we are doing vdw decoupling, modify the pairs section to explicitly state interactions rather than having them determined from the parameters
   if vdw_decoupling:
        indices = extract_section(  lines, 'pairs' )
        for index in indices:
            # Extract line
            line, comments = stripcomments( lines[index] )
            # Parse line
            elements = line.split()
            nelements = len(elements)
            # skip if not all elements found
            if (nelements < 3): continue
            #parse line
            pair = dict()
            pair['i'] = int(elements[0])
            pair['j'] = int(elements[1])
            pair['type'] = int(elements[2])
            #Check for potential problems
            if not pair['type']==1:
                raise 'Error: ligandtools.py only knows how to do vdw decoupling for pairs type 1'
            #Look up atomtypes
            type_i = type_by_index[ pair['i']]
            type_j = type_by_index[ pair['j']]
            #Look up sigma and epsilon
            epsilon_i = epsilon_by_type[ type_i ]
            epsilon_j = epsilon_by_type[ type_j ]
            sigma_i = sigma_by_type[ type_i ]
            sigma_j = sigma_by_type[ type_j ]
            #Compute total
            sigma = 0.5*(sigma_i + sigma_j)
            fudgeLJ = 0.5
            epsilon = fudgeLJ*sqrt(epsilon_i*epsilon_j)
            #Build new line
            line = line.replace('\n','')
            line = line + '   ' + str(sigma) + '   ' + str(epsilon) + '\n'
            # replace the line
            lines[index] = line

   # Process [ bonds ] section
   indices = extract_section(lines, 'bonds')
   for index in indices:
      # extract the line
      line,comments = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 5): continue
      # parse line
      bond = dict()
      bond['i'] = int(elements[0])      
      bond['j'] = int(elements[1])
      bond['function'] = int(elements[2])
      bond['Req'] = float(elements[3])
      bond['Keq'] = float(elements[4])
      # skip if an atom range is specified, and this is not among the atoms we're looking for
      if perturb_atom_indices:
         if (bond['i'] not in perturb_atom_indices) and (bond['j'] not in perturb_atom_indices): continue      
      # construct a new line
      line = "%(i)5d %(j)5d %(function)5d%(Req)12.4e%(Keq)12.4e" % bond
      line += " %(Req)12.4e%(Keq)12.4e" % bond + comments+'\n'
      # replace the line
      lines[index] = line

   # Process [ angles ] section
   indices = extract_section(lines, 'angles')
   for index in indices:
      # extract the line
      line,comments = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 6): continue
      # parse line
      angle = dict()
      angle['i'] = int(elements[0])      
      angle['j'] = int(elements[1])
      angle['k'] = int(elements[2])
      angle['function'] = int(elements[3])
      angle['theta'] = float(elements[4])
      angle['cth'] = float(elements[5])
      # skip if an atom range is specified, and this is not among the atoms we're looking for
      if perturb_atom_indices:
         if (angle['i'] not in perturb_atom_indices) and (angle['j'] not in perturb_atom_indices) and (angle['k'] not in perturb_atom_indices): continue      
      # construct a new line
      line = "%(i)5d %(j)5d %(k)5d %(function)5d%(theta)12.4e%(cth)12.4e" % angle
      line += " %(theta)12.4e%(cth)12.4e" % angle + comments+'\n'
      # replace the line
      lines[index] = line

   # Set rotatable bond torsions in B state to zero, if desired.
   if perturb_torsions:
      # Determine list of rotatable bonds to perturb.
      rotatable_bonds = list()
      if not molecule:
          raise 'Error: OE molecule must be provided corresponding to the topology must be providded when using pertub_torsions' #Quick error check

      for bond in molecule.GetBonds():
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
         line,comments = stripcomments(lines[index])         
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

         # skip if an atom range is specified, and this is not among the atoms we're looking for
         if perturb_atom_indices:
            if (i not in perturb_atom_indices) and (j not in perturb_atom_indices) and (k not in perturb_atom_indices) and (l not in perturb_atom_indices): continue      

         # function number
         function = int(elements[4])
         if function != 3: raise "Only [dihedrals] function = 3 is supported."
         # C0 - C5 parameters
         C = list()
         for element in elements[5:11]:
            C.append(float(element))

         # reconstruct perturbed line
         line = "    %-4s %-4s %-4s %-4s %3d%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f" % (i, j, k, l, function, C[0], C[1], C[2], C[3], C[4], C[5])         
         if (j,k) in rotatable_bonds:
            # perturb rotatable bonds
            line += " %12.5f%12.5f%12.5f%12.5f%12.5f%12.5f %s perturbed" % (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, comments) + "\n"
         else:
            # don't perturb 
            line += " %12.5f%12.5f%12.5f%12.5f%12.5f%12.5f" % (C[0], C[1], C[2], C[3], C[4], C[5]) + comments + "\n"

         # replace the line
         lines[index] = line
         
   # Replace topology file.
   outfile = open(topology_filename, 'w')
   for line in lines:
      outfile.write(line)
   outfile.close()

   return
#=============================================================================================
def totalCharge(topology_filename):
   """Determine the total charge of the first molecule molecule in a .top file.
   
   ARGUMENTS
     topology_filename (string) - the name of the topology file to parse

   NOTE
     Only the first [ atoms ] section is read, and only charges from A state are read.
   """
    
   # Read the contents of the topology file.
   infile = open(topology_filename, 'r')
   lines = infile.readlines()
   infile.close()

   total_charge = 0.0

   # Process [ atoms ] section
   indices = extract_section(lines, 'atoms')
   for index in indices:
      # extract the line
      line,comments = stripcomments(lines[index])
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

      total_charge += atom['charge']

   # return the total charge
   return total_charge
     
#=============================================================================================
def add_ligand_to_gro(targetgro, liggro, outgro, resname = 'TMP', add_after_resnum = None):
    """Append ligand coordinates to the end of an existing gromacs .gro file.

   ARGUMENTS
     targetgro (string) - gromacs .gro file to which ligand coordinates are to be appended (not modified)
     liggro (string) - gromacs .gro file containing ligand coordinates (not modified)
     outgro (string) - filename of gromacs .gro file to contain merged target and ligand coordinates (created)

   OPTIONAL ARGUMENTS
     resname (string) - name for ligand residue (default: 'TMP')
     add_after_resnum (integer) -- Default None. If this is specified, ligand will be added to gro file after a specified residue number (i.e., after the protein, for example) rather than at the end of the file. 

   NOTES
     No effort is made to adjust box size.
     Ligand must be a single residue.
     New residue number is derived from last residue in target .gro file.
     In the case of add_after_resnum, everything is done as normally in terms of residue number and atom number calculation, EXCEPT that the ligand is inserted into the gro file after the specified residue number and then trjconv is used to correct the residue and atom numbering to be consecutive.
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

    if not add_after_resnum:
        # Create new gromacs .gro file in memory.
        outtext = [ targetlines[0] ]
        outtext.append(' %s\n' % newatomnum)
        for line in targetlines[2:-1]:
            outtext.append(line)
        #Compute residue number
        thisres = line.split()[0]
        i = len(thisres)
        while not thisres[0:i].isdigit():
            i-=1
        resnum = int( thisres[0:i])
        atomnum = int(line[15:20])
    else: #If we want to add the ligand after a specified residue number rather than at the end
        #Find the line corresponding to residue number we want to add
        residueline = 1
        resnum = 1
        #Scan until we get to where we want to be
        while resnum <= add_after_resnum:
            residueline+=1
            #Get number of residue
            thisres = targetlines[ residueline ].split()[0]
            i = len(thisres)
            while not thisres[0:i].isdigit():
                i-=1
            resnum = int(thisres[0:i])
            atomnum = int(targetlines[residueline][15:20]) #DLM add 10/21/09
        #Create new gromacs .gro file in memory
        outtext = [ targetlines[0] ]
        outtext.append(' %s\n' % newatomnum)
        for line in targetlines[2:residueline]:
            outtext.append(line)

    # Append the ligand coordinate lines, renumbering atom and residue numbers.
    resnumname='%4s%-4s' % (resnum+1, resname)
    for line in liglines[2:-1]:
        anum = int( line[15:20].split()[0] )
        newanum = atomnum+anum
        line = ' '+resnumname+line[9:15]+('%5s' % newanum)+line[20:]
        outtext.append(line)

    #If we are not adding at the end, add the rest of the file
    if add_after_resnum:
        for line in targetlines[residueline:-1]:
            outtext.append(line)

    # Add box line from target .gro file to end of file.
    outtext.append(targetlines[-1])

    # Write modified .gro file.
    file = open(outgro, 'w')
    file.writelines(outtext)
    file.close()

    #If we wrote it into the middle of the gro file, use trjconv to correct residue/atom numbering
    if add_after_resnum!=None:
        commands.getoutput('echo 0 | trjconv -f %(outgro)s -s %(outgro)s -o %(outgro)s' % vars() )

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
     The transformation is only applied to the first molecule in the mol2 file.

   TODO
     This function is still difficult to read.  It should be rewritten to be comprehensible by humans.
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
def merge_protein_ligand_topologies(prottop,protgro,ligtop,liggro,complextop,complexgro):
   """Merge gromacs topology and coordinate files for a ligand into those of a protein to form a complex where protein and ligand are merged into a single molecule.

   ARGUMENTS
     prottop (string) - filename of protein topology file, potentially solvated (gromacs .top file)
     protgro (string) - filename of protein coordinate file (gromacs .gro file)
     ligtop (string) - filename of ligand topology file
     liggro (string) - filename of ligand coordinate file
     complextop (string) - filename of topology file for complex to be written
     complexgro (string) - filename of coordinate file for complex to be written
     
   NOTES
     Ordinarily, we want to use restraints between the ligand and the protein, and GROMACS currently requires that all atoms which will be involved in restraints be in the same [ moleculetype ] section. THIS requirement means that they must be sequential in the .gro file. Therefore, this script will not only add the ligand in the same molecules section of the topology file, but it will also edit the protein .gro file to renumber all of the solvent to add the ligand before it, immediately following the protein.

   """
   
   # Merge topology files.
   add_ligand_to_topology(prottop,ligtop,complextop)
   # Merge coordinate files.  
   add_ligand_to_gro(protgro,liggro,complexgro)

   return
#=============================================================================================
def read_next_topo_section(textarray):
   """Reads the [ named ] section in the GROMACS topology file from the passed textarray and 
   returns an array containing only that section (stopping before the subsequent section). 
   len(section) will be the number of lines in that section. Also returns a second array 
   containing that portion of textarray prior to the beginning of the next section."""
   
   #Pattern to match for start
   secstart=re.compile(r'\s*\[.+\]')

   #I think this is actually robust for comments: Comments shouldn't be recognized
   #by the search string above (nor should any bracketed object not preceded by
   #spaces or nothing and so commented out sections will be treated as 'headers' and left untouched.
   
   #Find matches
   section=[]
   foundstart=False
   startline=0
   ctr=0
   for line in textarray:
     #search for match
     m=secstart.match(line)
     #If there is a match
     if m:
       #If this is the first match
       if not foundstart:
          #We have found the start of the section we want 
          foundstart=True
          #Save this line
          section.append(line)
          startline=ctr
          continue #Go back to the top of the loop (the next line)
       #If this is not the first match
       elif foundstart:
          #Don't keep going into the next section
          break
     #Otherwise if this is not a match, so it doesn't start or end a section. 
     #Only save if we already got to the section.
     elif foundstart:
       section.append(line)
     ctr+=1
   #End loop over array; we've now got the section.
   
   #Return the extra portion also
   extra=textarray[0:startline]
   return section,extra
#=============================================================================================
def add_ligand_to_topology(prottop,ligtop,complextop):
   """Append a ligand gromacs topology file to a protein topology file, merging into a single molecule, producing a complex topology file.

   ARGUMENTS
     prottop (string) - filename of protein topology file to read
     ligtop (string) - filename of ligand topology file to read
     complextop (string) - filename of topology file to create of merged complex     

   """

   ####
   #READ INPUT
   #####
   #Read protein and ligand topology files into arrays.
   infile=open(prottop,'r')
   prottext=infile.readlines()
   infile.close()
   infile=open(ligtop,'r')
   ligtext=infile.readlines()
   infile.close()

   ####
   #SET UP TO PARSE
   ####
   #Section names:
   secname=re.compile(r'\[\s*(?P<name>\w*)\s*\]')
   #To strip lines beginning with comments:
   comment=re.compile('\s*;.*')
   #To parse protein atoms section
   atomline=re.compile(r'(?P<lspc>\s*)(?P<nr>\d+)(?P<nr_spc>\s+)(?P<type>\w+\d*[_]*\d*\s+)'
   r'(?P<resnr>\d+)(?P<middle>\s+\w+\s+\w+\d*[*]*\d*\s+)(?P<cgnr>\d+)(?P<end>\s+.*)')
  
   #To parse bond lines:
   bondline=re.compile(r'(?P<lspc>\s*)(?P<ai>\d+)\s+(?P<aj>\d+)(?P<end>\s+\d+.*)')
   #To parse angle lines:
   angleline=re.compile(r'(?P<lspc>\s*)(?P<ai>\d+)\s+(?P<aj>\d+)\s+(?P<ak>\d+)(?P<end>\s+\d+.*)')
   #To parse dihedral lines:
   dihline=re.compile(r'(?P<lspc>\s*)(?P<ai>\d+)\s+(?P<aj>\d+)\s+(?P<ak>\d+)\s+(?P<al>\d+)'
   r'(?P<end>\s+\d+.*)')
   #To parse pairs lines
   pairline=re.compile(r'\s*(?P<ai>\d+)\s+(?P<aj>\d+)(?P<end>\s+.*)')

   #Store information from ligand
   ligsection={}
   #Sections to copy from ligand topology file verbatim (i.e., not merge with those in protein topology)
   #Below these are assumed to precede the [ moleculetypes ] section when added to the protein
   #topology so if additional sections are add for which this is not true, that will need to be 
   #generalized.
   copylist=['atomtypes','nonbond_params']  
   #Sections to merge with those in protein topology
   #Adding additional sections here will extract information from ligand topology, but additional information is required below on how to merge these with protein topology.
   mergelist=['atoms','bonds','angles','dihedrals','pairs']
   savelist=mergelist+copylist
   #Those not listed in either of these locations will be discarded (i.e. moleculetype) 

   ####
   #READ LIGAND FILE
   ####
   endindex=len(ligtext)
   loc=0
   while loc<endindex:
      #Read next chunk
      (sec,headers)=read_next_topo_section(ligtext[loc:])
      #Find match for name
      m=secname.match(sec[0])
      name=m.group('name')
      #Find location of end of section
      endloc=loc+len(sec)+len(headers)

      #Now decide what to do.
      #Ignore all headers for the ligand.

      #Now extract information from the sections we need.
      #Is this a section we want to save?
      if name in savelist:
        #Make it an empty array for starters, unless we have already created it.
        #This may sometimes be the case, for example if there are multiple [dihedrals] sections
        if not ligsection.has_key(name):
           ligsection[name]=[]
        #Save lines in the section except those beginning with comments
        for line in sec:
          m=comment.match(line)
          if not m:
            ligsection[name].append(line)
        
      #Done saving ligand information
      loc=endloc

   ####
   #READ PROTEIN FILE AND APPEND LIGAND STUFF
   ####
   
   #protein properties
   resnum=0
   atomnum=0
   cgnrnum=0
   #Make sure only copy ligand dihedrals once (there are sometimes two dihedrals sections for the protein)
   copiedDih=False
   #store new output
   newtop=[]
   #read
   endindex=len(prottext)
   loc=0
   while loc < endindex:
      #Read next chunk
      (sec,headers)=read_next_topo_section(prottext[loc:])
      #Find match for name
      m=secname.match(sec[0])
      name=m.group('name')
      #Find location of end of section
      endloc=loc+len(sec)+len(headers)

      #First do headers
      for line in headers:
         newtop.append(line)

      #Now decide what to do with section
      #Any atomtypes and nonbonded params come before the moleculetypes section.
      if name=='moleculetype':
         #For each ligand section to copy
         for ligname in copylist:
           #If I've found this section (i.e. there will be no nonbonded params for charge topology)    
           if ligsection.has_key(ligname):
              #For each line here, append the line
              for line in ligsection[ligname]:
                newtop.append(line)
         #Now add moleculetypes section.
         for line in sec:
            newtop.append(line)
      #Otherwise, go ahead and append the protein section
      else:
         for line in sec:
           newtop.append(line)
         #Also, if this is one of the sections for which we also have ligand information, we need to 
         #combine that information. This usually takes some renumbering.
         if name in mergelist:
             #Pop blank end-of-section line off list
             newtop.pop()
             if name=='atoms':
               #First get how many atoms there are in protein
               lastprotline=sec[len(sec)-2]
               m=atomline.match(lastprotline)
               #For debugging: Print line
               resnum=int(m.group('resnr'))
               atomnum=int(m.group('nr'))
               cgnrnum=int(m.group('cgnr'))

               #Change atom number, residue number, and cgnr in ligand section to correct values
               for line in ligsection[name]:
                 m=atomline.match(line)
                 if m:
                   newresnr=str(int(m.group('resnr'))+resnum)
                   newatnr=str(int(m.group('nr'))+atomnum)
                   newcgnr=str(int(m.group('cgnr'))+cgnrnum)
                   newline=m.group('lspc')+newatnr+m.group('nr_spc')+m.group('type')+newresnr+m.group('middle')+newcgnr+m.group('end')+'\n'
                   #Append this line
                   newtop.append(newline)
          
             elif name=='bonds':
               for line in ligsection[name]:
                 m=bondline.match(line)
                 if m:
                   newatomi=str(int(m.group('ai'))+atomnum)
                   newatomj=str(int(m.group('aj'))+atomnum)
                   newline=m.group('lspc')+newatomi+'   '+newatomj+m.group('end')+'\n'
                   newtop.append(newline)         

             elif name=='angles':
                for line in ligsection[name]:
                  m=angleline.match(line)
                  if m:
                    newatomi=str(int(m.group('ai'))+atomnum)
                    newatomj=str(int(m.group('aj'))+atomnum)
                    newatomk=str(int(m.group('ak'))+atomnum)
                    newline=m.group('lspc')+newatomi+'   '+newatomj+'   '+newatomk+m.group('end')+'\n'
                    newtop.append(newline)

             elif name=='dihedrals':
                if not copiedDih:
                  #Some molecules (i.e. methane) have no dihedrals, in which case these do not need to be added and trying to add them will result in a key error, so check:
                  if ligsection.has_key(name):
                    for line in ligsection[name]:
                      m=dihline.match(line)
                      if m:
                        newatomi=str(int(m.group('ai'))+atomnum)
                        newatomj=str(int(m.group('aj'))+atomnum)
                        newatomk=str(int(m.group('ak'))+atomnum)
                        newatoml=str(int(m.group('al'))+atomnum)
                        newline=m.group('lspc')+newatomi+'   '+newatomj+'   '+newatomk+'   '+newatoml+m.group('end')+'\n'
                        newtop.append(newline)
                  #Prevent writing ligand dihedrals twice
                  copiedDih=True;

             elif name=='pairs':
                for line in ligsection[name]:
                  m=pairline.match(line)
                  if m:
                    newatomi=str(int(m.group('ai'))+atomnum)
                    newatomj=str(int(m.group('aj'))+atomnum)
                    newline=newatomi+'  '+newatomj+m.group('end')+'\n'
                    newtop.append(newline)


             #Done with if statements for sections
             #Write blank end-of-section line since I popped the others off.
             newtop.append('\n')
      #OK. Now we're done copying, etc. Update location
      loc=endloc
     
   #Exit loop over protein topology
   #Write new topology
   file=open(complextop,'w')
   file.writelines(newtop)
   file.close()
   #DONE.
   return
#=============================================================================================
# TEST DRIVER
if __name__ == '__main__':

   # Test all capabilities of ligandtools.   

   # Create a molecule.
   molecule = createMoleculeFromIUPAC('phenol')
   # molecule = createMoleculeFromIUPAC('biphenyl')
   # molecule = createMoleculeFromIUPAC('4-(2-Methoxyphenoxy) benzenesulfonyl chloride')
   # molecule = createMoleculeFromIUPAC('n-pentane')   
   # molecule = createMoleculeFromIUPAC('n-decane')

   # Write mol2 file for the molecule.
   writeMolecule(molecule, 'phenol.mol2')

   # Expand conformations
   expanded_molecule = expandConformations(molecule, verbose = True)
   writeMolecule(expanded_molecule, 'phenol-expanded.mol2')

   # Charge the molecule with OpenEye tools.
   charged_molecule = assignPartialCharges(molecule, charge_model = 'am1bcc', verbose = True)
   writeMolecule(charged_molecule, 'phenol-am1bcc-openeye.mol2')

   # Charge the molecule with OpenEye tools using average charges over multiple conformations.
   multiconformer_charged_molecule = assignPartialCharges(molecule, charge_model = 'am1bcc', multiconformer = True, verbose = True)
   writeMolecule(multiconformer_charged_molecule, 'phenol-am1bcc-openeye-multiconformer.mol2')

   # Charge the molecule with OpenEye tools using trick from C. Bayly to minimize intramolecular contacts
   minimizedcontacts_charged_molecule = assignPartialCharges(molecule, charge_model = 'am1bcc', minimize_contacts = True, verbose = True)
   writeMolecule(minimizedcontacts_charged_molecule, 'phenol-am1bcc-openeye-minimizedcontacts.mol2')   

   # Write GAFF parameters for AMBER
   parameterizeForAmber(charged_molecule, topology_filename = 'phenol.prmtop', coordinate_filename = 'phenol.crd', resname = 'PHE')

   # Write GAFF parameters for gromacs.
   parameterizeForGromacs(charged_molecule, topology_filename = 'phenol.top', coordinate_filename = 'phenol.gro', resname = 'PHE')

   # Write GAFF parameters for gromacs, using antechamber to generate AM1-BCC charges.
   parameterizeForGromacs(molecule, topology_filename = 'phenol-antechamber.top', coordinate_filename = 'phenol-antechamber.gro', charge_model = 'bcc', resname = 'PHE')

   # Modify gromacs topology file for alchemical free energy calculation after storing unperturbed copy
   os.system('cp phenol.top phenol_unperturbed.top')
   perturbGromacsTopology('phenol.top', molecule = molecule)

   # Convert .top file to .itp file
   top_to_itp('phenol.top', 'phenol.itp', moleculetype = 'phenol')

