import tempfile
import commands
from openeye.oechem import *
from openeye.oeomega import *
from openeye.oeiupac import *
from openeye.oeshape import *
from openeye.oeproton import *
import os

"""Ligtools:
- get_ligand: Takes a pdb file, extracts specified ligand and tries to protonate and assign bond types. Optionally writes out a output file (of type specified by the filename, i.e. mol2 or pdb) of it; also returns it as an OE mol.
- add_ligand_to_gro: Add a ligand at the end of an existing gro file.
- ligmol2_to_gromacs: Convert a ligand mol2 file to gromacs top and gro files using antechamber and amb2gmx.pl. 
- generate_conf_from_file: Generates (and optionally writes to file) a ligand conformation (or more than one) for a molecule file of arbitrary (OE readable) type; returns it.
- fit_mol_to_refmol: Fit a OE molecule (multi-conformer) to a reference molecule (single conformer, i.e. a ligand structure from a pdb file); write out an output file of the best N matches, where N is specified.
- EnumerateProtonation to enumerate possible protonation states.

Requirements: 
- Amber and Antechamber installations in your path
- MMTOOLSPATH environment variable set to the location of mmtools (for amb2gmx.pl).

By D. Mobley, 6/28/2007."""

def get_ligand(pdbfile, resnum = None, resname = None, outfile = None, chain = None):
   """Open specified pdb file (from provided path), and output the HETATM entries and associated connect entries for the specified ligand. Ligand can be specified by residue number, residue name, or both. Return an OE molecule containing the ligand. If outfile is provided, writes the molecule to the output file also. Also optionally provide a chain ID which, if specified, will search for the specified residue number/name together with the chain ID."""

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
   print workdir
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
   print tleapout
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
       #This will be the end of file, so make sure we don't let a last line straggle in here
       break;
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
       linenum+=1
     else:
       outtext.append(line)
       linenum +=1    

     file = open(outputitp,'w')
     file.writelines(outtext)
     file.close()

def set_subst_name(mol2file, name):
   """Edit specified mol2 file to change the subst_name to name; the OE tools fail to modify this."""
   file = open(mol2file, 'r')
   text = file.readlines()
   file.close()
   
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

   for line in atomsec:
     #If we blindly search and replace, we'll tend to clobber stuff, as the subst_name might be "1" or something lame like that that will occur all over. 
     #If it only occurs once, just replace it
     if line.count(repltext)==1:
       outtext.append( line.replace(repltext, name) )
     else:
       #Otherwise grab the string left and right of the subst_name and sandwich the new subst_name in between. This can probably be done easier in Python 2.5 with partition, but 2.4 is still used someplaces.
       #Loop through the line and tag locations of every non-space entry
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
                  

   for line in text[atomend:]:
     outtext.append(line)
   #Write new mol2
   file = open(mol2file,'w')
   file.writelines(outtext)
   file.close()


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


#Function definition for splitting mol2 files by pose
def split_mol2_poses(infile,outpath,outname):
  """Pass this an input multi-conformer mol2 file, and it will split it into single-pose mol2 files which it will write to the output directory. Relevant stuff:
- infile: Name of input file
- outpath: Path to output directory, which must already exist
- outname: Prefix (without .mol2) of output file names. i.e. "name" will lead to output files "name0.mol2", "name1.mol2", etc., in the output directory."""
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
