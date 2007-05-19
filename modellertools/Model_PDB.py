from modeller import *
from modeller.automodel import *
import os
import os.path
import shutil
import tempfile

class Model_PDB:
  def one_letter_code(self, three_letter_code):
    """Return the one-letter amino acid symbol, given a three-letter code.

    REQUIRED ARGUMENTS:
      three_letter_code - The three-letter code.
      
    RETURN VALUES:
      one_letter_code - The one-letter code.
      
    """
    
    _res3to1 = { 'ALA' : 'A',
                 'ASN' : 'N',
                 'ASP' : 'D',
                 'ARG' : 'R',
                 'CYS' : 'C',
                 'GLN' : 'Q',
                 'GLU' : 'E',               
                 'GLY' : 'G',
                 'HIS' : 'H',                
                 'ILE' : 'I',
                 'LEU' : 'L',
                 'NLE' : 'n',
                 'LYS' : 'K',
                 'MET' : 'M',
                 'PRO' : 'P',
                 'PHE' : 'F',
                 'SER' : 'S',
                 'THR' : 'T',
                 'TRP' : 'W',
                 'TYR' : 'Y',
                 'VAL' : 'V',
                 'UNK' : 'X' }

    return _res3to1[three_letter_code.upper()]

  def getCompleteSequence(self, pdbfilename, chain=' '):
    """Extract the complete sequence from a PDB file.

    sequence = getCompleteSequence(pdbfilename, chain=' ')
    
    REQUIRED ARGUMENTS
      pdbfilename - the filename of the PDB file to import from

    OPTIONAL ARGUMENTS
      chain - the one-character chain ID of the chain to import (default ' ')

    RETURN VALUES
      sequence (string) - the complete sequence in one-letter code

    """

    # Read the PDB file into memory.
    pdbfile = open(pdbfilename, 'r')
    lines = pdbfile.readlines()
    pdbfile.close()
       
    # Get the complete sequence from the SEQRES fields.
    complete_sequence = ""  # one-letter codes
    resid = None
    for line in lines:
        if line[0:6] == 'SEQRES':
            # Parse line into fields.
            field = { }
            field['serNum'] = int(line[8:10])
            field['chainID'] = line[11:12]
            field['numRes'] = line[13:17]
            field['resNames'] = line[19:70].split()

            # Add these residues to the sequence if they belong to the chain of interest.
            if chain == field['chainID']:
                for resName in field['resNames']:                    
                    complete_sequence += self.one_letter_code(resName)

    # Return first residue and complete sequence.
    return complete_sequence


  def make_model(self, pdbfilename, seq_fn, outfile, proj_name=' ', chain=' '):
    """Creates a PDB file for the sequence by fitting it to the PDB file using Modeller..

    myModel.make_model(pdb_file, seq_file, outfile)
    
    REQUIRED ARGUMENTS
      pdbfilename - the filename of the PDB file to fit to
      seq_fn - the filename of the sequence to model, should contain one sequence using one letter code without spaces
      outfile - filename to write final structure to

    OPTIONAL ARGUMENTS
      proj_name - if not specified a temporary working directory is created and deleted once the model is completed.  If specified, all temporary files will be stored in the directory name specified.
      chain - the one-character chain ID of the chain to import (default ' ')

    RETURN VALUES
     None

    """
    
    # Ensure files exists.
    if not os.path.exists(pdbfilename):
        raise ParameterException, "Specified PDB file %s not found." % pdbfilename
    if not os.path.exists(seq_fn):
        raise ParameterException, "Specified sequence file %s not found." % seq_fn

    # Append full path to pdbfilename.
    pdbfilename = os.path.abspath(pdbfilename)

    # Create a directory for running MODELLER.
    del_proj = False
    if proj_name == ' ':
      proj_name = tempfile.mkdtemp()
      del_proj = True
    else:
      if os.path.exists(proj_name):
          for filename in os.listdir(proj_name):
              os.remove(os.path.join(proj_name,filename))
          os.rmdir(proj_name)
      os.mkdir(proj_name)
    proj_name = os.path.abspath(proj_name)

    # Get the complete sequence.
    complete_sequence = self.getCompleteSequence(pdbfilename, chain)
    nresidues = len(complete_sequence)
    print "PDB sequence: " + complete_sequence

    # read in the sequence
    seqfile = open(seq_fn, 'r')
    seq = seqfile.readline().strip()
    seqfile.close()
    print "target sequence: " + seq

    # Change working directory to temporary directory.
    olddir = os.getcwd()
    os.chdir(proj_name)

    # create the PIR format sequence file for Modeller
    mod_seq_fn = 'target.ali'
    mod_seq_file = open(mod_seq_fn, 'w')
    print >> mod_seq_file, ">P1;target"
    print >> mod_seq_file, "sequence:target:::::::0.00: 0.00"
    print >> mod_seq_file, seq + "*"
    mod_seq_file.close()

    # Create a new environemnt.
    env = environ()

    # do alignment
    aln = alignment(env)
    mdl = model(env, file=pdbfilename)
    aln.append_model(mdl, align_codes='template', atom_files=pdbfilename)    aln.append(file='target.ali', align_codes='target')    aln.align2d()    aln.write(file='alignment.ali', alignment_format='PIR')    aln.write(file='alignment.pap', alignment_format='PAP')

    # make a single model
    a = automodel(env, alnfile='alignment.ali', knowns='template', sequence='target')
    a.starting_model = 1
    a.ending_model = 1
    a.make()

    # go back to old directory
    os.chdir(olddir)
   
    # copy to specified file (outfile) in next directory up
    shutil.copyfile(os.path.join(proj_name, "target.B99990001.pdb"), outfile)

    # remove working directory if just using temp dir
    if del_proj:
      for filename in os.listdir(proj_name):
        os.remove(os.path.join(proj_name, filename))
      os.rmdir(proj_name)
      
      