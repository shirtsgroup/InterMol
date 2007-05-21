#=============================================================================================# modelPDB.py## Uses Modeller to build PDB starting structures for molecular simulations.## Written 2007-05-20 by# Gregory R. Bowman <gbowman@stanford.edu># Biophysics Program, Stanford University# Pande group#=============================================================================================# TODO:#=============================================================================================# CVS VERSION CONTROL INFORMATION#=============================================================================================# GLOBAL IMPORTS
from modeller import *
from modeller.automodel import *
import os
import os.path
import shutil
import tempfile
#=============================================================================================

class ModelPDB(object):
    """
    A class for making starting PDB structures for simulations.

    Example usage:
import modelPDB
import sys

files = [["1YRI.pdb", "1YRI_wtv.seq"]]

for [pdb_file, seq_file] in files:
    print "PDB File: " + pdb_file
    print "Sequence File: " + seq_file

    ext_ind = seq_file.rindex(".")
    pdbout = seq_file[0:ext_ind] + ".pdb"

    myModel = modelPDB.ModelPDB()
    myModel.makeModel(pdb_file, seq_file, pdbout)
    """

    def getOneLetterCode(self, threeLetterCode):
        """Return the one-letter amino acid symbol, given a three-letter code.

        REQUIRED ARGUMENTS:
          threeLetterCode - The three-letter code.
      
        RETURN VALUES:
          oneLetterCode - The one-letter code.
      
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

        return _res3to1[threeLetterCode.upper()]

    def getCompleteSequence(self, pdbFilename, chain=' '):
        """Extract the complete sequence from a PDB file.

        sequence = getCompleteSequence(pdbfilename, chain=' ')
    
        REQUIRED ARGUMENTS
          pdbFilename - the filename of the PDB file to import from

        OPTIONAL ARGUMENTS
          chain - the one-character chain ID of the chain to import (default ' ')

        RETURN VALUES
          sequence (string) - the complete sequence in one-letter code

        """
    
        # Read the PDB file into memory.
        pdbFile = open(pdbFilename, 'r')
        lines = pdbFile.readlines()
        pdbFile.close()
       
        # Get the complete sequence from the SEQRES fields.
        completeSequence = ""  # one-letter codes
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
                        completeSequence += self.getOneLetterCode(resName)

        # Return first residue and complete sequence.
        return completeSequence


    def makeModel(self, pdbFilename, seqFn, outFile, projName=' ', chain=' '):
        """Creates a PDB file for the sequence by fitting it to the PDB file using Modeller..

        myModel.make_model(pdbFile, seqFile, outFile)
    
        REQUIRED ARGUMENTS
          pdbFilename - the filename of the PDB file to fit to
          seqFn - the filename of the sequence to model, should contain one sequence using one letter code without spaces
          outFile - filename to write final structure to

        OPTIONAL ARGUMENTS
          projName - if not specified a temporary working directory is created and deleted once the model is completed.  If specified, all temporary files will be stored in the directory name specified.
          chain - the one-character chain ID of the chain to import (default ' ')

        RETURN VALUES
         None

        """
    
        # Ensure files exists.
        if not os.path.exists(pdbFilename):
            raise ParameterException, "Specified PDB file %s not found." % pdbFilename
        if not os.path.exists(seqFn):
            raise ParameterException, "Specified sequence file %s not found." % seqFn
    
        # Append full path to pdbfilename.
        pdbFilename = os.path.abspath(pdbFilename)

        # Create a directory for running MODELLER.
        delProj = False # whether to delete working directory at end
        if projName == ' ':
            projName = tempfile.mkdtemp()
            delProj = True
        else:
            if os.path.exists(projName):
                for filename in os.listdir(projName):
                    os.remove(os.path.join(projName,filename))
                os.rmdir(projName)
            os.mkdir(projName)
        projName = os.path.abspath(projName) # full path to working directory

        # Get the complete sequence.
        completeSequence = self.getCompleteSequence(pdbFilename, chain)
        nResidues = len(completeSequence)
        print "PDB sequence: " + completeSequence

        # read in the sequence
        seqFile = open(seqFn, 'r')
        seq = seqFile.readline().strip()
        seqFile.close()
        print "target sequence: " + seq

        # Change working directory to temporary directory.
        oldDir = os.getcwd()
        os.chdir(projName)

        # create the PIR format sequence file for Modeller
        modSeqFn = 'target.ali'
        modSeqFile = open(modSeqFn, 'w')
        print >> modSeqFile, ">P1;target"
        print >> modSeqFile, "sequence:target:::::::0.00: 0.00"
        print >> modSeqFile, seq + "*"
        modSeqFile.close()

        # Create a new environemnt.
        env = environ()

        # do alignment
        aln = alignment(env)
        mdl = model(env, file=pdbFilename)
        aln.append_model(mdl, align_codes='template', atom_files=pdbFilename)        aln.append(file='target.ali', align_codes='target')        aln.align2d()        aln.write(file='alignment.ali', alignment_format='PIR')        aln.write(file='alignment.pap', alignment_format='PAP')

        # make a single model
        a = automodel(env, alnfile='alignment.ali', knowns='template', sequence='target')
        a.starting_model = 1
        a.ending_model = 1
        a.make()

        # go back to old directory
        os.chdir(oldDir)
   
        # copy to specified file (outfile) in next directory up
        shutil.copyfile(os.path.join(projName, "target.B99990001.pdb"), outFile)

        # remove working directory if just using temp dir
        if delProj:
          for filename in os.listdir(projName):
            os.remove(os.path.join(projName, filename))
          os.rmdir(projName)
      
      
