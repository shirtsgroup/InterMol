#=============================================================================================
# modelPDB.py
#
# Uses Modeller to build PDB starting structures for molecular simulations.
#
# Written 2007-05-20 by
# Gregory R. Bowman <gbowman@stanford.edu>
# Biophysics Program, Stanford University
# Pande group
#=============================================================================================
# TODO:
# - Merge modelMissingAtoms and makeModel.
# - Add ability to model missing residues by reading residues that are present in template PDB file and minding chain breaks.
# - Add option (as default) to generate multiple models (maybe 50?) and score them with DOPE to generate higher-quality results
# - Add option to have target sequence read from SEQRES/DBREF blocks in template PDB if no sequence is specified. Use getCompleteSequence.
# - Allow automatic selection of modeling protocol: Copy coordinates if possible, or not if a mutation is modeled.
# - Move some static methods outside the class.
# - What should initializer do?
# - Check to make sure disulfide bonds are handled correctly.
# - Turn off MODELLER debug output.
#=============================================================================================
# CHANGELOG:
# 2007-08-20 VAV Added ACE and NH2 cappipn
# 2007-06-20 JDC Fixed a bug in processing of SEQADV, but still not sure it's robust.
# 2007-06-19 JDC Added new method modelMissingAtoms to model missing atoms and residues.
# 2007-05-20 GRB Created, debugged.
#=============================================================================================
# CVS VERSION CONTROL INFORMATION
#=============================================================================================
# GLOBAL IMPORTS
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb 

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

    def getCompleteSequenceSimple(self, pdbFilename, chain=' '):
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

    def getCompleteSequence(self, pdbfilename, chain=' '):
        """Extract the complete sequence and starting residue index from a PDB file.

        first_residue_id, sequence = getCompleteSequence(pdbfilename, chain=' ')

        REQUIRED ARGUMENTS
            pdbfilename - the filename of the PDB file to import from

        OPTIONAL ARGUMENTS
            chain - the one-character chain ID of the chain to import (default ' ')

        RETURN VALUES
            first_residue (integer) - the index of the first residue
            sequence (string) - the complete sequence in one-letter code

        """

        # Read the PDB file into memory.
        pdbfile = open(pdbfilename, 'r')
        lines = pdbfile.readlines()
        pdbfile.close()

        # Extract the DBREF field from the header to determine the span of residues provided in the SEQRES fields.
        first_residue = None
        for line in lines:
            if line[0:6] == "DBREF ":
                # Parse line into fields.
                field = { }
                field["idCode"] = line[7:11]
                field["chainID"] = line[12:13]
                field["seqBegin"] = int(line[14:18])
                field["insertBegin"] = line[18:19]
                field["seqEnd"] = int(line[20:24])
                field["insertEnd"] = line[24:25]
                field["database"] = line[26:32]
                field["dbAccession"] = line[33:41]
                field["dbIdCode"] = line[42:54]
                field["dbseqBegin"] = line[55:60]
                field["idbnsBeg"] = line[60:61]

                # Store initial and final residues.
                first_residue = field["seqBegin"]
                last_residue = field["seqEnd"]

        # Raise an exception if the PDB file does not contain a DBREF field.
        if not first_residue:
            raise RuntimeError, "No DBREF field found in PDB file -- noncompliant with PDB 2.1 format."

        # Process SEQADV records, if present.
        # COLUMNS       DATA TYPE       FIELD      DEFINITION
        # -----------------------------------------------------------------
        #  1 -  6       Record name     "SEQADV"
        #  8 - 11       IDcode          idCode    ID code of this entry.
        # 13 - 15       Residue name    resName   Name of the PDB residue in conflict.
        # 17            Character       chainID   PDB chain identifier.
        # 19 - 22       Integer         seqNum    PDB sequence number.
        # 23            AChar           iCode     PDB insertion code.
        # 25 - 28       LString         database  
        # 30 - 38       LString         dbIdCode  Sequence database accession number.
        # 40 - 42       Residue name    dbRes     Sequence database residue name.
        # 44 - 48       Integer         dbSeq     Sequence database sequence number.
        # 50 - 70       LString         conflict  Conflict comment.        
        for line in lines:
            if line[0:6] == 'SEQADV':
                # DEBUG
                print line,
                # Parse line into fields.
                field = { }
                field['idCode'] = line[7:11]
                field['resName'] = line[12:15]
                field['chainID'] = line[16:17]
                field['seqNum'] = int(line[18:22])
                field['iCode'] = line[22:23]
                field['database'] = line[24:28]
                field['dbIdCode'] = line[29:38]
                field['dbRes'] = line[39:42]
                field['dbSeq'] = line[43:48]
                field['conflict'] = line[49:70]

                # If SEQADV has an earlier-numbered residue than DBREF, change the first or last residue number.
                if (chain == field['chainID']) and (field['seqNum'] < first_residue):
                    first_residue = field['seqNum']
                if (chain == field['chainID']) and (field['seqNum'] > last_residue):
                    last_residue = field['seqNum']                    

        # Get the complete sequence from the SEQRES fields and store in a dictionary.
        print "Processing SEQRES fields..."
        sequence = { }
        seqNum = first_residue
        resid = None
        for line in lines:
            if line[0:6] == 'SEQRES':
                # DEBUG
                print line,
                # Parse line into fields.
                field = { }
                field['serNum'] = int(line[8:10])
                field['chainID'] = line[11:12]
                field['numRes'] = line[13:17]
                field['resNames'] = line[19:70].split()

                # Add these residues to the sequence if they belong to the chain of interest.
                if chain == field['chainID']:
                    for resName in field['resNames']:
                        sequence[seqNum] = resName
                        seqNum += 1

        # DEBUG
        print "sequence:"
        ordered = sequence.keys()
        ordered.sort()
        for key in ordered:
            print "%5d %3s" % (key, sequence[key])

        # TODO: Check to ensure sequence is contiguous.

        # Find index of first residue.
        firstSeqNum = min(sequence.keys())
        lastSeqNum = max(sequence.keys())
        
        # Reconstruct sequence in one-letter-code
        complete_sequence = ""
        for seqNum in range(firstSeqNum, lastSeqNum+1):
            resName = sequence[seqNum]
            complete_sequence += self.getOneLetterCode(resName)
        first_residue = firstSeqNum

        # DEBUG
        print "complete sequence: %s" % complete_sequence

        # Return first residue and complete sequence.
        return first_residue, complete_sequence

    def getPresentSequence(self, pdbfilename, chain=' '):
        """Extract the sequence of residues for which any ATOM records define coordinates coordinates in a PDB file.
        
        present_sequence = getPresentSequence(pdbfilename, chain=' ')
    
        REQUIRED ARGUMENTS
          pdbfilename - the filename of the PDB file from which the sequence is to be extracted

        OPTIONAL ARGUMENTS
          chain - the one-character chain ID of the chain to import (default ' ')

        RETURN VALUES
          sequence (dictionary) - sequence[residue_id] is the one-letter code corresponding to residue index residue_id

        The ATOM records are read, and the sequence for which there are atomic coordinates is stored.

        """
        
        # Read the PDB file into memory.
        pdbfile = open(pdbfilename, 'r')
        lines = pdbfile.readlines()
        pdbfile.close()

        # Extract the sequence for which there are defined atomic coordinates.
        sequence = { }
        for line in lines:
            if line[0:6] == "ATOM  ":
                # Parse line into fields.
                field = { }
                field["serial"] = int(line[6:11])
                field["name"] = line[12:16]
                field["altLoc"] = line[16:17]
                field["resName"] = line[17:20]
                field["chainID"] = line[21:22]
                field["resSeq"] = int(line[22:26])
                field["iCode"] = line[26:27]
                
                # Add these residues to the sequence if they below to the chain of interest.
                if (chain == field['chainID']):
                    sequence[field["resSeq"]] = self.getOneLetterCode(field["resName"])

        # Return dictionary of present residues.
        return sequence

    def modelMissingAtoms(self, pdbFilename, outputFilename, chain=' ', debug = False):
        """Model missing atoms/residues in a specified PDB file using MODELLER.
      
        REQUIRED ARGUMENTS
          pdbFilename - the filename of the PDB file to model missing atoms and residues for
          outputFilename - the filename for the desired final model

        OPTIONAL ARGUMENTS
          chain - the one-character chain ID of the chain to model (default ' ')
          debug - flag to print extra debug output and leave temporary directory (default False)

        NOTES

        The specified chain from pdbFilename is processed through MODELLER to build missing
        atoms and residues specified in the SEQRES entry of the PDB file but not present in
        the PDB file.
        
        This procedure is loosely based on the protocol appearing at
        
        http://salilab.org/modeller/wiki/Missing_residues
        
        The complete sequence is read from the SEQRES fields, and the DBREF field used to
        determine the span of residues described in the SEQRES fields.  A heavy-atom topology
        as constructed in MODELLER for the complete sequence, coordinates present in the PDB file
        transferred, and the remaining heavy-atom coordinates built from ideal geometry.
        Finally, a single standard simulated-annealing-based modeling step is performed using
        the standard automodel protocol but allowing only the atoms and residues that were undefined in
        the PDB file to move.
        
        """
        
        # Ensure specified PDB file exists.
        import os.path
        if not os.path.exists(pdbFilename):
            raise ParameterException, "Specified PDB file %s not found." % pdbFilename
        
        # Append full path to pdbFilename and outputFilename
        import os.path
        pdbFilename = os.path.abspath(pdbFilename)
        outputFilename = os.path.abspath(outputFilename)
                
        # Create a temporary directory for running MODELLER.
        import tempfile
        import os.path
        tmpdir = tempfile.mkdtemp()
        if debug: print "tmpdir = %s" % tmpdir
        
        # Get the complete sequence without chain breaks from the SEQRES/DBREF fields of the source PDB file.
        first_residue_id, complete_sequence = self.getCompleteSequence(pdbFilename, chain)
        nresidues = len(complete_sequence)
        last_residue_id = first_residue_id + nresidues - 1
        
        # Get the sequence of residues that are at least partially present in the PDB file as a dictionary.
        # present_sequence_dict[residue_id] is the one-letter-code of the residue residue_id, if there are any ATOM records for this residue.
        present_sequence_dict = self.getPresentSequence(pdbFilename, chain)
                
        # Generate alignment of the template sequence (residues for which any coordinates are defined) against the target (complete sequence from SEQRES/DBREF)
        present_sequence = ""
        for residue_id in range(first_residue_id, first_residue_id + nresidues):
            if present_sequence_dict.has_key(residue_id):
                # TODO: Check integrity against complete_sequence.            
                present_sequence += present_sequence_dict[residue_id]
            else:
                present_sequence += '-'

        # Change working directory to temporary directory.
        import os
        olddir = os.getcwd()
        os.chdir(tmpdir)

        # Generate alignment file for MODELLER.
        import os
        alignment_filename = os.path.join(tmpdir, 'model.ali')
        alignment_file = open(alignment_filename, 'w')
        print >> alignment_file, ">P1;%s" % "template"
        print >> alignment_file, "%s:%s:%d:%s:%d:%s:%s:%s:%s:%s" % ( "structure", pdbFilename, min(present_sequence_dict.keys()), chain, max(present_sequence_dict.keys()), chain, " ", " ", " ", " " )
        print >> alignment_file, "%s*" % present_sequence
        print >> alignment_file, ""    
        print >> alignment_file, ">P1;%s" % "target"
        print >> alignment_file, "%s:%s:%d:%s:%d:%s:%s:%s:%s:%s" % ( "sequence", "target", first_residue_id, chain, last_residue_id, chain, " ", " ", " ", " " )
        print >> alignment_file, "%s*" % complete_sequence
        alignment_file.close()
        if debug:
            import commands
            print "alignment file:"
            print commands.getoutput('cat %(alignment_filename)s' % vars())
        
        # Call MODELLER to generate topology, transfer coordinates, and build from internal coordinates.
        import modeller
        import modeller.automodel
        
        # Create a new environemnt.
        env = modeller.environ()
        
        # Specify the topology and parameters to use.
        # TODO: Is this necessary, or can we rely on the defaults?
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        
        # Read in alignment.
        aln = modeller.alignment(env)
        print alignment_filename
        aln.append(file=alignment_filename, align_codes='all')
        
        # Create a model.
        model = modeller.model(env)
        
        # Generate the topology from the target sequence.
        model.generate_topology(aln['target'])
        
        # Transfer defined coordinates from template.
        model.transfer_xyz(aln)
        
        # Determine which atoms are undefined because they are missing in the template, and create a selection from them.
        missing_atom_indices = []
        for atom_index in range(len(model.atoms)):
            atom = model.atoms[atom_index]
            if atom.x == -999:
                missing_atom_indices.append(atom_index)
                
        # DEBUG: Write model coordinates to a PDB file.
        model.write(file=os.path.join(tmpdir,'transferred.pdb'))
        
        # Build the remaining undefined atomic coordinates from ideal internal coordinates stored in residue topology files.
        model.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
        
        # DEBUG: Write model coordinates to a PDB file.
        if debug: model.write(file=os.path.join(tmpdir,'built.pdb'))
        
        # Override the 'select_atoms' routine in the 'automodel' class to select only the atoms with undefined atomic coordinates in template PDB.
        class mymodel(modeller.automodel.automodel):
            def select_atoms(self):
                missing_atoms = modeller.selection()
                for atom_index in missing_atom_indices:
                    missing_atoms.add(self.atoms[atom_index])
                return missing_atoms

        # Ensure selected atoms feel all nonbonded interactions.
        env.edat.nonbonded_sel_atoms = 1
        
        # Set up automodel.
        #a = mymodel(env, inifile='built.pdb', alnfile=alignment_filename, knowns='template', sequence='target')
        a = mymodel(env, alnfile=alignment_filename, knowns='template', sequence='target')
        
        # Set parameters for automodel.
        # Build only one model.
        # TODO: Have more models built by default (perhaps 50?)
        a.starting_model = 1
        a.ending_model = 1
        
        # Generate model(s).
        a.make()

        # TODO: Rescore models and select the best one.
        # For now, we only use the first model.
        final_model_summary = a.outputs[0]
        
        # Copy resulting model to desired output PDB filename.
        import shutil
        shutil.copy(final_model_summary['name'], outputFilename)
                
        # Restore working directory.
        os.chdir(olddir)
        
        # Clean up temporary directory.
        if (not debug):
            for filename in os.listdir(tmpdir):
                os.remove(os.path.join(tmpdir,filename))
            os.rmdir(tmpdir)

        return

    def makeModel(self, pdbFilename, sequenceFilename, outFile, projName=' ', chain=' ', captermini=False):
        """Creates a PDB file for the sequence by fitting it to the PDB file using Modeller.

        myModel.make_model(pdbFile, sequenceFilename, outFile)
    
        REQUIRED ARGUMENTS
          pdbFilename - the filename of the PDB file to fit to
          sequenceFilename - the filename of the sequence to model, should contain one sequence using one letter code without spaces
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
        if not os.path.exists(sequenceFilename):
            raise ParameterException, "Specified sequence file %s not found." % sequenceFilename
    
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
        completeSequence = self.getCompleteSequenceSimple(pdbFilename, chain)
        nResidues = len(completeSequence)
        print "PDB sequence: " + completeSequence

        # read in the sequence
        seqFile = open(sequenceFilename, 'r')
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
	
	
        if (captermini):
	    class mymodel (automodel):
		def special_patches(self, aln):
		    # Acylate the N-terminal 
		    self.patch(residue_type='ACE ', residues=(self.residues[0])) 				
		    self.patch(residue_type='CT2 ', residues=(self.residues[-1])) 				
	

        # Create a new environemnt.
        env = environ()

        # Override automatic patching for NTER an CTER 
        if (captermini):
	        env.patch_default = False

        # do alignment
        aln = alignment(env)
        mdl = model(env, file=pdbFilename)
        aln.append_model(mdl, align_codes='template', atom_files=pdbFilename)
        aln.append(file='target.ali', align_codes='target')
        aln.align2d()
        aln.write(file='alignment.ali', alignment_format='PIR')
        aln.write(file='alignment.pap', alignment_format='PAP')

        # make a single model
        if (captermini):
            a = mymodel(env, alnfile='alignment.ali', knowns='template', sequence='target')
        else:
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
      
