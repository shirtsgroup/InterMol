
# ----------------------------------------------------------------------
# MODULE DATA
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# METHODS
# ----------------------------------------------------------------------

def one_letter_code(three_letter_code):
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

def getCompleteSequence(pdbfilename, chain=' '):
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
            # Store initial residue.
            first_residue = field["seqBegin"]
    # Raise an exception if the PDB file does not contain a DBREF field.
    if not first_residue:
        raise RuntimeError, "No DBREF field found in PDB file -- noncompliant with PDB 2.1 format."
            
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
                    complete_sequence += one_letter_code(resName)

    # DEBUG
    print "complete sequence: %s" % complete_sequence

    # Return first residue and complete sequence.
    return first_residue, complete_sequence

def getPresentSequence(pdbfilename, chain=' '):
    """Extract the sequence for which there are atomic coordiantes defined from a PDB file.

    present_sequence = getPresentSequence(pdbfilename, chain=' ')
    
    REQUIRED ARGUMENTS
      pdbfilename - the filename of the PDB file to import from

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
                sequence[field["resSeq"]] = one_letter_code(field["resName"])

    # Return dictionary of present residues.
    return sequence

def readAtomsFromPDB(pdbfilename):
    """Read atom records from the PDB and return them in a list.

    present_sequence = getPresentSequence(pdbfilename, chain=' ')
    
    REQUIRED ARGUMENTS
      pdbfilename - the filename of the PDB file to import from

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

    # Read atoms.
    atoms = [ ]
    for line in lines:
        if line[0:6] == "ATOM  ":
            # Parse line into fields.
            atom = { }
            atom["serial"] = int(line[6:11])
            atom["name"] = line[12:16].strip()
            atom["altLoc"] = line[16:17]
            atom["resName"] = line[17:20]
            atom["chainID"] = line[21:22]
            atom["resSeq"] = int(line[22:26])
            atom["iCode"] = line[26:27]
            atom["x"] = float(line[30:38])
            atom["y"] = float(line[38:46])
            atom["z"] = float(line[46:54])
            atom["occupancy"] = float(line[54:60])
            atom["tempFactor"] = float(line[60:66])
            atom["segID"] = line[72:76]
            atom["element"] = line[76:78]
            atom["charge"] = line[78:80]
            
            atoms.append(atom)
            
    # Return list of atoms.
    return atoms

# Model missing atoms/residues using MODELLER.
def modelMissingAtoms(pdbfilename, chain=' '):
    """Model missing atoms/residues of current selection using MODELLER.
      
    REQUIRED ARGUMENTS
      pdbfilename - the filename of the PDB file to import from

    OPTIONAL ARGUMENTS
      chain - the one-character chain ID of the chain to import (default ' ')

    The specified chain from pdbfilename is processed through MODELLER to build missing
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
    the PDB file to move.  The resulting structure is loaded into PyMOL.
    
    """

    # Ensure file exists.
    import os.path
    if not os.path.exists(pdbfilename):
        raise ParameterException, "Specified PDB file %s not found." % pdbfilename

    # Append full path to pdbfilename.
    import os.path
    pdbfilename = os.path.abspath(pdbfilename)

    # Create a temporary directory for running MODELLER.
    import tempfile
    import os.path
    tmpdir = tempfile.mkdtemp()
    print "tmpdir is: %s" % tmpdir

    # Get the complete sequence.
    first_residue_id, complete_sequence = getCompleteSequence(pdbfilename, chain)
    nresidues = len(complete_sequence)
    last_residue_id = first_residue_id + nresidues - 1

    # Get the sequence that is present in the PDB file as a dictionary.
    # present_sequence_dict[residue_id] is the one-letter-code of the present residue, if present
    present_sequence_dict = getPresentSequence(pdbfilename, chain)

    # Generate alignment of the present sequence against the target.
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
    print >> alignment_file, "%s:%s:%d:%s:%d:%s:%s:%s:%s:%s" % ( "structure", pdbfilename, min(present_sequence_dict.keys()), chain, max(present_sequence_dict.keys()), chain, " ", " ", " ", " " )
    print >> alignment_file, "%s*" % present_sequence
    print >> alignment_file, ""    
    print >> alignment_file, ">P1;%s" % "target"
    print >> alignment_file, "%s:%s:%d:%s:%d:%s:%s:%s:%s:%s" % ( "sequence", "target", first_residue_id, chain, last_residue_id, chain, " ", " ", " ", " " )
    print >> alignment_file, "%s*" % complete_sequence
    alignment_file.close()

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

    # Get the remaining undefined coordinates from internal coordinates.
    model.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

    # DEBUG: Write model coordinates to a PDB file.
    model.write(file=os.path.join(tmpdir,'built.pdb'))
    
    # Override the 'select_atoms' routine in the 'automodel' class
    class mymodel(modeller.automodel.automodel):
        def select_atoms(self):
            missing_atoms = modeller.selection()
            for atom_index in missing_atom_indices:
                missing_atoms.add(self.atoms[atom_index])
            return missing_atoms

    # Ensure selected atoms feel all nonbonbded interactions.
    env.edat.nonbonded_sel_atoms = 1

    # Set up automodel.
    #a = mymodel(env, inifile='built.pdb', alnfile=alignment_filename, knowns='template', sequence='target')
    a = mymodel(env, alnfile=alignment_filename, knowns='template', sequence='target')

    # Set parameters for automodel.
    a.starting_model = 1
    a.ending_model = 1

    # Generate models.
    a.make()

    # Get the resulting model.
    final_model_summary = a.outputs[0]
    
    # Copy resulting model to 'final.pdb'.
    import shutil
    shutil.copy(final_model_summary['name'], 'final.pdb')

    # Restore working directory.
    os.chdir(olddir)
    
    # Clean up temporary directory.
    for filename in os.listdir(tmpdir):
        os.remove(os.path.join(tmpdir,filename))
    os.rmdir(tmpdir)

