# hacks.py 
#
# An protein topology is constructed in the extended conformation in AMBER LEaP, and an implicit solvent simulation
# is conducted at high temperature with a Langevin thermostat to ensure rapid sampling of thermally denatured states.
#
# We start with a PDB file that comes from MCCE with the appropriate AMBER residues names already substituted.
#
# John D. Chodera
# 6 Jun 2007
#
# TODO:
# - Determine how to use substitution with both class methods and local variables.

def getPresentSequence(pdbfilename, chain=' ', allow_four_letters = False):
    """Extract the sequence for which there are atomic coordiantes defined from a PDB file.

    present_sequence = getPresentSequence(pdbfilename, chain=' ')
    
    REQUIRED ARGUMENTS
      pdbfilename - the filename of the PDB file to import from

    OPTIONAL ARGUMENTS
      chain - the one-character chain ID of the chain to import (default ' ')
      allow_four_letters (logical) - if True, will allow four-letter amino acid names (default: False)

    RETURN VALUES
      sequence (dictionary) - sequence[residue_id] is the one-letter code corresponding to residue index residue_id

    The ATOM records are read, and the sequence for which there are atomic coordinates is stored.

    """

    # Read the PDB file into memory.
    pdbfile = open(pdbfilename, 'r')
    lines = pdbfile.readlines()
    pdbfile.close()

    # Extract the sequence for which there are defined atomic coordinates.
    sequence = ""
    last_resSeq = None
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
            if (chain == field['chainID'] and field['resSeq'] != last_resSeq):
                if(len(sequence) > 0):
                    sequence += " " # add space before three-letter code
                sequence += field["resName"]
                last_resSeq = field['resSeq']
                # Add fourth letter, if present.
                if(allow_four_letters and line[20] != ' '):
                    sequence += line[20]

    # Return dictionary of present residues.
    return sequence

def getCompleteSequence(pdbFilename, chain=' '):
    """Extract the complete three-letter-code sequence from a PDB file.

    sequence = getCompleteSequence(pdbfilename, chain=' ')
    
    REQUIRED ARGUMENTS
        pdbFilename - the filename of the PDB file to import from
    
    OPTIONAL ARGUMENTS
        chain - the one-character chain ID of the chain to import (default ' ')

    RETURN VALUES
        sequence (string) - the complete sequence in space-separated three-letter code

    """
    
    # Read the PDB file into memory.
    pdbFile = open(pdbFilename, 'r')
    lines = pdbFile.readlines()
    pdbFile.close()
       
    # Get the complete sequence from the SEQRES fields.
    completeSequence = ""  # sequence in concatenated space-delimited three-letter codes
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
                    if(len(completeSequence) > 0):
                        completeSequence += " " # add space before three-letter code
                    completeSequence += resName # add three-letter code
                    
    # Return first complete space-separated three-letter code sequence
    return completeSequence
