import tempfile


"""Tools for performing common cleaning operations on pdb files prior to setting up molecular simulations -- such as removing cosolvents, ligands, and waters.


REQUIREMENTS:

    None so far

TODO:

    Nothing so far

AUTHORS:

    D. Mobley, 8/13/2009

CHANGELOG: 

    Instituted 8/13/2009.

"""


def removeHETATM( inpdb, outpdb, retainwater = None, retainresnames = [], retainresnumbers = []):
    """Open specified input pdb file, remove all HETATM coordinate and CONECT entries, and write to specified output PDB file. Optional water manipulation options provided.

ARGUMENTS:
    - inpdb: Pathname of input pdb file
    - outpdb: Pathname of output pdb file
OPTIONAL:
    - retainwater: Default None (remove all water (in HETATM section)). Also may be "True", which will retain all water, or a list specifying the numbers (integer format) (as in pdb file residue column) of water molecules to retain.
    - retainresnames: List of names of residues to retain; default empty list (none)
    - retainresnumbers: List of residue numbers to retain; default empty list (none)

"""

    #Read input from file
    file = open( inpdb, 'r')
    lines = file.readlines()
    file.close()

    #Storage for output
    outpdbtext = []
    #Store atoms to remove
    removeatoms = []


    for line in lines:
        fieldtype = line[0:6].split()[0]
        #If we're in the coordinate section
        if fieldtype=='ATOM' or fieldtype=='HETATM':
            #Extract some other info
            tresname = line[17:20].split()[0] #Residue name
            taltloc = line[16:17] #alt loc indicator
            tresnum = int(line[22:26].split()[0]) #residue number
            tatom = int(line[6:11].split()[0]) #Atom number
            tchain = line[21] #Chain ID

            #Handle water
            if tresname =='HOH':
                if retainwater: #If it is True or a list
                    if retainwater==True:
                        outpdbtext.append(line)
                    elif retainwater.count( tresnum ) >0:
                        outpdbtext.append(line)
                    else:
                        removeatoms.append( tatom ) #Track which atoms we are removing
                else:
                    removeatoms.append( tatom ) #Track atoms we're removing
            #Handle all other HETATM types
            elif fieldtype == 'HETATM':
                #Retain residue names or numbers if specified
                if retainresnumbers.count(tresnum)>0 or retainresnames.count(tresname)>0:
                    outpdbtext.append(line)
                else: #Otherwise remove
                    removeatoms.append( tatom )
            else: #If it is an ATOM entry, store it
                outpdbtext.append(line)

        elif fieldtype=='CONECT':  #Remove CONECT entries for specified atoms
            anums = line.replace('CONECT','').split()
            anums_i = [int(atom) for atom in anums]
            found = False
            for atom in anums_i:
                if removeatoms.count(atom)>0:
                    found = True
            if not found: #Store only if this is not a connect entry involving some of the atoms we are removing
                outpdbtext.append(line)


        else: #If it is not ATOM/HETATM/CONECT, store.
            outpdbtext.append(line)


    file = open(outpdb, 'w')
    file.writelines(outpdbtext)
    file.close()

def removeResidue( inpdb, outpdb, resnum = None, resname = None):
    """Remove a specified residue from an input pdb file, write output to output pdb file.

ARGUMENTS:
    - inpdb: Pathname of input pdb file
    - outpdb: Pathname of output pdb file
AND ONE OF:
    - resnum: Residue number to remove; default None
    - resname: Residue name to remove; default None
First checks for matches to resnum; if resnum is not specified, checks for resname. Specifying both will result in only resnum being used.

"""


    file = open(inpdb, 'r')
    lines = file.readlines()
    file.close()

    removeatoms = []
    outpdbtext=[]

    for line in lines:
        fieldtype = line[0:6].split()[0]
        if fieldtype=='HETATM' or fieldtype=='ATOM':
            tresname = line[17:20].split()[0]
            taltloc = line[16:17]
            tresnum = int(line[22:26].split()[0])
            tatom = int(line[6:11].split()[0])
            tchain = line[21]
            match = False
            
            #Look to see if it matches
            if resnum:
                if resnum.count(tresnum)>0:
                    match = True
            elif resname:
                if resname.count(tresname)>0:
                    match = True
            if not match:
                outpdbtext.append(line)
            else:
                removeatoms.append( tatom )
    
        #Remove corresponding CONECT entries
        elif fieldtype=='CONECT':  #Remove CONECT entries for specified atoms
            anums = line.replace('CONECT','').split()
            anums_i = [int(atom) for atom in anums]
            found = False
            for atom in anums_i:
                if removeatoms.count(atom)>0:
                    found = True
            if not found: #Store only if this is not a connect entry involving some of the atoms we are removing
                outpdbtext.append(line)

        #Store other lines
        else:
            outpdbtext.append(line)


    #Write output
    file = open(outpdb, 'w')
    file.writelines( outpdbtext)
    file.close()
