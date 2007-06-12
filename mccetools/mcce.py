import tempfile
import os
import sys
import commands
import re
import math

"""
Calculation of residue pKa and likely protonation states using the MCCE (multiconformer continuum electrostatics) package.

DESCRIPTION
  This package contains various tools for using MCCE to predict protonation states and set up output pdb file using these protonation states (and Pande group GROMACS AMBER port naming conventions) for simulation. Originally written in Feb./Mar. 2007 by Imram Haque; D. Mobley edited in May 2007. Variously edited by John D. Chodera, Vincent Voelz, and Greg Bowman thereafter.

REQUIREMENTS
- MCCE (available from http://www.sci.ccny.cuny.edu/~mcce)
- MCCE requires a working copy of the continuum electrostatics package 'delphi' from the Honig group.

FUNCTIONALITY
- Main function seems to be protonation_state, which takes as input a pdb file, a pH, and the path to MCCE, and returns an array of lines of a pdb file, generated using predicted protonation states. Currently seems to require absolute path to the input pdb file. Reads input PDB, runs MCCE, which protonates the pdb and generates various temporary files. Calls the remaining functions in here to parse the MCCE output and put the appropriate protonation states into the output pdb, which is returned as a text array; make naming in output pdb be consistent with AMBER naming. Also checks histidine protonation and names histidines accordingly; checks for disulfide bonds and renames; gives terminal residues correct names. Note that this can be rather slow (hours+?) due to MCCE speed and depending on size of protein. Note also that the returned pdb file is a text array without newline characters. 

EXAMPLES

KNOWN LIMITATIONS
- C terminal residues are determined by the last residue in the pdb file, which means the last residue in the file needs to not be a HET atom or similar, but rather the terminal residue of the protein.

KNOWN BUGS
- There may be a bug with the determination of protonation states because a Lysine-rich 13-residue N-terminal segment of 3gb1 came out entirely neutral.

TODO
- We should be sure to strip out waters and possibly HETATMs first, though we are unlikely to get these from MODELLER.
- Move methods that rename residues to match AMBER convention based on protonation state to ambertools.
- Make some of these functions private if they are just going to confuse users when they do help(mcce).
- Get default MCCE_HOME from environment variables, if present.
- Add class constructor to store MCCE path and prm file contents internally and run a self-test to make sure it is functioning correctly.  Turn existing methods into class methods rather than static methods.
- Update according to current Pande lab style guide.
- Add examples
- Do something less confusing with the parameter dictionary.
- Test on proteins with known pKas
- Select most likely global protonation state, rather than on a per-residue basis.

REVISION LOG
- 6-06-2007: JDC: Made to conform more closely to Pande group style guide.  Added/modified documentation.  Fixed bug in doubly-protonated HIS renaming.
- 6-05-2007: VAV;  Fixed xtraprms={} bugs
- 5-21-2007: JDC: Updating according to current Pande lab style guide.
- 5-8-2007: DLM: Adding preliminary documentation (above) based on my known knowledge and perusing some of the below. Lots more needs to be done.
- 5-9-2007: DLM: A bit more documentation, and fixed protonation_state routine so as not to require an absolute path to the pdb file (obtains an absolute path from the file name plus current directory).
- 5-9-2007: DLM: Fixed a bug in oxygen naming on the C terminal residue (used O and OXT, though it was supposed to (and claimed to) use OC1 and OC2.)

- 5-14-2005: VV: Fixed another terminii naming problem -->  ffamber residues do not include a CLYN (C-terminal neutral lysine), only CLYP.
- 5-15-2005: VV: Added protonatePDB, which writes a pdb outputfile file with the right protomation state
-                Added read_paramfile().  paramgen() can now also take default params from a MMCE-style prm file

"""

def paramgen(pdbpath, mcce_location, fromfile=None, xtraprms={}):
    """Generate a dictionary of parameters used by MCCE.
    
    REQUIRED ARGUMENTS
        pdbpath - The path to the PDB file used as input for MCCE.  Must be ABSOLUTE path  
        mcce_location - The base directory in which MCCE is installed.

    OPTIONAL ARGUMENTS
        fromfile - if provided, parameters are read from the given MCCE parameter file.  Note that MCCE location and PDB path parameters are overridden.   
                
    RETURNS
        A dictionary of MCCE parameters
                   
    TODO
    - Replace hard-coded defaults with default parameter filename.
    
    """
    
    # Store MCCE parameters in a dictionary.
    mcce_params={} 

    # This was just generated from the sample .prm files in the MCCE distribution.
    # TODO: Load these in automatically from a default .prm file instead of having them hardcoded in here.
    mcce_params['DO_PREMCCE']='t'
    mcce_params['DO_ROTAMERS']='t'
    mcce_params['DO_ENERGY']='t'
    mcce_params['DO_MONTE']='t'
    mcce_params['EPSILON_PROT']='8.0'
    mcce_params['TITR_TYPE']='ph'
    mcce_params['TITR_PH0']='7.0'
    mcce_params['TITR_PHD']='1.0'
    mcce_params['TITR_EH0']='0.0'
    mcce_params['TITR_EHD']='30.0'
    mcce_params['TITR_STEPS']='1'
    mcce_params['MINIMIZE_SIZE']='t'
    mcce_params['TERMINALS']='t'
    mcce_params['CLASH_DISTANCE']='2.0'
    mcce_params['H2O_SASCUTOFF']='0.05'
    mcce_params['ROT_SPECIF']='f'
    mcce_params['ROT_SWAP']='t'
    mcce_params['PACK']='f'
    mcce_params['ROTATIONS']='6'
    mcce_params['SAS_CUTOFF']='1.00'
    mcce_params['VDW_CUTOFF']='10.0'
    mcce_params['REPACKS']='5000'
    mcce_params['REPACK_CUTOFF']='0.01'
    mcce_params['HDIRECTED']='f'
    mcce_params['HDIRDIFF']='1.0'
    mcce_params['HDIRLIMT']='36'
    mcce_params['HV_RELAX_NCYCLE']='0'
    mcce_params['HV_RELAX_DT']='10'
    mcce_params['HV_RELAX_NITER']='10'
    mcce_params['HV_RELAX_VDW_THR']='2.'
    mcce_params['HV_RELAX_HV_VDW_THR']='5.'
    mcce_params['HV_TORS_SCALE']='20.'
    mcce_params['HV_RELAX_N_SHAKE']='10000'
    mcce_params['HV_RELAX_CONSTRAINT']='1.0'
    mcce_params['HV_RELAX_CONSTRAINT_FRC']='20.'
    mcce_params['RELAX_H']='t'
    mcce_params['RELAX_E_THR']='-5.0'
    mcce_params['RELAX_NSTATES']='50'
    mcce_params['RELAX_N_HYD']='6'
    mcce_params['RELAX_CLASH_THR']='5.'
    mcce_params['RELAX_PHI']='1.0'
    mcce_params['RELAX_NITER']='300'
    mcce_params['RELAX_TORQ_THR']='0.5'
    mcce_params['EPSILON_SOLV']='80.0'
    mcce_params['GRIDS_DELPHI']='65'
    mcce_params['GRIDS_PER_ANG']='2.0'
    mcce_params['RADIUS_PROBE']='1.4'
    mcce_params['IONRAD']='2.0'
    mcce_params['SALT']='0.15'
    mcce_params['BIG_PAIRWISE']='5.0'
    mcce_params['MONTE_SEED']='-1'
    mcce_params['MONTE_T']='298.15'
    mcce_params['MONTE_FLIPS']='3'
    mcce_params['MONTE_NSTART']='100'
    mcce_params['MONTE_NEQ']='300'
    mcce_params['MONTE_REDUCE']='0.001'
    mcce_params['MONTE_RUNS']='6'
    mcce_params['MONTE_NITER']='2000'
    mcce_params['MONTE_TRACE']='50000'
    mcce_params['NSTATE_MAX']='1000000'
    mcce_params['DELPHI_START']='1'
    mcce_params['DELPHI_END']='99999'
    
    # Read parameters from a .prm file, if provided.
    if fromfile != None:
        if os.path.exists(fromfile):
            print 'Reading parameters from file',fromfile,'...'
            fileparams = read_paramfile(fromfile)
            for key in fileparams.keys():
                print 'key',key,'fileparams[key]',fileparams[key]
                mcce_params[key] = fileparams[key]
        else:
            print 'Can\'t find parameter file',fromfile,'. Exiting...'
            # TODO: Throw exception here instead.
            exit
                
    # Override PDB- and MCCE-LOCATION-dependent fields.
    mcce_params['MCCE_HOME'] = mcce_location
    mcce_params['INPDB'] = pdbpath
    mcce_params['EXTRA'] = mcce_params['MCCE_HOME']+'/extra.tpl'
    mcce_params['RENAME_RULES'] = mcce_params['MCCE_HOME']+'/name.txt'
    mcce_params['DELPHI_EXE'] = mcce_params['MCCE_HOME']+'/bin/delphi'

    # Enable the four stages of MCCE so that it will perform all necessary steps to sample protonation states.
    # Note that these are typically set to 'f' in the provided run.prm.* template files.
    mcce_params['DO_PREMCCE']  = 't'
    mcce_params['DO_ROTAMERS'] = 't'
    mcce_params['DO_ENERGY']   = 't'
    mcce_params['DO_MONTE']    = 't'
        
    # Override any extra parameters specified.
    for key in xtraprms.keys():
        mcce_params[key]=xtraprms[key]
        
    # Return the updated dictionary of MCCE parameters.
    return mcce_params

def print_prm(mcceparams):
    """Print the contents of the dictionary of MCCE parameters.
        
    ARGUMENTS
        mcceparams - the dictionary of MCCE parameters to be printed
        
    """
    
    sorted = mcceparams.keys()
    print sorted
    sorted.sort()
    for i in sorted:
        print '%-40s%s'%(mcceparams[i]," ("+i+")")
        
    return

def read_paramfile(paramfile):
        """Read the contents of an MCCE parameter file, returning the key-value pairs in a dictionary.
        
        ARGUMENTS
            paramfile - absolute pathname of the MCCE parameter file to be read
            
        RETURNS
            params - a dictionary of the key-value pairs contained in the MCCE parameter file
                
        """
        
        # Create a dictionary to store MCCE parameters.
        params = {}  # params['key'] contains the value corresponding to 'key' in the MCCE file

        # Read and parse the MCCE file.
        fin = open(paramfile,'r')
        for line in fin.readlines():
            fields = line.split()
            if len(fields) > 1:
              if (fields[-1][0] == '(' ) & (fields[-1][-1] == ')'):
                params[ fields[-1][1:-1] ] = fields[0] 
        fin.close()                    
        
        # Return the dictionary of MCCE parameters.
        return params

def write_paramfile(mcceparams, filename):
    """Write MCCE parameters (stored as key-value dictionary) to a specified MCCE parameter file.
    
    ARGUMENTS
        mcceparams - dictionary of MCCE parameters stored as key-value pairs
        filename - name of MCCE parameter file to write to (overwriting if already exists)
    """

    # Open the MCCE parameter file for writing.
    runfile = open(filename,'wt')        

    # Sort keys.
    sorted = mcceparams.keys()
    sorted.sort()

    # Write key-value pairs in MCCE format.
    for key in sorted:
        runfile.write('%-40s (%s)\n'%(mcceparams[key],key) )

    # Close the MCCE parameter file.
    runfile.close()
        
def run_mcce(mcceparams):
    """Run MCCE executable using specified parameter dictionary.
    
    ARGUMENTS
        mcceparams - Dictionary of MCCE parameters
    
    RETURNS 
        output - stdout capture from MCCE run
        
    Runs MCCE in the current working directory, using run parameters given in variable mcceparams (from the paramgen function).
    This method will clobber file run.prm in the current directory, as well as any previous MCCE output files.
    It will generate all sorts of output files, which it is the caller's responsibility to clean up.
    """
        
    os.chdir(os.getcwd())
    write_paramfile(mcceparams, 'run.prm')
    output = commands.getoutput( os.path.join(mcceparams['MCCE_HOME'],'bin/mcce') )
    return output
    
def ps_mostlikely(fort38path):
    """Finds the set of most likely protonation states from the file fort38path/fort.38. Assumes one set of pH values in the given file.

    ARGUMENTS
        fort38path - absolute pathname of the directory in which the fort.38 file produced by MCCE appears
        
    RETURNS
        searchstring - regular expression query string for most likely states
        neutrallines - regular expression query string for neutral states
        poslines - regular expression query string for positive states
        neglines - regular expression query string for negative states
    """
    
    f=open(fort38path+"/fort.38","rt")
    lines=f.readlines()[1:]         #Skip the header
    lines2=lines[:] # Make a copy of file contents because the extraction procedure will destroy this copy.
    f.close()
    stack=[]
    # We want any _000 line - they're common to all protonation states
    searchstring="_000 "
    
    # Overview:
    # Grab a line and compare to the last one.
    # If it's the same residue, keep the more likely one.
    # If a new residue, add the stored line to the master "keep" list and store what we just found.
    while (len(lines)>0):
        nextline=lines.pop(0)
        if (len(stack)==0):
            stack.append(nextline)
        else:
            stackline=stack.pop()
            nextident=nextline[0:3]+nextline[6:10]
            stackident=stackline[0:3]+stackline[6:10]
            if (nextident!=stackident):
                # Add the identifier to the list of atoms to pull
                searchstring=searchstring+"|"+stackline[0:3]+r" "+stackline[5:15]
                #print "Keeping "+stackline[0:15]
                stack.append(nextline)
            else:
                # Next option for the same residue
                stackprob=float(stackline.split(" ")[1])
                nextprob=float(nextline.split(" ")[1])
                if (nextprob > stackprob):
                    stack.append(nextline)
                else:
                    stack.append(stackline)
    # Handle the last element on the stack
    stackline=stack.pop()
    searchstring=searchstring+"|"+stackline[0:3]+r" "+stackline[5:15]

    # Generate regular expressions for neutral/pos/neg states based on charge labeling in our duplicate array.
    # Select the lines for any particular state, then transform each line into a regular expression, then join lines with | operator.
    neutrallines = "_000 |"+("|".join(map(lambda x: x[0:3]+r" "+x[5:15],filter(lambda x:re.search(r'^...0',x),lines2))))
    poslines = "|".join(map(lambda x: x[0:3]+r" "+x[5:15],filter(lambda x:re.search(r'^...\+',x),lines2)))
    neglines = "|".join(map(lambda x: x[0:3]+r" "+x[5:15],filter(lambda x:re.search(r'^...-',x),lines2)))

    # Return regular expression search strings.
    return (searchstring,neutrallines,poslines,neglines)        
        
def renumber_atoms(pdbarr):
    """Renumbers atom entries (list items) of a PDB file so they are in sequence starting from 1.
    
    ARGUMENTS
        pdbarr - list of PDB atoms to be renumbered; mutated by call
    """
    
    # Renumber the atom entries in place.
    for i in range(len(pdbarr)):
        pdbarr[i]=pdbarr[i][0:6]+str(i+1).rjust(5)+pdbarr[i][11:]        

def rename_termini(npdb):
    """Renames the 'NTR' and 'CTR' residues generated by MCCE.
    Renames the O and OXT atoms at the C-terminus to OC1 and OC2 for AMBER.
    
    ARGUMENTS
        npdb - 
    
    """

    # Get three-letter residue names of N- and C-termini.
    ntrname = npdb[1][0][17:20]
    ctrname = npdb[-2][0][17:20]

    # Rename N-terminal residue 'XXX' to 'NXXX'.
    for j in [0,1]:
        for i in range(len(npdb[j])):
            npdb[j][i] = npdb[j][i][0:17]+'N'+ntrname+npdb[j][i][21:]
            
    # Rename C-terminal residue from 'XXX' to 'CXXX', with some exceptions (noted below).
    # Also rename terminal oxygen atoms.
    for j in [-1,-2]:
        for i in range(len(npdb[j])):
            # VAV 5/14/2007: gmx ffamber does not have CLYN, only CLYP 
            if ctrname=='LYN':      
                npdb[j][i] = npdb[j][i][0:17]+'CLYP'+npdb[j][i][21:]
            else:   
                npdb[j][i] = npdb[j][i][0:17]+'C'+ctrname+npdb[j][i][21:]
                            
            # DLM 5/9/2007: Though this says that it changes the O and OXT names, it doesn't
            # Adding the below to do the renaming.
            if npdb[j][i][13:16].split()[0]=='O':
                npdb[j][i] = npdb[j][i][0:13]+'OC1'+npdb[j][i][16:]
            elif npdb[j][i][13:16].split()[0]=='OXT':
                npdb[j][i] = npdb[j][i][0:13]+'OC2'+npdb[j][i][16:]
                
            #VAV 6/5/2007: gmx ffamber doesn't have "CD1" as an ILE/CILE atomname, only "CD"
            if ctrname=='ILE':      
                if npdb[j][i][13:16].split()[0]=='CD1':
                    npdb[j][i] = npdb[j][i][0:13]+'CD '+npdb[j][i][16:]
        
def nest_pdb(pdbarr):
    """Collect lines of the PDB file by residue.
    
    ARGUMENTS:
        pdbarr - list of lines from PDB file
        
    RETURNS
        nestedpdb - nested PDB file, such that nestedpdb[i][j] will be line j from residue i.
    """
    
    nestedpdb=[]
    residue=[]
    for line in pdbarr:
        if (len(residue) == 0):
            residue.append(line)
        else:
            if (line[17:27] == residue[-1][17:27]):
                residue.append(line)
            else:
                nestedpdb.append(residue)
                residue=[line]
    nestedpdb.append(residue)
    
    # Return the nexted PDB file.
    return nestedpdb

def unnest_pdb(npdb):
    """Expand the lines from the "nested" PDB file into a flat list of lines.
    
    ARGUMENTS
        npdb - "nested" PDB file; mutated in place
        
    RETURNS
        pdbarr - un-nested PDB file, such that pdbarr[i] will be line i from the PDB file
    """
    
    pdbarr=[]
    for res in npdb:
        for atm in res:
            pdbarr.append(atm)
    return pdbarr
    
def rename_residues(pdbarr):
    """Convert residue names (and some heavy atom names) to AMBER (ffamber) format.
    
    ARGUMENTS
        pdbarr - array of PDB atoms
        
    RETURNS
        pdbarr - updated array of PDB atoms with AMBER appropriate residue names
    """

    # DEBUG check.
    #print "At start of renaming, pdbarr has",len(pdbarr),"items"

    # Transform the PDB file into a nested format that is easier to manipulate.
    npdb = nest_pdb(pdbarr)
    
    # DEBUG check.
    #print "Nest/unnest leads to",len(unnest_pdb(npdb)),"items"

    # Convert residue and atom names to be appropriate for AMBER.
    rename_charged(npdb)
    histidine_search(npdb)
    disulfide_search(npdb)
    rename_termini(npdb)

    # Un-nest PDB
    pdbarr = unnest_pdb(npdb)
    
    # DEBUG check.
    #print "At end of renaming, pdbarr has ", len(pdbarr), "items"

    # Return the PDB file.
    return pdbarr

def rename_charged(npdb):
    """Generate AMBER-specific residue names for charged residues from MCCE residue names.
        
    ARGUMENTS
        npdb - "nested" PDB data structure generated by nest_pdb. Modified to reflect AMBER names.
    """
        
    resname_and_state = map(lambda x:x[-1][17:20]+x[-1][-1],npdb)
    #print resname_and_state
    for i in range(len(npdb)):
        if (resname_and_state[i]=='HIS+'):
            npdb[i]=map(lambda x:x.replace('HIS','HIP'),npdb[i])
        elif (resname_and_state[i]=='LYS0'):
            npdb[i]=map(lambda x:x.replace('LYS','LYN'),npdb[i])
        elif (resname_and_state[i]=='LYS+'):
            npdb[i]=map(lambda x:x.replace('LYS','LYP'),npdb[i])
        elif (resname_and_state[i]=='CYS0'):
            npdb[i]=map(lambda x:x.replace('CYS','CYN'),npdb[i])
        elif (resname_and_state[i]=='CYS-'):
            npdb[i]=map(lambda x:x.replace('CYS','CYM'),npdb[i])
        elif (resname_and_state[i]=='ASP0'):
            npdb[i]=map(lambda x:x.replace('ASP','ASH'),npdb[i])
        #Aspartate requires no sub
        elif (resname_and_state[i]=='GLU0'):
            npdb[i]=map(lambda x:x.replace('GLU','GLH'),npdb[i])
        # Glutamate requires no sub

    return

def get_coords(atomname,residue):
    """Extract coordinates from given atom in list of PDB lines (such as all those lines for a residue in a nested PDB file).
    
    ARGUMENTS
        atomname - the name of the atom
        residue - list of lines from PDB file

    RETURNS
        (x,y,z) - tuple of atom coordinates, as floats
    
    """
    for n in range(len(residue)):
        if (residue[n][12:16].strip() == atomname.strip()):
            iX=float(residue[n][30:38])
            iY=float(residue[n][38:46])
            iZ=float(residue[n][46:54])
            return (iX,iY,iZ)

    raise "Atom not found!"        

def disulfide_search(npdb, min_dist = 1.8, max_dist = 2.2):
    """Rename CYS to CYX if it participates in a disulfide bond, as judged by distance range (inclusive).
    
    ARGUMENTS
        npdb - nested PDB file
        
    OPTIONAL ARGUMENTS
        min_dist - minium distance cutoff for perceiving disulfide bond (default 1.8 A)
        max_dist - maximum distance cutoff for perceiving disulfide bond (default 2.2 A)
    """
        
    residues_to_rename=set([])
    for i in range(len(npdb)):
        if (npdb[i][0][17:20] != 'CYS'):
            continue
        # Found a cysteine, now track down the sulfur
        iX,iY,iZ=get_coords('SG',npdb[i])
                
        for j in range(i+1,len(npdb)):
            if (npdb[j][0][17:20]!= 'CYS'):
                continue
            (jX,jY,jZ) = get_coords('SG',npdb[j])
                
            dX=iX-jX
            dY=iY-jY
            dZ=iZ-jZ

            distance = math.sqrt(dX*dX+dY*dY+dZ*dZ)
            if (distance >= min_dist and distance <= max_dist):
                residues_to_rename = residues_to_rename | set([i,j])
        
    # Rename the residues we selected
    for i in residues_to_rename:
        npdb[i]=map(lambda x:x.replace('CYS','CYS2'),npdb[i])           
        
    return

def atom_is_present(pdblines, atomname):
    """Returns TRUE if the given atom is present in the given PDB atom lines.
    
    ARGUMENTS
        pdblines - list of PDB lines
        atomname - the name of the atom to check the existence of

    RETURNS
        is_present - True if the given atom name is present, False otherwise

    """

    is_present = False
    for pdbline in pdblines:
        if (pdbline[13:16] == atomname):
            is_present = True

    return is_present

def histidine_search(npdb):
    """Rename HIS residues to HID, HIE, or HIP by examining which protons are present.
    
    ARGUMENTS
        npdb - nested  PDB
    """
    
    for i in range(len(npdb)):
        resname = npdb[i][0][17:20]
        if (resname == 'HIS'):
            HE_present = atom_is_present(npdb[i], 'HE2') # bonded to NE2
            HD_present = atom_is_present(npdb[i], 'HD1') # bonded to ND1
            
            if (HD_present and HE_present):
                npdb[i]=map(lambda x:x.replace('HIS','HIP'),npdb[i])                
            elif (HE_present):
                npdb[i]=map(lambda x:x.replace('HIS','HIE'),npdb[i])
            elif (HD_present):
                npdb[i]=map(lambda x:x.replace('HIS','HID'),npdb[i])
            else:
                raise "No protons found for histidine."

    return

def pdb_cleanup(pdbarr):
    """Remove extraneous entries from MCCE PDB files.
    
    ARGUMENTS
        pdbarr - list of PDB lines
        
    RETURNS
        updated_pdbarr - updated list of PDB lines
    """
    # Nuke everything after the coordinates
    # Use a default occupancy of 1 and a default B-factor of 0
    for i in range(len(pdbarr)):
        atom=pdbarr[i][0:54]
        # The atom symbol will be the first character after stripping
        # whitespace and numbers
        # NOTE: this will fail for atoms with more than one letter in their symbol
        #            eg, counterions
        atomsymbol = atom[12:16].strip(" 0123456789")[0]
        atom = atom + "1.00".rjust(6) + "0.00".rjust(6) + " "*10 + atomsymbol.rjust(2)+ " "*2
        pdbarr[i]=atom

    npdb=nest_pdb(pdbarr)

    for i in range(len(npdb)):
        # Set the residue index numbers
        # NOTE: We set chain identifier to blank, so this won't work for multi-chain PDBs
        for j in range(len(npdb[i])):
            atom = npdb[i][j]
            atom = atom[:21]+ " " + str(i+1).rjust(4)+ " "*4 + atom[30:]
            #print len(atom)
            npdb[i][j] = atom
            
    updated_pdbarr = unnest_pdb(npdb)
        
    return updated_pdbarr

def protonation_state(pdbfile, pH, mccepath, cleanup=True, prmfile=None, xtraprms={}):
    """Performs a pH titration on all titratable residues in a PDB file."""
    
    #Convert PDB file name to path, as paramgen expects an absolute path
    pdbpath=os.path.join(os.getcwd(),pdbfile)
    print 'pdbpath', pdbpath
    
    # Set up the parameters for the MCCE run
    params=paramgen(pdbpath,mccepath, fromfile=None, xtraprms=xtraprms)
    params['TITR_PH0']=str(pH)
    params['TITR_PHD']="1.0" #pH titration interval (i.e. stepsize)
    params['TITR_STEPS']="1"
        
    # Create a temporary directory with the run.prm file and run MCCE
    tempdir=tempfile.mkdtemp();
    print "Running MCCE in temporary directory %s..." % tempdir

    # remember our curent working directory
    thisdir = os.getcwd()
        
    # run mcce in the temporary directory
    os.chdir(tempdir)
    output = run_mcce(params)
    print output
        
    # chdir back to where we came from
    os.chdir(thisdir)

    pdbarr = ps_processmcce(tempdir)

    # Generate a breakpoint for interactive testing...
    #       raw_input("about to clean house...")
    
    # clean up the temp dir
    # In the OLD version of MMCE: There should only be files, no subdirs --Imran(?)
    # BUT, in mcce2.2:            There is an "energies" subdir -- added lines to remove these --vv
    if (cleanup):
        cwd = os.getcwd()
        os.chdir(tempdir)
        for i in os.listdir(tempdir):
            thisdir = os.path.join(tempdir,i) 
            if os.path.isdir(thisdir):
                for j in os.listdir(thisdir):
                    os.unlink(os.path.join(thisdir,j))
                os.rmdir(thisdir) 
            else:
                print 'removing', i, '...'
                os.unlink(thisdir)
        os.chdir(cwd)
        os.rmdir(tempdir)
        
    return pdbarr

def titrate(pdbfile, pHstart, pHstep, pHiters, mccepath, cleanup=True, prmfile=None, xtraprms={}):
    """Performs a pH titration on all titratable residues in a PDB file."""
    
    #Convert PDB file name to path, as paramgen expects an absolute path
    pdbpath=os.path.join(os.getcwd(),pdbfile)
    print 'pdbpath', pdbpath
    
    # Set up the parameters for the MCCE run
    params=paramgen(pdbpath, mccepath, fromfile=prmfile, xtraprms=xtraprms)
    
    # Calculate a titration curve with the following parameters
    params['TITR_PH0']=str(pHstart)
    params['TITR_PHD']=str(pHstep)
    params['TITR_STEPS']=str(pHiters)

    # Create a temporary directory with the run.prm file and run MCCE
    tempdir=tempfile.mkdtemp();
    print "Running MCCE in temporary directory %s..." % tempdir

    # remember our curent working directory
    thisdir = os.getcwd()
        
    # run mcce in the temporary directory
    os.chdir(tempdir)
    output = run_mcce(params)
    print output

    # chdir back to where we came from
    os.chdir(thisdir)
    
    # Read data from the temporary directory.
    pdbarr = ps_processMCCETitration(tempdir)
    
    # clean up the temp dir
    # In the OLD version of MMCE: There should only be files, no subdirs --Imran(?)
    # BUT, in mcce2.2:            There is an "energies" subdir -- added lines to remove these --vv
    if (cleanup):
        # Remove contents of temporary directory.
        for i in os.listdir(tempdir):
            thisdir = os.path.join(tempdir,i) 
            if os.path.isdir(thisdir):
                for j in os.listdir(thisdir):
                    os.unlink(os.path.join(thisdir,j))
                    os.rmdir(thisdir) 
            else:
                print 'removing', i, '...'
                os.unlink(thisdir)
        # Delete the temporary directory.
        os.rmdir(tempdir)
        
    return pdbarr

def protonatePDB(pdbfile, outfile, pH, mccepath, cleanup=True, prmfile=None, xtraprms={}):
    """Determine the most likely protonation state for a given protein structure.

    REQUIRED ARGUMENTS
        pdbfile           The path to the PDB file used as input for MCCE.  Must be ABSOLUTE path.  
        outfile           Name of PDB file to write in most likely protonation state.
        pH                Solvent pH at which to determine protonation state.
        mccepath          Path to MCCE executable.
                             
    OPTIONAL ARGUMENTS
        cleanup           If cleanup=True, ccleanup the temporary working dir for the MCCE calculations
        prmfile           If specified, will read MCCE parameters from file instead of default 
        xtraprms          Any additional extra MCCE parameters to be defined (overriding those from file or default)
                          These are specified as a dictionary of strings, e.g.:   {'MONTE_NEQ':'300'}
    """

    # Store current working directory.
    thisdir = os.getcwd()
    
    # Determine most likely protonation state, storing result in 'pdbarr'.
    pdbarr = protonation_state(pdbfile,pH,mccepath,cleanup=cleanup, prmfile=prmfile, xtraprms=xtraprms)
    
    # Write PDB file name to absolute path
    outpath=os.path.join(thisdir,outfile)
    fout = open(outpath,'w')
    for line in pdbarr:
        fout.write(line+'\n')
    fout.close()       

    # return from the directory we came from
    os.chdir(thisdir)
    
    return
        
def titratePDB(pdbfile, outfile, pHstart, pHstep, pHiters, mccepath, cleanup=True, prmfile=None, xtraprms={}):
    """Perform a pH titration on all titratable residues in a PDB file to determine their pKa values.

    REQUIRED ARGUMENTS
        pdbfile           The path to the PDB file used as input for MCCE.  Must be ABSOLUTE path  
        outfile           File to write the resilts of the pH titration calculation
        pHstart           pH value to start the titration
        pHstep            \Delta(pH) step size to increment the pH during titration
        pHiters           Number of total pH values to calculate,  including the startng pH and ending pH
                             
    OPTIONAL ARGUMENTS
        cleanup           If cleanup=True, ccleanup the temporary working dir for the MCCE calculations
        prmfile           If specified, will read MCCE parameters from file instead of default 
        xtraprms          Any additional extra MCCE parameters to be defined (overriding those from file or default)
                          These are specified as a dictionary of strings, e.g.:   {'MONTE_NEQ':'300'}
    """

    print 'in titratePDB: xtraprms', xtraprms
    thisdir = os.getcwd()
    pdbarr = titrate(pdbfile, pHstart, pHstep, pHiters, mccepath, cleanup=cleanup, prmfile=prmfile, xtraprms=xtraprms)

    # Write pK.out file name to absolute path
    outpath=os.path.join(thisdir,outfile)
    fout = open(outpath,'w')
    for line in pdbarr:
        fout.write(line.strip()+'\n')
    fout.close()       
        
    return

def ps_processmcce(tempdir):
    """Handles the file processing work for protonation_state
    
    ARGUMENTS
        tempdir - the temporary directory to run MCCE in.
    """

    # Build and use a regex to grab the appropriate entries from the MCCE PDB
    sstr,neut,pos,neg=ps_mostlikely(tempdir)
    rx=re.compile(sstr)
    print "searchstring: \"",sstr
    print "neut: \"",neut
    print "pos: \"",pos
    print "neg: \"",neg
    f=open(tempdir+"/step2_out.pdb","rt")
    pdbarr=filter(lambda x:rx.search(x),f.readlines())
    f.close()
        
    # Label the lines of the PDB file with charge state
    rxn=re.compile(neut)
    rxp=re.compile(pos)
    rxm=re.compile(neg)
    for i in range(len(pdbarr)):
        if (rxn.search(pdbarr[i])):
            pdbarr[i]=pdbarr[i]+"0"
        elif (rxp.search(pdbarr[i])):
            pdbarr[i]=pdbarr[i]+"+"
        elif (rxm.search(pdbarr[i])):
            pdbarr[i]=pdbarr[i]+"-"
        else:
            print "ERROR!"
            raise "Missed a case in charge state labeling"
        
    renumber_atoms(pdbarr)
    pdbarr=rename_residues(pdbarr)
    pdbarr=pdb_cleanup(pdbarr)

    return pdbarr

def ps_processMCCETitration(tempdir):
    """For now, handles the file processing work for titrate()
    
    ARGUMENTS
        tempdir - the temporary directory in which MCCE had been run
        
    RETURNS
        lines - concatenated output of pK.out and fort.38
    """
    
    fin = open(os.path.join(tempdir,'pK.out'),'r')
    lines = fin.readlines()
    fin.close()
    
    fin = open(os.path.join(tempdir,'fort.38'),'r')
    lines += fin.readlines()
    fin.close()
    
    return lines
            
            
