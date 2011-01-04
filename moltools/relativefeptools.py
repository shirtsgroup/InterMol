from openeye.oechem import *
from mmtools.moltools.ligandtools import *
from mmtools.gromacstools.MdpFile import *
import commands
import os.path
import shutil
from numpy import *


"""Tools for setting up relative binding free energy calculations with AMBER and gromacs.

REQUIREMENTS

  AMBER and Antechamber installations (in PATH).  Be sure to download the latest version of Antechamber separately.

TODO

AUTHORS

  Written by John D. Chodera <jchodera@stanford.edu>, Stanford, 10 Feb 2008.

  Some edits and generalization by David L. Mobley, dmobley@gmail.com, University of New Orleans, 15 November 2010.
"""

#=============================================================================================
# METHODS FOR READING, EXTRACTING, OR CREATING MOLECULES
#=============================================================================================

def write_file(filename, contents):
    """Write the specified contents to a file.

    ARGUMENTS
      filename (string) - the file to be written
      contents (string) - the contents of the file to be written

    """

    outfile = open(filename, 'w')

    if type(contents) == list:
        for line in contents:
            outfile.write(line)
    elif type(contents) == str:
        outfile.write(contents)
    else:
        raise "Type for 'contents' not supported: " + repr(type(contents))

    outfile.close()

    return

def read_file(filename):
    """Read contents of the specified file.

    ARGUMENTS
      filename (string) - the name of the file to be read

    RETURNS
      lines (list of strings) - the contents of the file, split by line
    """

    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    return lines

def select_off_section(lines, label):
   """Return set of lines associated with a section.

   ARGUMENTS
      lines (list of strings) - the lines in the .off file
      label (string) - the section label to locate

   RETURNS
      section_lines (list of strings) - the lines in that section
   """

   section_lines = list()

   nlines = len(lines)
   for start_index in range(nlines):
       # get line
       line = lines[start_index].strip()
       # split into elements
       elements = line.split()
       # see if keyword is matched
       if (len(elements) > 0):
           if (elements[0] == '!' + label):
               # increment counter to start of section data and abort search
               start_index += 1
               break
           
   # throw an exception if section not found
   if (start_index == nlines):
       raise "Section %(label)s not found." % vars()

   # Locate end of section.
   for end_index in range(start_index, nlines):
       # get line
       line = lines[end_index].strip()
       # split into elements
       elements = line.split()
       # see if keyword is matched
       if (len(elements) > 0):
           if (elements[0][0]=='!'):
               break

   # return these lines
   return lines[start_index:end_index]

def loadGAFFMolecule(ligand_basepath, ligand_name, ligand_mol2_suffix = '.openeye', debug = False):
    """Load molecule from OpenEye mol2 file and AMBER LEaP .off file with GAFF atom typenames and integer bond types.

    ARGUMENTS
      ligand_basepath (string) - name of directory containing both %(ligand_name)s.openeye.mol2 and %(ligand_name)s.off files.
      ligand_name (string) - prefix for filenames

    RETURNS
      ligand (OEMol) - the OEMol molecule with atom type names and integer bond types replaced by GAFF values read from LEaP .off file

    OPTIONAL:
        ligand_mol2_suffix: mol2 files may have an additional string before their extension, such as phenol.openeye.mol2 has the string '.openeye'. Default: '.openeye'. 


    NOTES
      Uses LEaP .off file to extract GAFF atom names, where entry for ligand is named 'ligand'.

    TODO
      Create molecule from scratch from only .off file, using NewAtom()?
      
    """

    # construct filenames
    ligand_mol2_filename = os.path.join(ligand_basepath, '%(ligand_name)s%(ligand_mol2_suffix)s.mol2' % vars())
    ligand_off_filename = os.path.join(ligand_basepath, '%(ligand_name)s.off' % vars())    

    if debug: print "Attempting to read files %s and %s..." % (ligand_mol2_filename, ligand_off_filename)

    if not os.path.exists(ligand_mol2_filename) or not os.path.exists(ligand_off_filename):
        print "Ligand %(ligand_name)s incomplete." % vars()
        return None  #NOTE: This should probably raise an exception instead, or at least optionally raise an exception -- clearly failing and just proceeding is generally bad form. 

    # load file
    ligand = readMolecule(ligand_mol2_filename, normalize = True)

    # override title
    ligand.SetTitle(ligand_name)

    # read .off file
    off_lines = read_file(ligand_off_filename)

    # build a dictionary of atoms
    atoms = dict()
    for atom in ligand.GetAtoms():
        name = atom.GetName()
        atoms[name] = atom

    # replace atom names with GAFF atom names
    section_lines = select_off_section(off_lines, 'entry.ligand.unit.atoms')
    off_atom_names = dict() # off_atom_names[i] is name of atom i, starting from 1
    atom_index = 1
    for line in section_lines:
        # parse section
        # !entry.ligand.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg
        # "C1" "c3" 0 1 131072 1 6 -0.073210 
        elements = line.split()
        name = elements[0].strip('"')
        type = elements[1].strip('"')
        charge = float(elements[7])
        # store changes
        atom = atoms[name]
        atom.SetType(type)
        atom.SetPartialCharge(charge)
        # store atom ordering
        off_atom_names[atom_index] = name
        atom_index += 1
    
    # build a dictionary of bonds
    bonds = dict()
    for bond in ligand.GetBonds():
        begin_atom_name = bond.GetBgn().GetName() # name of atom at one end of bond
        end_atom_name = bond.GetEnd().GetName() # name of atom at other end of the bond
        # store bond in both directions
        bonds[(begin_atom_name,end_atom_name)] = bond
        bonds[(end_atom_name,begin_atom_name)] = bond

    # replace bond types with GAFF integral bond types
    section_lines = select_off_section(off_lines, 'entry.ligand.unit.connectivity')
    for line in section_lines:
        # parse section
        elements = line.split()
        i = int(elements[0])
        j = int(elements[1])
        bondtype = int(elements[2])
        # get atom names
        begin_atom_name = off_atom_names[i]
        end_atom_name = off_atom_names[j]
        # store changes
        bond = bonds[(begin_atom_name,end_atom_name)]
        bond.SetIntType(bondtype)

    # DEBUG
#    print "atoms"
#    for atom in ligand.GetAtoms():
#        print '%4s %4s' % (atom.GetName(), atom.GetType())
#    print "\nbonds"
#    for bond in ligand.GetBonds():
#        print "%4s %4s %d" % (bond.GetBgn().GetName(), bond.GetEnd().GetName(), bond.GetIntType())
#    print ""

    return ligand

def totalPartialCharge(ligand):
    """Compute the total partial charge of all atoms in molecule.

    ARGUMENTS
      ligand (OEMol) - the molecule

    RETURNS
      total_charge (float) - the total charge
    """

    total_charge = 0.0
    for atom in ligand.GetAtoms():
        total_charge += atom.GetPartialCharge()

    return total_charge

def totalIntegralPartialCharge(ligand):
    """Compute the total partial charge of all atoms in molecule and return nearest integer.

    ARGUMENTS
      ligand (OEMol) - the molecule

    RETURNS
      total_charge (int) - the total charge rounded to the nearest integer
    """

    total_charge = 0.0
    for atom in ligand.GetAtoms():
        total_charge += atom.GetPartialCharge()

    total_charge = int(round(total_charge))

    return total_charge

def adjustCharges(q_n, target_net_charge):
    """Adjust charges q_n to have desired net integral charge but closely match original charges.

    New charged qnew_n are defined by

      qnew_n = alpha (q_n - mu) + b

    where
      mu = Q/N
      Q = \sum_n q_n

    and alpha and b are chosen to minimize RMS(qnew_n - q_n) while guaranteeing \sum_n qnew_n = net_charge.      

    ARGUMENTS
      q_n (numpy float array of length N) - original charges
        CHANGELOG:
            11/15/10: DLM, switched exception handling from raising a string to using a class, as string exceptions are no longer supported in newer Python versions.target_net_charge (integer or float) - the desired net charge

    RETURNS
      qnew_n (numpy float array of length N) - adjusted charges
    """

    # get number of charges
    N = q_n.size

    # compute original total charge
    Q = sum(q_n)

    # compute average charge
    mu = Q / float(N)

    # determine scaling constant alpha and shift value b to minimize RMS and guarantee new net charge
    alpha = sum((q_n - mu)**2 * (1.0 - mu/q_n)) / sum((q_n-mu)**2)
    b = 0

    # compute new charges
    qnew_n = alpha * (q_n - mu) + b

    # check constraints
    Qnew = sum(qnew_n)
    if (abs(Qnew - target_net_charge) > 0.01):
        print "%6s %6s" % ('old', 'new')
        for n in range(N):
            print "%6.3f %6.3f" % (q_n[n], qnew_n[n])                                  
        raise "adjustCharges failed: original net charge %f, new net charge is %f" % (Q, Qnew)
    
    # compute RMS
    rms = sqrt( (1.0/float(N)) * sum((qnew_n - q_n)**2) )
    print "rms fit = %f" % rms

    # return new charges
    return qnew_n

def determineCommonSubstructure(ligands, debug = False, min_atoms = 4):
    """Find a common substructure shared by all ligands.

    The atom type name strings and integer bond types are used to obtain an exact match.

    ARGUMENTS
      ligands (list of OEMol) - the set of ligands for which the common substructure is to be determined.

    OPTIONAL ARGUMENTS
      debug (boolean) - if True, debug information is printed (default: False)
      min_atoms (int) - minimum number of atoms for substructure match (default: 4)

    RETURN VALUES
      common_substructure (OEMol) - a molecule fragment representing the common substructure
    
    CHANGELOG
      11/15/10 - DLM: Removing the attempt to convert OEMols to OEMol with OEMol( as this raises an exception in my (newer) OEChem. 
    """
    # Determine number of ligands
    nligands = len(ligands)

    # First, initialize with first ligand.
    common_substructure = ligands[0].CreateCopy() #DLM modification 11/15/10 -- this is how copies should now be made

    # DEBUG: Show original atoms
    if debug:
        atom_index = 1
        for atom in common_substructure.GetAtoms():
            print "%5d %6s %12s" % (atom_index, atom.GetName(), atom.GetType())
            atom_index += 1
        print ""

    # Now delete bits that don't match every other ligand.
    for ligand in ligands[1:]:
        # get ligand name
        ligand_name = ligand.GetTitle()
        
        # Create an OEMCSSearch from this molecule.
        mcss = OEMCSSearch(ligand, OEExprOpts_StringType, OEExprOpts_IntType)

        # ignore substructures smaller than 4 atoms
        mcss.SetMinAtoms(min_atoms)

        # perform match
        for match in mcss.Match(common_substructure):
            nmatched = match.NumAtoms()
            
            if debug: print "%(ligand_name)s : match size %(nmatched)d atoms" % vars()
            
            # build list of matched atoms in common substructure
            matched_atoms = list()
            for matchpair in match.GetAtoms():
                atom = matchpair.target
                matched_atoms.append(atom)

            # delete all unmatched atoms from common substructure
            for atom in common_substructure.GetAtoms():
                if atom not in matched_atoms:
                    common_substructure.DeleteAtom(atom)

            if debug:
                # print all retained atoms
                atom_index = 1
                for atom in common_substructure.GetAtoms():
                    print "%5d %6s %12s" % (atom_index, atom.GetName(), atom.GetType())
                    atom_index += 1
                print ""

            # we only need to consider one match
            break
    
    # return the common substructure
    return common_substructure

def determineMinimumRMSCharges(common_substructure, ligands, debug = False, min_atoms = 4):
    """Determine charges for common substructure that minimize root-mean-squared (RMS) deviation to all ligands while having same net charge as ligands.

    ARGUMENTS
      common_substructure (OEMol) - the common substructure (partial molecule) that matches the ligands
      ligands (list of OEMol) - the set of ligands for which the common_substructure is a shared component

    OPTIONAL ARGUMENTS
      debug (boolean) - if True, debug information is printed (default: False)
      min_atoms (int) - minimum number of atoms for MCSS match (default: 4)
      
    RETURNS
      charged_substructure (OEMol) - the common substructure with atomic charges updated

    NOTES
      An exception will be thrown if the ligands do not have the same total integral charge.
    
    REFERENCES
      [1] Betts JT. A compact algorithm for computing the stationary point of a quadratic function subject to linear constraints. ACM Trans Math Soft 6(3):391-397, 1980.
      
    CHANGELOG:
        11/15/10: DLM, switched exception handling from raising a string to using a class, as string exceptions are no longer supported in newer Python versions.

"""

    if debug: print "Determining minimum RMS charges..."

    # maximum allowed total charge deviation before raising an exception
    CHARGE_DEVIATION_TOLERANCE = 0.01

    # determine total integral charge on first ligand
    Q = totalIntegralPartialCharge(ligands[0])
    if debug: print "Q = %f" % Q


    class ChargeError(Exception):
        def __init__(self, msg, expr):
            self.expr = expr
            self.msg = msg
            print msg
            print "Exiting..."
        

    # check to make sure all ligands have same net charge
    for ligand in ligands:
        if abs(totalPartialCharge(ligand) - Q) > CHARGE_DEVIATION_TOLERANCE:
            #raise "Ligand %s has charge (%f) different from target charge (Q = %f) - tolerance is %f." % (ligand.GetTitle(), totalPartialCharge(ligand), Q, CHARGE_DEVIATION_TOLERANCE)
            raise ChargeError( "Ligand %s has charge (%f) different from target charge (Q = %f) - tolerance is %f." % (ligand.GetTitle(), totalPartialCharge(ligand), Q, CHARGE_DEVIATION_TOLERANCE), 'abs(totalPartialCharge(ligand) - Q) > CHARGE_DEVIATION_TOLERANCE' )

    # determine number of ligands
    K = len(ligands)
    if debug: print "K = %d" % K
          
    # determine number of atoms in common substructure
    N = common_substructure.NumAtoms()
    if debug: print "N = %d" % N

    # build index of which atoms correspond to which atom indices
    atom_indices = dict() # atom_indices[atom] is the index (0...(N-1)) of atom 'atom'
    n = 0
    for atom in common_substructure.GetAtoms():
        atomname = atom.GetName()
        atom_indices[atomname] = n
        n += 1

    # extract charge on all ligands for common set
    q_kn = zeros([K,N], float64) # q_kn[k,n] is the charge on the atom of ligand k corresponding to atom n of the common substructure
    k = 0 # ligand index
    mcss = OEMCSSearch(common_substructure, OEExprOpts_StringType, OEExprOpts_IntType)        
    mcss.SetMinAtoms(min_atoms)
    for ligand in ligands:
        print ligand.GetTitle()
        # perform match
        for match in mcss.Match(ligand):
            # only store the first match
            for matchpair in match.GetAtoms():
                atomname = matchpair.pattern.GetName()
                n = atom_indices[atomname]
                q_kn[k,n] += matchpair.target.GetPartialCharge()                
            break
        # increment ligand counter
        k += 1
    
    # compute average charge on common substructure over all ligands
    qavg_n = mean(q_kn, 0)
    print qavg_n.shape

    # solve linear programming problem to determine minimum RMS charge set subject to constraint that total charge equals Q
    # 
    # See reference [1] for information on how the matrix is set up.
    #
    # We solve the linear programming problem A x = b where
    # 
    # A = [I 1; 1' 0]
    # b = [qavg_n ; Q]
    # x = [q_n ; lambda]

    A = zeros([N+1,N+1], float64)
    A[0:N,0:N] = eye(N)
    A[0:N,N] = ones([N])
    A[N,0:N] = ones([N])    

    b = zeros([N+1,1], float64)
    b[0:N,0] = qavg_n[0:N]
    b[N,0] = Q

    import numpy.linalg
    Ainv = numpy.linalg.pinv(A)
    x = matrix(Ainv) * matrix(b)

    # extract charges
    q_n = array(x[0:N])
    lagrange_multiplier = x[N]
    if debug: print "lagrange multiplier is %f" % lagrange_multiplier

    # DISPLAY CHARGES
    if debug:
        # header
        print "%6s" % 'fit',
        for k in range(K):
            print " %6d" % k, 
        print ""
        # charges
        for n in range(N):
            print "%6.3f" % q_n[n],
            for k in range(K):
                print " %6.3f" % q_kn[k,n],
            print ""
        print ""
        print ""
        # print total charges
        print "%6.3f" % sum(q_n[:]),
        for k in range(K):
            print " %6.3f" % sum(q_kn[k,:]),
        print ""        

    # make a copy of the common substructure to assign new charges to
    charged_common_substructure = OEMol(common_substructure)

    # assign charges
    n = 0
    for atom in charged_common_substructure.GetAtoms():
        charge = float(q_n[n])
        atom.SetPartialCharge(charge)
        n += 1

    # return copy of common substructure with partial charges assigned
    return charged_common_substructure


def pairwise_common_substructures( ligand_basepath, ligand_list, Debug = False ):
    """Functionalized version of pairwise-common-substructures script written by J. Chodera; this version written by D. Mobley. Takes a list of ligands file names and returns (PRESENTLY NOTHING -- TEXT OUTPUT ONLY) 

    ARGUMENTS:
        ligand_basepath: Base directory where these ligands are found; will have names as given in ligand_list, and extensions are assumed to be .off, Amber .off format
        ligand_list: List of ligand name prefixes (extensions assumed to be .off)

    OPTIONAL ARGUMENTS:
        Debug: If set to true, print extra debugging information. Default: False

    """

    nligands = len(ligand_list)

    for ligand1_index in range(len(ligand_list)):
        # get ligand name
        ligand1_name = ligand_list[ligand1_index]

        # Load molecule with GAFF atom names and integer bondtypes.
        ligand1 = loadGAFFMolecule(ligand_basepath, ligand1_name)
        if not ligand1: 
            print "Error reading %s..." % ligand1_name
            continue

        # Create an OEMCSSearch from this molecule.
        mcss = OEMCSSearch(ligand1, OEExprOpts_StringType, OEExprOpts_IntType)
        # ignore substructures smaller than 4 atoms
        mcss.SetMinAtoms(4)
            
        # Consider all other ligands.
        for ligand2_index in range(len(ligand_list)):
            # skip self
            if ligand1_index == ligand2_index: continue

            # get name
            ligand2_name = ligand_list[ligand2_index]

            # Load the molecule
            ligand2 = loadGAFFMolecule(ligand_basepath, ligand2_name)
            if not ligand2: 
                print "Error reading %s..." % ligand2_name
                continue


            # Compare atomsets
            matchcount = 0
            maxmatch = 0
            for match in mcss.Match(ligand2):
                if Debug: #Print lots of extra debugging information if desired
                    print "Match:", matchcount, "Size:", match.NumAtoms(), " atoms"
                    # determine mutated charges and atomtypes
                    (mutated_charges_1, mutated_charges_2) = determineMutatedCharges(ligand1, ligand2, match)
                    (mutated_atomtypes_1, mutated_atomtypes_2) = determineMutatedAtomtypes(ligand1, ligand2, match)

                    # DEBUG
                    print "ligand1: %s" % ligand1.GetTitle()
                    for atom in ligand1.GetAtoms():
                        atomname = atom.GetName()
                        print "%6s : %16s %6.3f -> %16s %6.3f" % (atom.GetName(), atom.GetType(), atom.GetPartialCharge(), mutated_atomtypes_1[atomname], mutated_charges_1[atomname])
                    print "ligand2: %s" % ligand2.GetTitle()
                    for atom in ligand2.GetAtoms():
                        atomname = atom.GetName()
                        print "%6s : %16s %6.3f -> %16s %6.3f" % (atom.GetName(), atom.GetType(), atom.GetPartialCharge(), mutated_atomtypes_2[atomname], mutated_charges_2[atomname])
                    print "%4s %6s        %4s %6s   %6s" % ('atom', 'charge', 'atom', 'charge', 'charge')
                    index = 0
                    for mp in match.GetAtoms():
                        atom1 = mp.pattern
                        atom2 = mp.target
                        print "%4s %6.3f <----> %4s %6.3f : %6.3f" % (atom1.GetName(), atom1.GetPartialCharge(), atom2.GetName(), atom2.GetPartialCharge(), common_charge[index])
                        index +=1
                    print
                matchcount +=1
                maxmatch = max(maxmatch, match.NumAtoms() )

            print '%(ligand1_name)12s : %(ligand2_name)12s    %(maxmatch)8d atoms' % vars()

    return #Presently returns nothing
