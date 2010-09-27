#=============================================================================================
# Parameterize ligands.
# 
# Written by John D. Chodera <jchodera@gmail.com> 2008-02-01
#=============================================================================================

#=============================================================================================
# Imports
#=============================================================================================

from mmtools.moltools.ligandtools import *
from mmtools.gromacstools.MdpFile import *
import commands
import os
import os.path
import shutil

#=============================================================================================

#=============================================================================================
# PARAMETERS
#=============================================================================================

ligand_basepath = './ligands/' # base path for all modelsets
work_basepath = './ligands-parameterized/' # base path to set up all calculations in

# list of ligands to consider
ligand_name_list = [ ('LG%d-0-eq' % index) for index in (2, 5, 6, 8, 9) ]
extension = "mol2" # ligand filename extension (could change this to "mol2" if desired)

#=============================================================================================
# SUBROUTINES
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
#=============================================================================================
def parameterize_ligand(ligand_filename, ligand_name, work_path):
    """Parameterize ligand.

    ARGUMENTS
      ligand_filename (string) - name of mol2 file describing ligand
      work_path (string) - name of directory to place files in

    """

    # get current directory
    current_path = os.getcwd()

    # create the work directory if it doesn't exist
    if not os.path.exists(work_path):
        os.makedirs(work_path)

    # SET UP LIGAND TOPOLOGY

    print "\nCONSTRUCTING LIGAND TOPOLOGY"

    # Read the ligand into OEMol.
    print "Reading ligand..."
    ligand = readMolecule(ligand_filename, normalize = True)

    # Set name
    print "ligand_name = %s" % ligand_name
    ligand.SetTitle(ligand_name)

    # Get formal charge of ligand.
    ligand_charge = formalCharge(ligand)
    print "formal charge is %d" % ligand_charge

    # Write molecule with explicit hydrogens to mol2 file.
    print "Writing ligand mol2 file..."
    ligand_mol2_filename = os.path.abspath(os.path.join(work_path, '%(ligand_name)s.openeye.mol2' % vars()))
    writeMolecule(ligand, ligand_mol2_filename)

    # Set substructure name (which will become residue name).
    print "Modifying molecule name..."
    modifySubstructureName(ligand_mol2_filename, 'MOL')

    # Run antechamber to assign GAFF atom types.
    print "Running antechamber..."
    os.chdir(work_path)
    gaff_mol2_filename = os.path.join(work_path, '%(ligand_name)s.gaff.mol2' % vars())
    charge_model = 'bcc'
    command = 'antechamber -i %(ligand_mol2_filename)s -fi mol2 -o %(gaff_mol2_filename)s -fo mol2 -c %(charge_model)s -nc %(ligand_charge)d >& antechamber.out' % vars()
    print command
    output = commands.getoutput(command)    
    os.chdir(current_path)

    # skip if antechamber didn't work
    if not os.path.exists(gaff_mol2_filename):
        return

    # Generate frcmod file for additional GAFF parameters.
    ligand_frcmod_filename = os.path.join(work_path, '%(ligand_name)s.frcmod' % vars())
    output = commands.getoutput('parmchk -i %(gaff_mol2_filename)s -f mol2 -o %(ligand_frcmod_filename)s' % vars())
    print output

    # skip if parmchk didn't work
    if not os.path.exists(ligand_frcmod_filename):
        return

    # Run LEaP to generate topology / coordinates.
    ligand_prmtop_filename = os.path.join(work_path,'%(ligand_name)s.prmtop' % vars())
    ligand_crd_filename = os.path.join(work_path,'%(ligand_name)s.crd' % vars())
    ligand_off_filename = os.path.join(work_path, '%(ligand_name)s.off' % vars())

    tleap_input_filename = os.path.join(work_path, 'setup-ligand.leap.in')
    tleap_output_filename = os.path.join(work_path, 'setup-ligand.leap.out')
    contents = """
# Load GAFF parameters.
source leaprc.gaff

# load antechamber-generated additional parameters
mods = loadAmberParams %(ligand_frcmod_filename)s

# load ligand
ligand = loadMol2 %(gaff_mol2_filename)s

# check the ligand
check ligand

# report net charge
charge ligand

# save AMBER parameters
saveAmberParm ligand %(ligand_prmtop_filename)s %(ligand_crd_filename)s

# write .off file
saveOff ligand %(ligand_off_filename)s

# exit
quit
""" % vars()
    write_file(tleap_input_filename, contents)
    command = 'tleap -f %(tleap_input_filename)s > %(tleap_output_filename)s' % vars()
    output = commands.getoutput(command)

    # extract total charge
    ligand_charge = commands.getoutput('grep "Total unperturbed charge" %(tleap_output_filename)s | cut -c 27-' % vars())
    ligand_charge = int(round(float(ligand_charge))) # round to nearest whole charge
    print "ligand charge is %d" % ligand_charge    

    return

#=============================================================================================
# MAIN
#=============================================================================================

# clean out leap.log
commands.getoutput('rm -f leap.log')

for ligand in ligand_name_list:
    # construct pathnames
    ligand_filename = os.path.join(ligand_basepath, '%(ligand)s.%(extension)s' % vars())

    if not os.path.exists(ligand_filename):
        print "%s not found" % ligand_filename
        continue

    # parameterize ligand
    parameterize_ligand(ligand_filename, ligand, work_basepath)

