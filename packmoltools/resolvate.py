#!/usr/bin/python

#=============================================================================================
# resolvate
#
# Resolvate a PDB file and produce an AMBER .crd file.
#=============================================================================================
# AUTHORS
#
# Written by John Chodera, Stanford University, 2007-03-05
#=============================================================================================
# REQUIREMENTS
# - Packmol: http://www.ime.unicamp.br/~martinez/packmol/
# - numarray
#=============================================================================================
# TODO
#=============================================================================================
# VERSION CONTROL INFORMATION
__version__ = "$Revision: $"                                                                                              
#=============================================================================================
# MOD HISTORY
# 7/10/07 - Can call resolvate.py from cmd line or use as resolvate.resolvate_from_template()
#         - Can now read in *.gro templates for box volumes
#         - Can now write solvated structures to *.pdb format
#         - Added --tol option to specify the minimum distance (A) that the solvent can be packed 
#         - PDB reading and writing now peeled off into mmtools/pdbtools.py
#         - atom["resName"] is now 4-letter code
#=============================================================================================
# IMPORTS
import sys
import os.path
from optparse import OptionParser # For parsing of command line arguments
import commands
import shutil

import tempfile
import math


from mmtools.gromacstools.GromacsData import *

from mmtools.pdbtools import *

#=============================================================================================

#=============================================================================================
# PARAMETERS
#=============================================================================================

# here, we assume that packmol is in your path
packmol = 'packmol'
# packmol = '${HOME}/local/src/packmol/src/packmol'

#=============================================================================================
# SUBROUTINES
#=============================================================================================



def getPresentSequence(pdbfilename, chain=' '):
    """Extract the sequence for which there are atomic coordiantes defined from a PDB file.

    present_sequence = getPresentSequence(pdbfilename, chain=' ')
    
    REQUIRED ARGUMENTS
      pdbfilename - the filename of the PDB file to import from

    OPTIONAL ARGUMENTS
      chain - the one-character chain ID of the chain to import (default ' ')

    RETURN VALUES
      sequence - array of residue names

    The ATOM records are read, and the sequence for which there are atomic coordinates is stored.

    """

    # Read the PDB file into memory.
    pdbfile = open(pdbfilename, 'r')
    lines = pdbfile.readlines()
    pdbfile.close()

    # Extract the sequence for which there are defined atomic coordinates.
    sequence = [ ]
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
        if (chain == field['chainID']):
          if(field["resSeq"] != last_resSeq):
            sequence.append(field["resName"])
            last_resSeq = field["resSeq"]

    # Return dictionary of present residues.
    return sequence


def readAmberCrd(crd_filename):
  """Read an AMBER format .crd file.

  REQUIRED ARGUMENTS
    crd_filename - name of AMBER .crd file to read

  RETURNS
    title - title string
    natoms - number of atoms
    coordinates - 3*natoms array of coordinates
    box_dimensions - box coordinates if present, or None if not
    box_angles - box angles if present, or None if not    
  
  EXAMPLE
    crd_filename = 'output.crd'
    (title, natoms, coordinates, box_dimensions, box_angles) = readAmberCrd(crd_filename)

  """
  
  # Read contents of file.
  infile = open(crd_filename, 'r')
  lines = infile.readlines()
  infile.close()

  # Parse header.
  title = lines.pop(0).rstrip(' \n')
  natoms = int(lines.pop(0).split()[0])

  # Parse coordinates.
  coordinates = [ ]
  for line in lines:
    elements = line.split()
    for element in elements:
      coordinates.append(float(element))

  # Extract box dimensions and coordinates, if present.
  if (len(coordinates) >= 3*natoms + 3):
    box_dimensions = coordinates[3*natoms:3*natoms+3]
  if (len(coordinates) == 3*natoms + 6):
    box_angles = coordinates[3*natoms+3:3*natoms+6]
  coordinates = coordinates[0:3*natoms]

  # Return results
  return (title, natoms, coordinates, box_dimensions, box_angles)  


def readGromacsGro(gro_filename):
  """Read a GROMACS format .gro file.

  REQUIRED ARGUMENTS
    gro_filename - name of GROMACS .gro file to read

  RETURNS
    title - title string
    natoms - number of atoms
    coordinates - 3*natoms array of coordinates
    box_dimensions - box coordinates if present, or None if not
    box_angles - box angles if present, or None if not    
  
  EXAMPLE
    crd_filename = 'output.crd'
    (title, natoms, coordinates, box_dimensions, box_angles) = readGromacsGro(gro_filename)

  """
  
  # create a GromacsStructure object (yet to be filled with data)
  g = GromacsStructure(name=gro_filename)

  # load the *.gro file data into the gromacs structure
  g.load(gro_filename)
  
  # Parse header.
  title = g.header
  natoms = g.natoms

  # Parse coordinates, and find the maximum dimensions of the *.gro coords' "brick"
  coordinates = [ ]
  min_xyz = [None,None,None]
  max_xyz = [None,None,None]
  for atom in g.atoms:
    elements = atom.position()
    for i in range(0,len(elements)):
      coordinates.append(float(elements[i])*10.)    # gro files are in nm, so we must multiply by 10!
      if (min_xyz[i] == None) or (min_xyz[i] > elements[i]*10.):
          min_xyz[i] = elements[i]*10.
      if (max_xyz[i] == None) or (max_xyz[i] < elements[i]*10.):
          max_xyz[i] = elements[i]*10.

  # Parse the boxstring to extract box dimensions and coordinates, if present.
  
  # NOTES on *.gro file tails
  # The last line of a *.gro file contain the box vectors (free format, space separated reals), values:
  #     v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
  # the last 6 values may be omitted (they will be set to zero). Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0.
  # The box_angles we need to calculate are [ angle(v2,v3), angle(v1,v3), angle(v1,v2) ]
  
  if (0):
    boxfields = g.boxstring.split()
    if len(boxfields) >= 3:
      box_dimensions = [ float(boxfields[0]), float(boxfields[1]), float(boxfields[2]) ]
      if len(boxfields) >= 6:
          v1 = [ float(boxfields[0]), float(boxfields[3]), float(boxfields[4]) ]
          v2 = [ float(boxfields[1]), float(boxfields[5]), float(boxfields[6]) ]
          v3 = [ float(boxfields[2]), float(boxfields[7]), float(boxfields[8]) ]
          angle23 = 180./math.pi * math.acos( dot(v2,v3)/( math.sqrt( dot(v2,v2) * dot(v3,v3) ) ) )           
          angle13 = 180./math.pi * math.acos( dot(v1,v3)/( math.sqrt( dot(v1,v1) * dot(v3,v3) ) ) )           
          angle12 = 180./math.pi * math.acos( dot(v1,v2)/( math.sqrt( dot(v1,v1) * dot(v2,v2) ) ) )           
      else:
          box_angles = None
    else:
      box_dimensions = None
      box_angles = None

  # MORE NOTES:   We are going to BYPASS this geometry, and just work with the brick coordinates.  GROMACS
  # writes the coords as a rectangular "brick" regardless of the unit cell, for numerical efficiency
  
  if (1):
    box_dimensions = [ (max_xyz[0]-min_xyz[0]), (max_xyz[1]-min_xyz[1]), (max_xyz[2]-min_xyz[2]) ]
    box_angles = None

  # Return results
  return (title, natoms, coordinates, box_dimensions, box_angles)  


def dot(v1, v2):
    """Returns dot product of two vectors v1, v2 that are a list of floats."""
    
    if len(v1) != len(v2):
        print 'Error: vectors must be the same size'
        raise ValueError    
    result = 0.
    for i in range(0,len(v1)):
        result = result + v1[i]*v2[i]
    return result


def writeAmberCrd(crd_filename, coordinates, title = '', box_dimensions = None, box_angles = None):
  """Write an AMBER format .crd file.
  
  REQUIRED ARGUMENTS
    crd_filename - the filename of the AMBER .crd file to be written
    coordinates - array or list of length 3*Natoms of atomic coordinates
    box_coordinates - box coordinates to be appended to end

  RETURNS
    none

  EXAMPLE
    crd_filename = 'output.crd'
    coordinates = [1.0, 2.0, 3.0]
    writeAmberCrd(crd_filename, natoms, coordinates)  
    
  """
  # Check to make sure number of coordinates is divisiable by 3.
  if(len(coordinates) % 3 != 0):
    raise ParameterError('Number of atomic cordinates is not a multiple of 3.')
  
  # Determine number of atoms.
  natoms = len(coordinates) / 3
  
  # Append box_coordinates.
  coordinates_to_write = coordinates
  if (box_dimensions):
    coordinates_to_write += box_dimensions
  if (box_angles):
    coordinates_to_write += box_angles
      
  # Open file to write.
  outfile = open(crd_filename, 'w')
  
  # Write header.
  outfile.write('%s\n' % title)
  outfile.write('%6d\n' % natoms)
  
  # Write coordinates.
  coordinates_this_line = 0
  for coordinate in coordinates_to_write:
    # Write a coordinate.
    outfile.write('%12.7f' % coordinate)
    coordinates_this_line += 1
    
    # Wrap if we have written 6 coordinates/line.
    if (coordinates_this_line == 6):
      outfile.write('\n')
      coordinates_this_line = 0

  outfile.close()
  return


def isgrofile(grofile):
    """Returns True if the file ends in *.gro."""
    if len(grofile) <= 4:
        return False    
    if grofile[-4:] == '.gro':
        return True
    else:
        return False
    
def iscrdfile(crdfile):
    """Returns True if the file ends in *.crd."""
    if len(crdfile) <= 4:
        return False    
    if crdfile[-4:] == '.crd':
        return True
    else:
        return False
    
def ispdbfile(pdbfile):
    """Returns True if the file ends in *.crd."""
    if len(pdbfile) <= 4:
        return False    
    if pdbfile[-4:] == '.pdb':
        return True
    else:
        return False
    


def resolvate_from_template(reference_pdb, reference_crd, source_pdb, output_crd, tolerance=2.0):
  """Resolvates a PDB file based on an AMBER *.gro or a GROMACS *.gro template (which provides the box parameters)."""

  packmol = 'packmol'
  
  # Create temporary directory.
  tmpdir = tempfile.mkdtemp()
  print "tmpdir is %s" % tmpdir
  
  # Extract sequence from reference PDB file.
  reference_sequence = getPresentSequence(reference_pdb)
  
  # Extract seqence from source PDB file.
  source_sequence = getPresentSequence(source_pdb)
  
  # Determine number of water molecules as difference in sequence lengths.
  nwat = len(reference_sequence) - len(source_sequence)
  print "Will add %d water molecules." % nwat
  
  # Determine file type of the reference_crd file
  refcrdFileType = None
  if iscrdfile(reference_crd):
     refcrdFileType = 'crd'
  if isgrofile(reference_crd):
     refcrdFileType = 'gro'
  if refcrdFileType == None:
      print 'reference_crd =',reference_crd,'must either be *.crd or *.gro file!'
      raise IOError
      
  # Extract box coordinates
  if refcrdFileType=='crd':
      (title, natoms, coordinates, box_dimensions, box_angles) = readAmberCrd(reference_crd)
  else:
      (title, natoms, coordinates, box_dimensions, box_angles) = readGromacsGro(reference_crd)

  
  # Compute box parameters.
  print 'box_dimensions', box_dimensions
  (box_x, box_y, box_z) = box_dimensions
  # create a slightly-smaller box to prevent too-close waters on the pbc borders
  (box_x, box_y, box_z) = (box_x-tolerance, box_y-tolerance, box_z-tolerance)
  (center_x, center_y, center_z) = (box_x/2.0, box_y/2.0, box_z/2.0)
  
  # Get date.
  date = commands.getoutput('date')
  
  # Temporary output PDB filename.
  output_pdb = os.path.join(tmpdir, 'output.pdb')
  
  # Generate water reference file from atoms in first non-solute residue in reference PDB file.
  atoms = readAtomsFromPDB(reference_pdb)
  reference_water_residue_atoms = []
  first_water_residue = len(source_sequence) + 1
  if refcrdFileType=='crd':
    for atom in atoms:
      if (atom['resSeq'] == first_water_residue):
        reference_water_residue_atoms.append(atom)
  # for the gromacs files, the first non-solute residue is an ion, so we must skip these
  else:    
    # find the first residue with resName 'HOH' or 'SOL'       
    first_water_residue = None
    atom_index = 0
    while first_water_residue==None:
      atom = atoms[atom_index]
      if (atom['resName'].strip() == 'HOH') or (atom['resName'].strip() == 'SOL'):
          first_water_residue = atom['resSeq']
      atom_index += 1
    # chosose this water residue for wat.pdb
    for atom in atoms:
      if (atom['resSeq'] == first_water_residue):
        reference_water_residue_atoms.append(atom)
      
  reference_water_pdb = os.path.join(tmpdir, 'wat.pdb')
  writeAtomsToPDB(reference_water_pdb, reference_water_residue_atoms, renumber=True)
  
  # Copy source PDB file.
  shutil.copy(source_pdb, os.path.join(tmpdir, 'source.pdb'))
  
  # Generate title.
  title = "Resolvation of %(source_pdb)s in %(box_x)f A x %(box_y)f A x %(box_z)f A box by %(nwat)d water molecules" % vars()
  
  # Construct Packmol input file.
  packmol_input_filename = os.path.join(tmpdir, 'packmol.in')
  packmol_input = """\
  #
  # Generated by resolvate.py on %(date)s
  # %(title)s
  #
  
  # All atoms from diferent molecules will be at least 'tolerance' A apart at the solution
  
  tolerance %(tolerance)f
  
  # The type of the files will be pdb 
  
  filetype pdb
  
  # The name of the output file
  
  output output.pdb
  
  # The protein will be fixed with its center of mass at center of the
  # box, and no rotation (the first three zeros correspond to the position
  # of the center of mass and the last three correspond to the euler
  # angles of rotation, in radian, relative to the position in the input
  # file). 
  
  structure source.pdb
    number 1 
    resnumbers 1
    fixed %(center_x)f %(center_y)f %(center_z)f 0. 0. 0.
    centerofmass
  end structure
  
  # Water molecules will be put inside a box that contains the protein and ions.
  structure wat.pdb
    number %(nwat)d
    resnumbers 0
    inside box 0.0 0.0 0.0 %(box_x)f %(box_y)f %(box_z)f
  end structure
  
  """ % vars()
  
  outfile = open(packmol_input_filename, 'w')
  outfile.write(packmol_input)
  outfile.close()
  
  # Execute Packmol.
  print "Running packmol..."
  command = 'cd %(tmpdir)s ; %(packmol)s < %(packmol_input_filename)s' % vars()
  output = commands.getoutput(command)
  print output

  # Read the output PDB file.
  atoms = readAtomsFromPDB(os.path.join(tmpdir, 'output.pdb'))
  
  # Clean up temporary directory.
  for filename in os.listdir(tmpdir):
    os.remove(os.path.join(tmpdir,filename))
  os.rmdir(tmpdir)

  # write the solvated coordinate file to output_crd
  if ispdbfile(output_crd):
    # find which atoms are solute vs. solvent
    solvent_atoms = []
    solute_atoms = []
    for atom in atoms:
      # Solvent.
      if (atom['chainID'] == 'A'):
        solvent_atoms.append(atom)
      # Solute.
      if (atom['chainID'] == 'B'):
        solute_atoms.append(atom)
    # rearrange the order
    atoms = solute_atoms + solvent_atoms
    writeAtomsToPDB(output_crd, atoms, renumber = True)
  else:    # default output is as AMBER crd file
    # Get solvent coordinates (chain A) and solute coordinates (chain B)
    solvent_coordinates = []
    solute_coordinates = []
    
    for atom in atoms:
      # Solvent.
      if (atom['chainID'] == 'A'):
        solvent_coordinates += [ atom['x'], atom['y'], atom['z'] ]
      # Solute.
      if (atom['chainID'] == 'B'):
        solute_coordinates += [ atom['x'], atom['y'], atom['z'] ]
    
    # Write coordinates.
    coordinates = solute_coordinates + solvent_coordinates
    writeAmberCrd(output_crd, coordinates, title = title, box_dimensions = box_dimensions, box_angles = box_angles)
  


if __name__ == '__main__':
    
    #=============================================================================================
    # MAIN
    #=============================================================================================
    # Create command-line argument options.
    usage_string = """
    Resolvate a given PDB file to produce a solvated AMBER .crd file, or a PDB file.
    Numbers of water molecules are determined from reference PDB file, and box volume from reference AMBER .crd file.
    The reference file for the box volume can also be a GROMACS *.gro file. 
    
    usage: %prog --refpdb REFERENCE.pdb --refcrd REFERENCE.crd --source SOURCE.pdb --output OUTPUT.crd --tol 0.6
    
    example: %prog --refpdb system.pdb --refcrd system.crd --source extracted.pdb --output resolvated.crd --tol 0.6
    example: %prog --refpdb system.pdb --refcrd system.gro --source extracted.pdb --output resolvated.pdb --tol 0.6
    """
    
    version_string = "%prog %__version__"
    
    parser = OptionParser(usage=usage_string, version=version_string)
    
    parser.add_option("-p", "--refpdb", metavar='REFERENCE.pdb',
                      action="store", type="string", dest='reference_pdb', default=None,
                      help="Reference PDB file (to determine number of waters).")
    parser.add_option("-c", "--refcrd", metavar='REFERENCE.crd',
                      action="store", type="string", dest='reference_crd', default=None,
                      help="Reference AMBER .crd file of GROMACS *.gro (to determine box dimensions).")
    parser.add_option("-s", "--source", metavar='SOURCE.pdb',
                      action="store", type="string", dest='source_pdb', default=None,
                      help="PDB file to resolvate.")
    parser.add_option("-o", "--output", metavar='OUTPUT.crd',
                      action="store", type="string", dest='output_crd', default=None,
                      help="Name of solvated AMBER .crd file to generate.")
    parser.add_option("-t", "--tol", metavar='TOLERANCE',
                      action="store", type="float", dest='tolerance', default=2.0,
                      help="Minimum distance (in Angstroms) allowed when packing. Optional." )
    # Parse command-line arguments.
    (options,args) = parser.parse_args()
    
    # Perform minimal error checking.
    if ((not options.reference_pdb) or (not options.reference_crd) or (not options.source_pdb) or (not options.output_crd)):
      parser.print_help()
      parser.error("All options must be specified.\n")
    
    reference_pdb = options.reference_pdb
    reference_crd = options.reference_crd
    source_pdb = options.source_pdb
    output_crd = options.output_crd
    tolerance = options.tolerance

    resolvate_from_template(reference_pdb, reference_crd, source_pdb, output_crd, tolerance=tolerance)

  

