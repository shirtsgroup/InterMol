
import pdb
import profile
from mmtools.Topology.System import *

#print "Reading in a system ...",
#structureList = ['LG2-0-eq.gro']
#topologyList = ['LG2.top']
#sys = System("Test", structureList, topologyList, "Gromacs")
#print "Done."

# Test a simplified structure and topology
print "Reading a structure and preprocessed topology ...",
structureList = 'FullStructure.gro'
topologyList = 'FullTopology.top'
sys = System("Test", structureList, topologyList, "Gromacs")
print "Done."

# Test, join, and split multiple structures and topologies
#print "Adding multiple structures and topologies ...",
#structureList = ['Protein.gro', 'Ligand.gro', 'Solvent.gro']
#topologyList = ['Protein.top', 'Ligand.top', 'Solvent.top']
#sys = System("Test", structureList, topologyList, "Gromacs")
#print "Done."

#print "Adding a molecule to a system ...",
#sys.addMolecule('Ion.gro', 'Ion.gro', "Gromacs")
#print "Done."

# Test a single topology and multiple structures
#print "Reading multiple structures and a single topology ...",
#structureList = ['LG2-0-eq.gro']
#topologyList = ['LG2.top']
#sys = System("Test", structureList, topologyList, "Gromacs")
#print "Done."

# Write out as a single pair of structure and topology files
print "Writing multiple structures and topologies into a single structure and topology file ...",
sys.writeToFiles('Gromacs')
print "Done."

# Write out as seperate structure and topology files
#print "Writing multiple structures and topologies into multiple structure and topology files ..."
#ligandList = ['Protein', 'Ligand', 'Solvent']
#structureFileList = ['Protein-out.gro', 'Ligand-out.gro', 'Solvent-out.gro']
#topologyFileList = ['Protein-out.top', 'Ligand-out.top', 'Solvent-out.gro']
#sys.writeToFiles('Gromacs', topologyName=toplogyFileList, structureName=structureFileList, molecules=LigandFileList)
#print "Done"

#pdb.set_trace()
