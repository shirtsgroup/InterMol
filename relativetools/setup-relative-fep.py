#=============================================================================================
# Set up relative alchemical free energy calculations in gromacs using AMBER leap.
#
# Written by John D. Chodera <jchodera@gmail.com> 2008-02-01
# Modified by Michael Zhu <mnz3v@virginia.edu> 2010-09-27
#=============================================================================================

#=============================================================================================
# IMPORTS
#=============================================================================================

from mmtools.moltools.ligandtools import *
from mmtools.moltools.relativefeptools import *
from mmtools.gromacstools.MdpFile import *
import pdb
import commands
import os
import os.path
import shutil

#=============================================================================================
# PARAMETERS
#=============================================================================================

# Path to receptor PDB file.
# All heavy atoms must be modeled -- hydrogens are ignored. Residue names are translated.
protein_pdb_filename = "/home/mrs5pt/preptools/mmtools/relativetools/Examples/JNK3/receptor-models/1JNK.automodel.mcce_out.pdb" # pdb of the protein without solvent and ligand

# Path to parameterized ligands.
ligand_basepath = '/home/mrs5pt/preptools/mmtools/relativetools/Examples/JNK3/ligands-parameterized'

# Base path for various scripts needed.
script_basepath = "../../"

# Target ligand group.
#ligand_name_list = [ ('jnk.aff-%d' % index) for index in (8, 10, 11, 15, 22, 26, 32, 52) ]
#ligand_name_list = [ ('jnk.aff-%d' % index) for index in (13, 16, 53) ]
ligand_name_list = [ ('jnk.aff-%d' % index) for index in (13, 16) ]

# Base path to set up free energy calculations
work_basepath = 'fep/'

#gromacs_converter = 'amb2gmx.pl'
#gromacs_converter = 'acpypi.py'
gromacs_converter = 'acpype.py'
   
#=============================================================================================
# GROMACS PARAMETERS - Edit these to change simulation protocol.
#=============================================================================================

# number of lambda states - should the the same as the # of cores.
cores = 16

# running expanded ensemble or Hamiltonian?
#runtype = 'expanded'
runtype = 'hamilton' 

# Phases of free energy calculations to set up.
phases = ['solvent', 'complex']

# Header block
mdpblock = dict()

mdpblock['header']  = """
; VARIOUS PREPROCESSING OPTIONS = 
define                   =
"""

# Construct minimizer block.
mdpblock['minimization'] = """
; RUN CONTROL PARAMETERS = 
integrator               = cg
nsteps                   = 1000

; ENERGY MINIMIZATION OPTIONS = 
; Force tolerance and initial step-size = 
emtol                    = 10
emstep                   = 1.0e-4
; Max number of iterations in relax_shells = 
niter                    = 20
; Step size (1/ps^2) for minimization of flexible constraints = 
fcstep                   = 0
; Frequency of steepest descents steps when doing CG = 
nstcgsteep               = 1000
"""

# dynamics control block
mdpblock['dynamics-stochastic'] = """
; RUN CONTROL PARAMETERS = 
integrator               = sd
; start time and timestep in ps = 
tinit                    = 0
dt                       = 0.002
nsteps                   = 1000000 ; 2 ns
; mode for center of mass motion removal = 
comm-mode                = Linear
; number of steps for center of mass motion removal = 
nstcomm                  = 500
; group(s) for center of mass motion removal = 
comm-grps                = 
"""

# dynamics control block
mdpblock['dynamics-verlet'] = """
; RUN CONTROL PARAMETERS = 
integrator               = md-vv
; start time and timestep in ps = 
tinit                    = 0
dt                       = 0.002
nsteps                   = 1000000 ; 2 ns
; mode for center of mass motion removal = 
comm-mode                = Linear
; number of steps for center of mass motion removal = 
nstcomm                  = 1
; group(s) for center of mass motion removal = 
comm-grps                = 
"""

# output control block
mdpblock['output'] = """
; OUTPUT CONTROL OPTIONS = 
; Output frequency for coords (x), velocities (v) and forces (f) = 
nstxout                  = 5000 ; (must be integral multiple of nstdgdl for post-analysis)
nstvout                  = 0
nstfout                  = 0
; Output frequency for energies to log file and energy file = 
nstlog                   = 500 ; 
nstenergy                = 5000 ;
; Output frequency and precision for xtc file = 
nstxtcout                = 100000 ;
xtc-precision            = 10000
; This selects the subset of atoms for the xtc file. You can = 
; select multiple groups. By default all atoms will be written. = 
;xtc-grps                 = solute
; Selection of energy groups = 
energygrps               = System
"""

# constraints block
mdpblock['constraints-lincs'] = """
; OPTIONS FOR BONDS     = 
constraints              = hbonds ; constrain bonds to hydrogen
; Type of constraint algorithm = 
constraint-algorithm     = lincs
; Use successive overrelaxation to reduce the number of shake iterations = 
Shake-SOR                = no
; Relative tolerance of shake = 
shake-tol                = 1e-12
; Highest order in the expansion of the constraint coupling matrix = 
lincs-order              = 4
; Lincs will write a warning to the stderr if in one step a bond = 
; rotates over more degrees than = 
lincs-warnangle          = 90
"""

# constraints block
mdpblock['constraints-shake'] = """
; OPTIONS FOR BONDS     = 
constraints              = hbonds ; constrain bonds to hydrogen
; Type of constraint algorithm = 
constraint-algorithm     = shake
; Use successive overrelaxation to reduce the number of shake iterations =
Shake-SOR                = no
; Relative tolerance of shake =
shake-tol                = 1e-12
; Highest order in the expansion of the constraint coupling matrix =
"""

mdpblock['restraints'] = """
;Enable dihedral and distance restraints
dihre=yes
dihre-fc=1
disre=yes
disre_fc=1
"""

# thermal and pressure control block
mdpblock['thermostat-nosehoover'] = """
; GENERATE VELOCITIES FOR STARTUP RUN = 
gen_vel                  = yes
gen_temp                 = 0; K - will be set later
gen_seed                 = 12345

; OPTIONS FOR WEAK COUPLING ALGORITHMS = 

; Temperature coupling   = 
Tcoupl                   = Nose-Hoover
; Groups to couple separately = 
tc-grps                  = System
; Time constant (ps) and reference temperature (K) = 
tau_t                    = 5.0
ref_t                    = 300 ; K - will be set later
"""

# thermal and pressure control block
mdpblock['thermostat-stochastic'] = """
; GENERATE VELOCITIES FOR STARTUP RUN = 

gen_vel                  = yes
gen_temp                 = 300 ; K - will be set later
gen_seed                 = 12345

tc-grps                  = System
; OPTIONS FOR WEAK COUPLING ALGORITHMS = 
; Time constant (ps) and reference temperature (K) = 
tau_t                    = 1.0
ref_t                    = 300 ; K - will be set later
"""

# thermal and pressure control block
mdpblock['thermostat-equilibration'] = """
; GENERATE VELOCITIES FOR STARTUP RUN = 
gen_vel                  = yes
gen_temp                 = 300 ; K - will be set later
gen_seed                 = 12345

; OPTIONS FOR WEAK COUPLING ALGORITHMS = 

; Temperature coupling   = 
Tcoupl                   = Berendsen
; Groups to couple separately = 
tc-grps                  = System
; Time constant (ps) and reference temperature (K) = 
tau_t                    = 0.2
ref_t                    = 300 ; K - will be set later
"""

mdpblock['barostat'] = """
; Pressure coupling      = 
Pcoupl                   = MTTK
Pcoupltype               = isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar) = 
tau_p                    = 5.0 ; ps
compressibility          = 4.5e-5 ; 1/bar
ref_p                    = 1.0 ; bar
"""

mdpblock['barostat-equilibration'] = """
; Pressure coupling      = 
Pcoupl                   = berendsen
Pcoupltype               = isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar) = 
tau_p                    = 5.0 ; ps
compressibility          = 4.5e-5 ; 1/bar
ref_p                    = 1.0 ; bar
"""

mdpblock['nonbonded-solvent'] = """
; NEIGHBORSEARCHING PARAMETERS = 
; nblist update frequency = 
nstlist                  = 10
; ns algorithm (simple or grid) = 
ns_type                  = grid
; Periodic boundary conditions: xyz or no = 
pbc                      = xyz
; nblist cut-off         = 
rlist                    = 1.4

; OPTIONS FOR ELECTROSTATICS AND VDW = 
; Method for doing electrostatics = 
coulombtype              = PME
rcoulomb                 = 1.4
; Dielectric constant (DC) for cut-off or DC of reaction field = 
epsilon-r                = 1
; Method for doing Van der Waals = 
vdw-type                 = switch
; cut-off lengths        = 
rvdw-switch              = 1.25
rvdw                     = 1.3
; Apply long range dispersion corrections for Energy and Pressure = 
DispCorr                 = AllEnerPres
; Spacing for the PME/PPPM FFT grid = 
fourierspacing           = 0.10
; FFT grid size, when a value is 0 fourierspacing will be used = 
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
; EWALD/PME/PPPM parameters = 
pme_order                = 6
ewald_rtol               = 1e-06
ewald_geometry           = 3d
epsilon_surface          = 0
optimize_fft             = yes
"""

mdpblock['nonbonded-vacuum'] = """
; NEIGHBORSEARCHING PARAMETERS = 
; nblist update frequency = 
nstlist                  = 10 ; can we set this to 0 for vacuum?
; ns algorithm (simple or grid) = 
ns_type                  = simple
; Periodic boundary conditions: xyz or no = 
pbc                      = no
; nblist cut-off         = 
rlist                    = 25.0 ; to emulate no cutoff
domain-decomposition     = no

; OPTIONS FOR ELECTROSTATICS AND VDW = 
; Method for doing electrostatics = 
coulombtype              = cut-off
rcoulomb-switch          = 0
rcoulomb                 = 23.0 ; to emulate no cutoff
; Dielectric constant (DC) for cut-off or DC of reaction field = 
epsilon-r                = 1
; Method for doing Van der Waals = 
vdw-type                 = cut-off
; cut-off lengths        = 
rvdw                     = 23.0 ; to emulate no cutoff
; Apply long range dispersion corrections for Energy and Pressure = 
DispCorr                 = no
""" 

# free energy blocks
mdpblock['free-energy-common-relative'] = """
; OPTIONS FOR FIXED ENSEMBLE SIMULATIONS
; Free energy control stuff = 
free-energy               = yes 
sc-power                  = 1
sc-alpha                  = 0.5
fep-lambdas               = 0.0 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 0.0 0.0  0.0 0.0  0.0  0.0 0.0
coul-lambdas              = 0.0 0.13 0.25 0.5  1.0  1.0  1.0  1.0  1.0 1.0 1.0  1.0 1.0  1.0  1.0 1.0
vdw-lambdas               = 0.0 0.0  0.0  0.0  0.0  0.2  0.3  0.4  0.5 0.6 0.65 0.7 0.75 0.8  0.9 1.0
;bonded-lambdas            = 0.0 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 0.0 0.0  0.0 0.0  0.0  0.0 0.0
;restraint-lambdas         = 0.0 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 0.0 0.0  0.0 0.0  0.0  0.0 0.0
;init-lambda-weights       = 0.0 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 0.0 0.0  0.0 0.0  0.0  0.0 0.0
""" 

mdpblock['free-energy-fixed'] = """
lmc-stats                 = no
lmc-mc-move               = no
lmc-weights-equil         = no
dhdl-print-energy         = yes
nstdhdl                   = 100
nstfep                    = 20
""" 

mdpblock['free-energy-expanded'] = """
lmc-stats                 = wang-landau
lmc-mc-move               = gibbs
dhdl-print-energy         = yes
init-wl-delta             = 0.25
wl-scale                  = 0.85
wl-ratio                  = 0.85
lmc-weights-equil         = no

;lmc-weights-equil         = number-all-lambda
;lmc-weights-equil         = number-samples
;lmc-weights-equil         = number-steps
;lmc-weights-equil         = count-ratio
;lmc-weights-equil         = wl-delta

;weight-equil-number-all-lambda  =  50
;weight-equil-number-samples     =  1000
;weight-equil-number-steps       =  600
;weight-equil-count-ratio        =  0.80
;weight-equil-wl-delta           =  0.001

; schedule for switching off lambdas
; first, restraints are turned on as charges are switched off
; next, vdw and torsions are switched off
nstdhdl                   = 100
nstfep                    = 20
""" 
# mutate is default
mdpblock['mutate'] = """

"""

# free energy block for gromacs 4.5
mdpblock['decouple'] = """
couple-moltype           = MOL
couple-lambda0           = vdw-q
couple-lambda1           = no
couple-intramol          = no
""" % vars()

#=============================================================================================
def compose_blocks(blocknames):

    """Compose desired blocks into a list of strings

    ARGUMENTS
      blocknames (list of strings) - list of names of blocks to be concatenated

    RETURNS
      lines (list of strings) - concatenated blocks forming the contents of a gromacs .mdp file
    """

    contents = ""

    for blockname in blocknames:
        contents += mdpblock[blockname]

    lines = contents.split('\n')

    return lines

#=============================================================================================
# SUBROUTINES
#=============================================================================================
def write_file(filename, contents):
    """Write the specified contents to a file.

    ARGUMENTS
      filename (string) - the file to be written
      contents (string or list of strings) - the contents of the file to be written

    """

    outfile = open(filename, 'w')

    if type(contents) == list:
        for line in contents:
            outfile.write(line + "\n")
    elif type(contents) == str:
        outfile.write(contents)
    else:
        raise "Type for 'contents' not supported: " + repr(type(contents))

    outfile.close()

    return

def perturbGromacsTopologyToIntermediate(topology_filename, perturbed_resname, ligand, common_substructure, perturb_torsions = True):
   """Modify a gromacs topology file to add perturbed-state parameters to transform molecule into a common intermediate.

   ARGUMENTS
     topology_file (string) - the name of the Gromacs topology file to modify
     perturbed_resname (string) - residue name of residue in topology file to perturb (must be unique)
     ligand (OEMol) - ligand to be perturbed -- must be the same one used to generate the topology file with parameterizeForGromacs(), but with GAFF atom/bond types
     common_substructure (OEMol) - the intermediate to perturb to, with charges assigned and AMBER atom and bondtypes.

   OPTIONAL ARGUMENTS
     perturb_torsions (boolean) - if True, torsions containing at least two atoms in the annihilated region whose central bond is not in an aromatic ring will be turned off in B state (default: True)

   NOTES
     The common_intermediate MUST have AMBER atom and bondtypes.
     This code currently only handles the special format gromacs topology files produced by amb2gmx.pl -- there are allowed variations in format that are not treated here.
     Note that this code also only handles the first section of each kind found in a gromacs .top file.
     MRS: Moved to acpype.py
     
   TODO
     Perhaps this method should be combined with 'parameterizeForGromacs' as an optional second step, ensuring that the correct 'molecule' is used.
     Generalize this code to allow it to operate more generally on a specified moleculetype.

   """

   # DEFINE INTERNAL METHODS
   def stripcomments(line):
      """Return line with whitespace and comments stripped.
      """
      # strip comments
      index = line.find(';')
      if index > -1:
         line = line[0:index]
      # strip whitespace
      line = line.strip()
      # return stripped line
      return line

   def extract_section(lines, section):
      """Identify lines associate with a section.

      ARGUMENTS
        lines (list of strings) - the lines in the file
        section (string) - the section name to locate

      RETURNS
        indices (list of integers) - line indices within lines belonging to section
      """

      indices = list()

      nlines = len(lines)
      for start_index in range(nlines):
         # get line
         line = stripcomments(lines[start_index])
         # split into elements
         elements = line.split()
         # see if keyword is matched
         if (len(elements) == 3):
            if (elements[0]=='[') and (elements[1]==section) and (elements[2]==']'):
               # increment counter to start of section data and abort search
               start_index += 1
               break

      # throw an exception if section not found
      if (start_index == nlines):
         raise "Section %(section)s not found." % vars()

      # Locate end of section.
      for end_index in range(start_index, nlines):
         # get line
         line = stripcomments(lines[end_index])
         # split into elements
         elements = line.split()
         # see if keyword is matched
         if (len(elements) == 3):
            if (elements[0]=='['):
               break

      # compute indices of lines in section
      indices = range(start_index, end_index)

      # return these indices
      return indices
   # END INTERNAL METHODS

   # find correspondence between common substructure and ligand atoms
   mcss = OEMCSSearch(common_substructure, OEExprOpts_StringType, OEExprOpts_IntType)
   min_atoms = 4
   mcss.SetMinAtoms(min_atoms)

   # show common intermediate
   print "common_substructure:"
   for atom in common_substructure.GetAtoms():
       print "%6s %6s %6.3f" % (atom.GetName(), atom.GetType(), atom.GetPartialCharge())
   print "ligand:"
   for atom in ligand.GetAtoms():
       print "%6s %6s %6.3f" % (atom.GetName(), atom.GetType(), atom.GetPartialCharge())

   # store atomname -> charge and atomtype -> typename correspondence when atoms are present in substructure
   mutated_charges = dict() # map of charges
   retained_atoms = list() # list of which atoms are not turned into dummies
   perturbed_atomtypes = list() # list of perturbed atom types

   # Fix to put all atomtypes into perturbed_atomtypes because it can't hurt to have all of them in there--MZ
   for atom in ligand.GetAtoms():
       perturbed_atomtypes.append(atom.GetType())
   # End of the fix.

   # extract first match
   for match in mcss.Match(ligand):
       print "number of atoms in match: %d" % match.NumAtoms()
       for matchpair in match.GetAtoms():
           atomname = matchpair.target.GetName() # name of ligand atom
           charge = matchpair.pattern.GetPartialCharge() # charge on common substructure
           atomtype = matchpair.pattern.GetType() # GAFF atom type
           mutated_charges[atomname] = charge
           retained_atoms.append(atomname)
           #perturbed_atomtypes.append(atomtype)--MZ
       break
   print "mutated_charges"
   print mutated_charges
   print "retained_atoms"
   print retained_atoms
   print "perturbed_atomtypes"
   print perturbed_atomtypes

   # Read the contents of the topology file.
   lines = read_file(topology_filename)

   # Parse [ atomtypes ]
   atomtypes = list() # storage for atom types
   indices = extract_section(lines, 'atomtypes')
   for index in indices:
      # extract the line
      line = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip line if number of elements is less than expected
      if (nelements < 7): continue
      # parse elements
      atomtype = dict()
      atomtype['name'] = elements[0]
      atomtype['bond_type'] = elements[1]
      atomtype['mass'] = float(elements[2])
      atomtype['charge'] = float(elements[3])
      atomtype['ptype'] = elements[4]
      atomtype['sigma'] = float(elements[5])
      atomtype['epsilon'] = float(elements[6])
      # append
      atomtypes.append(atomtype)

   # Augment [ atomtypes ] list with any needed dummy atomtypes with zero LJ well depth.
   indices = extract_section(lines, 'atomtypes')
   for atomtype in atomtypes:
       if atomtype['name'] in perturbed_atomtypes:
           # make perturbed record
           perturbed_atomtype = atomtype
           # modify name and LJ well depth, leaving other quantities unperturb
           perturbed_atomtype['name'] += '_dummy'
           perturbed_atomtype['epsilon'] = 0.0
           # form the line
           line = "%(name)-10s%(bond_type)6s      0.0000  0.0000  A %(sigma)13.5e%(epsilon)13.5e ; perturbed\n" % perturbed_atomtype
           # insert the new line
           lines.insert(indices[-1], line)
           indices.append(indices[-1]+1)

   # Process [ atoms ] section
   atoms = list()
   atom_indices = dict()
   perturb_atom_indices = list() # list of atoms that are to be perturbed
   indices = extract_section(lines, 'atoms')
   for index in indices:
      # extract the line
      line = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 8): continue
      # parse line
      atom = dict()
      atom['nr'] = int(elements[0])
      atom['type'] = elements[1]
      atom['resnr'] = int(elements[2])
      atom['residue'] = elements[3]
      atom['atom'] = elements[4]
      atom['cgnr'] = int(elements[5])
      atom['charge'] = float(elements[6])
      atom['mass'] = float(elements[7])

      # skip if not in speicifed residue
      if not atom['residue'] == perturbed_resname: continue

      # set perturbation type
      atomname = atom['atom']
      # atom type
      if atomname in retained_atoms:
          atom['typeB'] = atom['type'] # atom type is retained
          atom['chargeB'] = mutated_charges[atomname] # common substructure charge is used
          atom['comment'] = 'common substructure'
      else:
          atom['typeB'] = atom['type'] + "_dummy" # atom is turned into a dummy
          atom['chargeB'] = 0.0 # atom is discharged
          atom['comment'] = 'annihilated'
          perturb_atom_indices.append(atom['nr']) # add atom index to list of atoms that will be perturbed

      # construct a new line
      line = "%(nr)6d %(type)10s %(resnr)6d %(residue)6s %(atom)6s %(cgnr)6d %(charge)10.5f %(mass)10.6f %(typeB)10s %(chargeB)10.5f %(mass)10.6f ; %(comment)s\n" % atom

      # replace the line
      lines[index] = line

      # store atoms
      atoms.append(atom)

      # store lookup of atom atom names -> atom numbers
      atom_indices[atom['atom']] = atom['nr']

   # Process [ bonds ] section
   indices = extract_section(lines, 'bonds')
   for index in indices:
      # extract the line
      line = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 5): continue
      # parse line
      bond = dict()
      bond['i'] = int(elements[0])
      bond['j'] = int(elements[1])
      bond['function'] = int(elements[2])
      bond['Req'] = float(elements[3])
      bond['Keq'] = float(elements[4])
      # skip if an atom range is specified, and this is not among the atoms we're looking for
      if perturb_atom_indices:
         if (bond['i'] not in perturb_atom_indices) and (bond['j'] not in perturb_atom_indices): continue
      # construct a new line
      line = "%(i)5d %(j)5d %(function)5d%(Req)12.4e%(Keq)12.4e" % bond
      line += " %(Req)12.4e%(Keq)12.4e\n" % bond
      # replace the line
      lines[index] = line

   # Process [ angles ] section
   indices = extract_section(lines, 'angles')
   for index in indices:
      # extract the line
      line = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 6): continue
      # parse line
      angle = dict()
      angle['i'] = int(elements[0])
      angle['j'] = int(elements[1])
      angle['k'] = int(elements[2])
      angle['function'] = int(elements[3])
      angle['theta'] = float(elements[4])
      angle['cth'] = float(elements[5])
      # skip if an atom range is specified, and this is not among the atoms we're looking for
      if perturb_atom_indices:
         if (angle['i'] not in perturb_atom_indices) and (angle['j'] not in perturb_atom_indices) and (angle['k'] not in perturb_atom_indices): continue
      # construct a new line
      line = "%(i)5d %(j)5d %(k)5d %(function)5d%(theta)12.4e%(cth)12.4e" % angle
      line += " %(theta)12.4e%(cth)12.4e\n" % angle
      # replace the line
      lines[index] = line

   # Set rotatable bond torsions in B state to zero, if desired.
   if perturb_torsions:
      # Determine list of rotatable bonds to perturb.
      rotatable_bonds = list()
      for bond in ligand.GetBonds():
         # This test is used because bond.IsRotor() doesn't seem to work correctly (e.g. phenol).
         if (not bond.IsAromatic()) and (bond.GetOrder() == 1) and (bond.GetBgn().GetDegree() > 1) and (bond.GetEnd().GetDegree() > 1):
            i = atom_indices[bond.GetBgn().GetName()]
            j = atom_indices[bond.GetEnd().GetName()]
            if (i < j):
               rotatable_bonds.append( (i,j) )
            else:
               rotatable_bonds.append( (j,i) )

      print "%d rotatable bond(s) found." % len(rotatable_bonds)

      # Search for [ dihedrals ] section.
      indices = extract_section(lines, 'dihedrals') # extract non-blank, non-comment lines
      for index in indices:
         # extract the line
         line = stripcomments(lines[index])
         # parse the line
         elements = line.split()
         nelements = len(elements)
         # skip unrecognized lines
         if (nelements < 11): continue
         # extract dihedral atom indices
         i = int(elements[0])
         j = int(elements[1])
         k = int(elements[2])
         l = int(elements[3])

         # skip if an atom range is specified, and this is not among the atoms we're looking for
         if perturb_atom_indices:
            if (i not in perturb_atom_indices) and (j not in perturb_atom_indices) and (k not in perturb_atom_indices) and (l not in perturb_atom_indices): continue

         # function number
         function = int(elements[4])
         if function != 3: raise "Only [dihedrals] function = 3 is supported."
         # C0 - C5 parameters
         C = list()
         for element in elements[5:11]:
            C.append(float(element))

         # reconstruct perturbed line
         line = "    %-4s %-4s %-4s %-4s %3d%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f" % (i, j, k, l, function, C[0], C[1], C[2], C[3], C[4], C[5])
         if (j,k) in rotatable_bonds:
            # perturb rotatable bonds
            line += " %12.5f%12.5f%12.5f%12.5f%12.5f%12.5f ; perturbed" % (0.0, 0.0, 0.0, 0.0, 0.0, 0.0) + "\n"
         else:
            # don't perturb 
            line += " %12.5f%12.5f%12.5f%12.5f%12.5f%12.5f" % (C[0], C[1], C[2], C[3], C[4], C[5]) + "\n"

         # replace the line
         lines[index] = line

   # Replace topology file.
   outfile = open(topology_filename, 'w')
   for line in lines:
      outfile.write(line)
   outfile.close()

   return
#=============================================================================================
def rename_pdb_for_amber(protein_pdb_filename, protein_amberpdb_filename):
    """Rename residues in a PDB file to match AMBER naming convention, omitting atoms AMBER will not recognize.

    ARGUMENTS
      protein_pdb_filename (string) - name of source PDB file to be read
      protein_amberpdb_filename (string) - name of PDB file to be written with AMBER naming conventions

    NOTE
      This conversion is still jnk3-specific at this time.
      Four-letter residue names aren't properly parsed, but a hack is in.

    """

    # read the PDB file
    lines = read_file(protein_pdb_filename)

    # allocate storage for target file contents
    outlines = list()

    # parse PDB file, creating new contents
    for line in lines:
        # only ATOM records will be retained
        if line[0:6] == "ATOM  ":
            # Parse line into fields.
            serial = line[6:11]
            name = line[12:16]
            altLoc = line[16:17]
            resname = line[17:21]
            chainid = line[21:22]
            resseq = line[22:26]
            idcode = line[26:27]
            remainder = line[27:]

            # fail without writing PDB file if early NTR or CTR detected
            if (resname == 'NTR ') or (resname == 'CTR '):
                print "Early chain break detected. Failing."
                return

            # drop the line it is a hydrogen
            if name[1] == 'H': continue
            # phosphorylated residues have a strange naming convention on their hydrogens
            if name[0] == 'H': continue

            # make residue name substitutions
            if resname == 'GYN ': resname = 'GUC '
            if resname == 'GLX ': resname = 'GLU '
            if resname == 'LYP ': resname = 'LYS '
            if resname == 'CYN ': resname = 'CYS '
            if resname == 'NSER': resname = 'SER '
            if (resname == 'CGLU') or (resname == 'CGLH'):
                resname = 'GLU '
                # rename terminal oxygens
                if name == ' OC1' or name == ' OC ':
                    name = ' OXT'
                if name == ' OC2':
                    name = ' O  '

            # renumber serial
            serial = '% 5d' % (len(outlines) + 1)

            # reconstitute line
            outline = 'ATOM  ' + serial + ' ' + name + altLoc + resname + chainid + resseq + idcode + remainder
            outlines.append(outline)

    write_file(protein_amberpdb_filename, outlines)

    return

#=============================================================================================
def writeBatchScript(path,jobname,runtype,cores,restraints):
    """        
    path (string) - absolute path for the script to be printed to
    runtype (string) - 'expanded' or 'hamilton'
    jobname (string) - jobname
    restraints (bool) - True if we are writing restraints, false otherwise.  Should only be true for restraints?

    """

    import re

    if (runtype == 'hamilton'):
        binpath = '/h3/n1/shirtsgroup/gromacs_4.5plus/MPI/bin'
    if (runtype == 'expanded'):    
        binpath = '/h3/n1/shirtsgroup/gromacs_4.5plus/NOMPI/bin'

    # write batch queue script
    print "Writing batch queue script..."
    path = os.path.abspath(path)

    cores_per_node = 8
    batchcommands = []
    batchcommands.append('#!/bin/sh')
    name = runtype[0:3] + jobname
    maxlen = min(15,len(name))
    name = name[0:maxlen]
    batchcommands.append('#PBS -N ' + name)
    batchcommands.append('#PBS -o ' + name + '.out')
    batchcommands.append('#PBS -e ' + name + '.err')

    if (runtype == 'hamilton'):
        nodes = cores/cores_per_node+1
    if (runtype == 'expanded'):
        nodes = 1

    batchcommands.append('#PBS -l select=' + str(nodes) + ':mpiprocs=' + str(cores_per_node) + ':ncpus=' + str(cores_per_node))

    general_setup = """\
#PBS -q shirtscluster
#PBS -W group_list=shirtsgroup
#PBS -r n                       
#PBS -l walltime=1536:00:00
source $GMXRC
module load MODULE
cat $PBS_NODEFILE
set MPI_THREADS = `cat $PBS_NODEFILE | wc -l`
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

NP=`wc -l < $PBS_NODEFILE`
NN=`sort -u $PBS_NODEFILE | wc -l`
echo Number of nodes is $NN
echo Number of processors is $NP
    """%vars()

    if (runtype == 'hamilton'): 
        general_setup = re.sub('MODULE ','mvapich2-intel',general_setup)
    if (runtype == 'expanded'):
        general_setup = re.sub('MODULE ','intel',general_setup)
        
    if (restraints):
        commands = """\
# create a tpr file from which to create restraints
GROMPP -f minimize.mdp -c system_GMX.gro -p system_GMX.top -o minimize.tpr -maxwarn 100 
# Integrate restraints into pre-minimization topfile
python compute_angles.py -f system_GMX.gro -s minimize.tpr -d $RANSEED
python restraint_topology.py -n system -p .

# constrained minimize with restraints
GROMPP -f minimize.mdp -c system_restr.gro -p system_restr.top -o minimize.tpr -maxwarn 100 
MDRUN -s minimize.tpr -x minimize.xtc -c minimize.gro -e minimize.edr -g minimize.log

# create a tpr file from which to create restraints
GROMPP -f equilibration.mdp -c minimize.gro -p system_restr.top -o equilibration.tpr -maxwarn 100 
# Integrate restraints into pre-equilibration topfile
cp system_restr.top minimize.top
#python compute_angles.py -f minimize.gro -s equilibration.tpr -d $RANSEED
python restraint_topology.py -n minimize -p .

# equilibration with restraints
GROMPP -f equilibration.mdp -c minimize_restr.gro -p minimize_restr.top -o equilibration.tpr -maxwarn 100 
MDRUN -s equilibration.tpr -o equilibration.trr -x equilibration.xtc -c equilibration.gro -e equilibration.edr -g equilibration.log

# create a TPR file, integrate restraints
GROMPP -f production.mdp -c equilibration.gro -p minimize_restr.top -o production.tpr -maxwarn 100 
cp minimize_restr.top equilibration.top
#python compute_angles.py -f equilibration.gro -s production.tpr -d $RANSEED
python restraint_topology.py -n equilibration -p .

# production with restraints
GROMPP -f production.mdp -c equilibration_restr.gro -p equilibration_restr.top -o production.tpr -maxwarn 100 
MDRUN -s production.tpr -o production.trr -x production.xtc -c production.gro -e production.edr -g production.log

# signal completion
mv run.sh run.sh.done
       """ %vars()
    else:
       commands = """\
# constrained minimize
GROMPP -f minimize.mdp -c system_GMX.gro -p system.top -o minimize.tpr -maxwarn 100 
MDRUN -s minimize.tpr -x minimize.xtc -c minimize.gro -e minimize.edr -g minimize.log

# equilibration
GROMPP -f equilibration.mdp -c minimize.gro -p system_GMX.top -o equilibration.tpr -maxwarn 100 
MDRUN -s equilibration.tpr -o equilibration.trr -x equilibration.xtc -c equilibration.gro -e equilibration.edr -g equilibration.log

# production
GROMPP -f production.mdp -c equilibration.gro -p system_GMX.top -o production.tpr -maxwarn 100 
    """ % vars()

    if (runtype == 'hamilton'):
        grompp = binpath + '/grompp_mpi_d '
        mdrun = 'mpiexec -comm mpich2-pmi' + binpath + '/mdrun_mpi_d -np ' + str(cores_per_node)
    if (runtype == 'expanded'):    
        grompp = binpath + '/grompp_d '
        mdrun = binpath + '/mdrun_d -nt ' + str(cores_per_node)

    commands = re.sub('GROMPP',grompp,commands)
    commands = re.sub('MDRUN',mdrun,commands)

    batchcommands.append(general_setup)
    batchcommands.append(commands)

    if (runtype == 'hamilton'): 
        for i in range(cores):
            grompp_run = grompp + ' -f production.' + str(i) + '.mdp -c equilibration.gro -p system_GMX.top -o production.' + str(i) + '.tpr -maxwarn 1000'            
            batchcommands.append(grompp_run)

    elif (runtype == 'expanded'):
        grompp_run = grompp + ' -f production.mdp -c equilibration.gro -p system_GMX.top -o production.tpr -maxwarn 100'
        batchcommands.append(grompp_run)

    if (runtype == 'hamilton'): 
        # eventually, replace the 20 with a multiple of nstfep to avoid issues.
        mdrun_run = mdrun + ' -multi ' + str(cores) + ' -replex 20 -s production..tpr -o production..trr -x production..xtc -c production..gro -e production..edr -g production..log -dhdl production.dhdl..xvg'

    elif (runtype == 'expanded'):
        mdrun_run = mdrun + ' -s production.tpr -o production.trr -x production.xtc -c production.gro -e production.edr -g production.log -dhdl production.dhdl.xvg'

    batchcommands.append(mdrun_run)

    write_file(os.path.join(path, 'run.sh'), batchcommands)

#=============================================================================================
def writeMdpFiles(path, runtype, cores, transform):
    """
    path (string) - file path for these mdpfiles
    runtype (string) - either 'expanded' or 'hamiltonian'
    cores (int) - number of cores that will be run on
    transform (string,optional) - 'decouple' or 'mutate'

    """

    # Constrained Minimization
    mdpfile = MdpFile(compose_blocks(['header', 'minimization', 'nonbonded-solvent', 'constraints-lincs', 'restraints', 'output'])) # 
    mdpfile.write(os.path.join(path, 'minimize.mdp'))

    # Equilibration
    mdpfile = MdpFile(compose_blocks(['header', 'dynamics-stochastic', 'nonbonded-solvent', 'constraints-lincs', 'thermostat-stochastic', 'barostat-equilibration', 'output']))
    mdpfile.setParameter('nsteps', '10000') # 20 ps equilibration
    mdpfile.randomizeSeed() # randomize velocities
    mdpfile.write(os.path.join(path, 'equilibration.mdp'))
    
    # Production
    if (runtype == 'expanded'):
        mdpfile = MdpFile(compose_blocks(['header', 'dynamics-stochastic', 'nonbonded-solvent', 'constraints-lincs', 'thermostat-stochastic', 'free-energy-common-relative', 'free-energy-expanded', transform, 'output']))
        mdpfile.randomizeSeed() # randomize velocities
        mdpfile.setTemperature(300);
        mdpfile.setParameter('init-lambda-state',0)
        mdpfile.write(os.path.join(path, 'production.mdp'))

    if (runtype == 'hamiltonian'):
        for lam in range(cores):  
            mdpfile = MdpFile(compose_blocks(['header', 'dynamics-stochastic', 'nonbonded-solvent', 'constraints-lincs', 'thermostat-stochastic', 'free-energy-common-relative', 'free-energy-fixed', 'output']))
            mdpfile.randomizeSeed() # randomize velocities
            mdpfile.setTemperature(300);
            mdpfile.setParameter('init-lambda-state',lam)
            mdpfile.write(os.path.join(solvent_path, 'production.',lam,'.mdp'))

#=============================================================================================
def setup_solvent_simulation(solvent_path, jobname, runtype, ligand, cores, ligand_off_filename, ligand_frcmod_filename, common_substructure):
    """Set up a solvent simulation.

    ARGUMENTS
      solvent_path (string) - the pathname where the simulation is to be set up (created if doesn't exist).
      jobname (string) - job name to be used to form batch queue job name
      runtype (string) - 'expanded' or 'hamilton'
      ligand (OEMol) - the ligand
      cores (int) - # of cores to run on
      ligand_off_filename (string) - leap library for ligand
      ligand_frcmod_filename (string) - additional ligand parameters
      common_substructure (OEMol) - 
    """

    print "\nPREPARING SOLVATED LIGAND"

    # create directory if it doesn't exist
    if not os.path.exists(solvent_path):
        os.makedirs(solvent_path)

    # get ligand formal charge
    ligand_charge = formalCharge(ligand)

    # set up the system with tLEaP
    print "Solvating the ligand with tleap..."
    system_prmtop_filename = os.path.join(solvent_path,'system.prmtop')
    system_crd_filename = os.path.join(solvent_path,'system.crd')
    tleap_input_filename = os.path.join(solvent_path, 'setup-system.leap.in')
    tleap_output_filename = os.path.join(solvent_path, 'setup-system.leap.out')
    clearance = 10.0 # clearance in A
    contents = """
source leaprc.ff99SB
source leaprc.gaff

# load antechamber-generated additional parameters
mods = loadAmberParams %(ligand_frcmod_filename)s

# Load ligand.
loadOff %(ligand_off_filename)s

# Create system.
system = combine { ligand }
""" % vars()
    # add counterions
    if (ligand_charge != 0):
        nions = abs(ligand_charge)
        if ligand_charge < 0: iontype = 'Na+'
        if ligand_charge > 0: iontype = 'Cl-'
        contents += """
# Add counterions.
addions system %(iontype)s %(nions)d
""" % vars()
    #
    contents += """
# Solvate in truncated octahedral box.
# solvateOct system TIP3PBOX %(clearance)f
# solvate in rectilinear box
solvateBox system TIP3PBOX %(clearance)f

# Check the system
check system

# Write the system
saveamberparm system %(system_prmtop_filename)s %(system_crd_filename)s
""" % vars()
    write_file(tleap_input_filename, contents)
    command = 'tleap -f %(tleap_input_filename)s > %(tleap_output_filename)s' % vars()
    output = commands.getoutput(command)

    # convert to gromacs
    print "Converting to gromacs..."

    converter = os.path.join(os.getenv('MMTOOLSPATH'), 'converters', gromacs_converter)
    system_prefix = os.path.join(solvent_path, 'system')
    current_path = os.getcwd()
    os.chdir(solvent_path)    

    if gromacs_converter == 'amb2gmx.pl':
        command = '%(converter)s --parmtop %(system_prmtop_filename)s --crd %(system_crd_filename)s --outname %(system_prefix)s' % vars()
    if gromacs_converter == 'acpype.py':    
        command = 'python2.6 %(converter)s -p system.prmtop -x system.crd -b system -o gmx' % vars()
    print command
    output = commands.getoutput(command)
    print output
    os.chdir(current_path)

    # set up perturbation
    print "Modifying topology file for perturbation..."
    system_top_filename = os.path.join(solvent_path, 'system_GMX.top')
    perturbGromacsTopologyToIntermediate(system_top_filename, 'MOL', ligand, common_substructure, perturb_torsions = True)

    # set up mdp files
    print "Writing mdp files..."
    writeMdpFiles(solvent_path, runtype, cores,'mutate')

    # write batch scripts
    print "Writing batch scripts..."
    writeBatchScript(solvent_path,jobname,runtype,cores,restraints=False)

    return
#=============================================================================================
def setup_vacuum_simulation(vacuum_path, jobname, runtype, ligand, cores, ligand_off_filename, ligand_frcmod_filename, common_substructure):
    """Set up a vacuum simulation.

    ARGUMENTS
      vacuum_path (string) - the pathname where the simulation is to be set up (created if doesn't exist).
      jobname (string) - job name to be used to form batch queue job name
      runtype (string) - 'expanded' or 'hamilton'
      ligand (OEMol) - the ligand
      ligand_off_filename (string) - leap library for ligand
      ligand_frcmod_filename (string) - additional ligand parameters
      common_substructure (OEMol) - 
    """

    print "\nPREPARING LIGAND IN VACUUM"

    # create directory if it doesn't exist
    if not os.path.exists(vacuum_path):
        os.makedirs(vacuum_path)

    # get ligand formal charge
    ligand_charge = formalCharge(ligand)

    # set up the system with tLEaP
    print "Solvating the ligand with tleap..."
    system_prmtop_filename = os.path.join(solvent_path,'system.prmtop')
    system_crd_filename = os.path.join(solvent_path,'system.crd')
    tleap_input_filename = os.path.join(solvent_path, 'setup-system.leap.in')
    tleap_output_filename = os.path.join(solvent_path, 'setup-system.leap.out')
    clearance = 10.0 # clearance in A
    contents = """
source leaprc.ff99SB
source leaprc.gaff

# load antechamber-generated additional parameters
mods = loadAmberParams %(ligand_frcmod_filename)s

# Load ligand.
loadOff %(ligand_off_filename)s

# Create system.
system = combine { ligand }
""" % vars()

    contents += """
# Check the system
check system

# Write the system
saveamberparm system %(system_prmtop_filename)s %(system_crd_filename)s
""" % vars()
    write_file(tleap_input_filename, contents)
    command = 'tleap -f %(tleap_input_filename)s > %(tleap_output_filename)s' % vars()
    output = commands.getoutput(command)

    # convert to gromacs
    print "Converting to gromacs..."

    converter = os.path.join(os.getenv('MMTOOLSPATH'), 'converters', gromacs_converter)
    system_prefix = os.path.join(vacuum_path, 'system')
    current_path = os.getcwd()
    os.chdir(solvent_path)    

    if gromacs_converter == 'amb2gmx.pl':
        command = '%(converter)s --parmtop %(system_prmtop_filename)s --crd %(system_crd_filename)s --outname %(system_prefix)s' % vars()
    if gromacs_converter == 'acpype.py':    
        command = 'python2.6 %(converter)s -r -p system.prmtop -x system.crd -b system -o gmx' % vars()
    print command
    output = commands.getoutput(command)
    print output
    os.chdir(current_path)

    # set up perturbation
    print "Modifying topology file for perturbation..."
    system_top_filename = os.path.join(vacuum_path, 'system_GMX.top')
    perturbGromacsTopologyToIntermediate(system_top_filename, 'MOL', ligand, common_substructure, perturb_torsions = True)

    # set up mdp files
    print "Writing mdp files..."
    writeMdpFiles(vacuum_path, runtype, cores,'mutate')

    # write batch scripts
    print "Writing batch scripts..."
    writeBatchScript(vacuum_path,jobname,runtype,cores,restraints=False)

    return

#=============================================================================================
def setup_complex_simulation(complex_path, jobname, runtype, ligand, cores, ligand_off_filename, ligand_frcmod_filename, protein_off_filename, protein_charge, common_substructure):
    """Set up free energy calculation for ligand in complex.

    ARGUMENTS
      complex_path (string) - the pathname where the simulation is to be set up (created if doesn't exist).
      jobname (string) - job name to be used to form batch queue job name
      runtype (string) - 'expanded' or 'hamilton' 
      ligand (OEMol) - the ligand
      cores (int) - # of cores to run
      ligand_off_filename (string) - leap library for ligand
      ligand_frcmod_filename (string) - additional ligand parameters
      common_substructure (OEMol) - 
    """

    print "\nPREPARING SOLVATED COMPLEX"

    # create path if it doesn't exist
    if not os.path.exists(complex_path):
        os.makedirs(complex_path)

    # get ligand formal charge
    ligand_charge = formalCharge(ligand)

    # set up the system in tLEaP
    print "Solvating the complex with tleap..."
    system_prmtop_filename = os.path.join(complex_path,'system.prmtop')
    system_crd_filename = os.path.join(complex_path,'system.crd')
    tleap_input_filename = os.path.join(complex_path, 'setup-system.leap.in')
    tleap_output_filename = os.path.join(complex_path, 'setup-system.leap.out')
    clearance = 10.0 # clearance in A to add around protein
    contents = """
source leaprc.ff99SB
source leaprc.gaff

# load antechamber-generated additional parameters
mods = loadAmberParams %(ligand_frcmod_filename)s

# Load protein.
loadOff %(protein_off_filename)s

# Load ligand.
loadOff %(ligand_off_filename)s

# Create system.
system = combine { protein ligand }
""" % vars()
    # add counterions for ligand
    if (ligand_charge != 0):
        nions = abs(ligand_charge)
        if ligand_charge < 0: iontype = 'Na+'
        if ligand_charge > 0: iontype = 'Cl-'
        contents += """
# Add counterions for ligand (will be annihilated with ligand).
addions system %(iontype)s %(nions)d
""" % vars()
    #
    # add counterions for protein
    if (protein_charge != 0):
        nions = abs(protein_charge)
        if protein_charge < 0: iontype = 'Na+'
        if protein_charge > 0: iontype = 'Cl-'
        contents += """
# Add counterions for protein.
addions system %(iontype)s %(nions)d
""" % vars()
    #
    contents += """
# Solvate in truncated octahedral box.
# solvateOct system TIP3PBOX %(clearance)f
# solvate in rectilinear box
solvateBox system TIP3PBOX %(clearance)f

# Check the system
check system

# Write the system
saveamberparm system %(system_prmtop_filename)s %(system_crd_filename)s
""" % vars()
    write_file(tleap_input_filename, contents)
    command = 'tleap -f %(tleap_input_filename)s > %(tleap_output_filename)s' % vars()
    output = commands.getoutput(command)

    # convert to gromacs
    print "Converting to gromacs..."

    converter = os.path.join(os.getenv('MMTOOLSPATH'), 'converters', gromacs_converter)
    system_prefix = os.path.join(complex_path, 'system')
    current_path = os.getcwd()
    os.chdir(complex_path)
    if gromacs_converter == 'amb2gmx.pl':
        command = '%(converter)s --parmtop %(system_prmtop_filename)s --crd %(system_crd_filename)s --outname %(system_prefix)s' % vars()
    if gromacs_converter == 'acpype.py':    
        command = 'python2.6 %(converter)s -r -p system.prmtop -x system.crd -b system -o gmx' % vars()
    print command
    output = commands.getoutput(command)
    print output
    os.chdir(current_path)

    # set up perturbation
    print "Modifying topology file for perturbation..."
    system_top_filename = os.path.join(complex_path, 'system_GMX.top')
    perturbGromacsTopologyToIntermediate(system_top_filename, 'MOL', ligand, common_substructure, perturb_torsions = True)

    # set up mdp files
    print "Writing mdp files..."
    writeMdpFiles(complex_path, runtype, cores,'mutate')

    # no restraints currently for relative calculations.
    # Copy restraint scripts to run directory
    #shutil.copy(os.path.join(script_basepath,'compute_angles.py'),complex_path)
    #shutil.copy(os.path.join(script_basepath,'restraint_topology.py'),complex_path)
    # write batch queue script
    #print "Writing batch queue script..."
    #complex_path = os.path.abspath(complex_path)
    #contents = """\

    # write batch scripts
    print "Writing batch scripts..."
    writeBatchScript(complex_path,jobname,runtype,cores,restraints=False)

    return
#=============================================================================================
def setup_system(protein_pdb_filename, ligand_path, ligand, common_substructure, work_path, jobname, runtype, phases):
    """Set up a system for alchemical free energy calculation in gromacs.

    ARGUMENTS
      protein_pdb_filename (string) - name of ffamber-named protein PDB file
      ligand_path (string) - pathname containing ligand files (.off, .frcmod, .mol2)
      ligand (OEMol) - ligand
      common_substructure (OEMol) - commmon substructure
      work_path (string) - name of directory to place files in
      jobname (string) - job name to prepend for batch queue system
      phases (string) - phases to set up -- any set of 'vacuum', 'solvent', and 'complex'

    """

    # get ligand name
    ligand_name = ligand.GetTitle()

    # get current directory
    current_path = os.getcwd()

    # create the work directory if it doesn't exist
    if not os.path.exists(work_path):
        os.makedirs(work_path)

    # SET UP PROTEIN TOPOLOGY

    print "\nCONSTRUCTING PROTEIN TOPOLOGY"

    # convert protein PDB file to AMBER naming conventions, dropping atoms that AMBER won't recognize (like protons)
    print "Converting PDB naming to AMBER for %s..." % protein_pdb_filename
    protein_amberpdb_filename = os.path.join(work_path, 'protein-ambernames.pdb')
    rename_pdb_for_amber(protein_pdb_filename, protein_amberpdb_filename)
    if not os.path.exists(protein_amberpdb_filename):
        print "AMBER naming conversion failed."
        return

    # run leap to set up protein and report on net charge
    print "Running LEaP to set up protein..."

    protein_prmtop_filename = os.path.join(work_path,'protein.prmtop')
    protein_crd_filename = os.path.join(work_path,'protein.crd')
    protein_off_filename = os.path.join(work_path, 'protein.off')

    tleap_input_filename = os.path.join(work_path, 'setup-protein.leap.in')
    tleap_output_filename = os.path.join(work_path, 'setup-protein.leap.out')
    contents = """
# Load AMBER ff99SB parameters.
source leaprc.ff99SB

# Load PDB file with AMBER naming conventions, stripped of hydrogens.
protein = loadpdb %(protein_amberpdb_filename)s

# check the system
check protein

# report net charge
charge protein

# write out parameters
saveAmberParm protein %(protein_prmtop_filename)s %(protein_crd_filename)s

# write as LEaP object
saveOff protein %(protein_off_filename)s

# exit
quit
""" % vars()
    write_file(tleap_input_filename, contents)
    command = 'tleap -f %(tleap_input_filename)s > %(tleap_output_filename)s' % vars()
    output = commands.getoutput(command)

    # extract total charge
    protein_charge = commands.getoutput('grep "Total unperturbed charge" %(tleap_output_filename)s | cut -c 27-' % vars())
    protein_charge = int(round(float(protein_charge))) # round to nearest whole charge
    print protein_charge

    # READ LIGAND INFORMATION

    # set ligand path names
    ligand_off_filename = os.path.join(ligand_path, '%(ligand_name)s.off' % vars())
    if not os.path.exists(ligand_off_filename):
        print "%(ligand_off_filename)s not found, skipping..." % vars()
        return
    ligand_frcmod_filename = os.path.join(ligand_path, '%(ligand_name)s.frcmod' % vars())
    if not os.path.exists(ligand_frcmod_filename):
        print "%(ligand_frcmod_filename)s not found, skipping..." % vars()
        return
    
    # SET UP SIMULATIONS
    if 'solvent' in phases:
        # PREPARE SOLVATED LIGAND
        solvent_path = os.path.join(work_path, 'solvent')
        setup_solvent_simulation(solvent_path, jobname, runtype, ligand, cores, ligand_off_filename, ligand_frcmod_filename, common_substructure)
        
    if 'complex' in phases:
        # PREPARE SOLVATED COMPLEX
        complex_path = os.path.join(work_path, 'complex')
        setup_complex_simulation(complex_path, jobname, runtype, ligand, cores, ligand_off_filename, ligand_frcmod_filename, protein_off_filename, protein_charge, common_substructure)

    if 'vacuum' in phases:
        # PREPARE vacuum simulation
        # !!!!not yet set up!!!!
        vacuum_path = os.path.join(work_path, 'vacuum')
        setup_vacuum_simulation(vacuum_path, jobname, runtype, ligand, cores, ligand_off_filename, ligand_frcmod_filename, common_substructure)

    return

#=============================================================================================
# MAIN
#=============================================================================================

# signal error if some files not found
if not os.path.exists(protein_pdb_filename):
    raise "%s not found" % protein_pdb_filename

# create output directory if it doesn't yet exist
if not os.path.exists(work_basepath):
    os.makedirs(work_basepath)

# read all ligands
ligands = list()
for ligand_name in ligand_name_list:
    # attempt to load the ligand with GAFF parameters
    ligand = loadGAFFMolecule(ligand_basepath, ligand_name)

    # append to list
    ligands.append(ligand)
if len(ligands) == 0:
    raise "No ligand structures found.\nligand_basepath = %s\nligand_name_list = %s\n" % (ligand_basepath, repr(ligand_name_List))
print "%d ligands loaded" % len(ligands)

# find common substructure
print "Determining common substructure using MCSS (this may take a couple minutes)..."
common_substructure = determineCommonSubstructure(ligands)

# DEBUG
print "\n%s" % common_substructure.GetTitle()
for atom in common_substructure.GetAtoms():
    atomname = atom.GetName()
    print "%6s : %16s %6.3f" % (atom.GetName(), atom.GetType(), atom.GetPartialCharge())

# find RMS-fit charges for common intermediate
print "Determining RMS-fit charges for common intermediate..."
common_substructure = determineMinimumRMSCharges(common_substructure, ligands)

# DEBUG
print "\n%s" % common_substructure.GetTitle()
for atom in common_substructure.GetAtoms():
    atomname = atom.GetName()
    print "%6s : %16s %6.3f" % (atom.GetName(), atom.GetType(), atom.GetPartialCharge())

# Write out common substructure in GAFF format.
intermediate_filename = os.path.join(work_basepath, 'common-substructure.gaff.mol2')
print "Writing common intermediate with GAFF atomtypes to %(intermediate_filename)s" % vars()
writeMolecule(common_substructure, intermediate_filename, preserve_atomtypes = True)
# Write out common substructure in Tripos format.
intermediate_filename = os.path.join(work_basepath, 'common-substructure.tripos.mol2')
print "Writing common intermediate with Tripos atomtypes to %(intermediate_filename)s" % vars()
writeMolecule(common_substructure, intermediate_filename)

# DEBUG
print "\n%s" % common_substructure.GetTitle()
for atom in common_substructure.GetAtoms():
    atomname = atom.GetName()
    print "%6s : %16s %6.3f" % (atom.GetName(), atom.GetType(), atom.GetPartialCharge())

# clean out leap.log
commands.getoutput('rm -f leap.log')

# copy job submission script into newly-created workign directory
print commands.getoutput('cp %(script_basepath)s/submit-lsf.py %(work_basepath)s/submit-lsf.py' % vars())

# keep track of all work directories set up for submit-lsf.py script
directories = list()
directories_file = '%(work_basepath)s/directories' % vars()

# iterate over ligands and set up each complex
for ligand in ligands:
    # construct pathnames
    ligand_name = ligand.GetTitle()
    work_path = '%(work_basepath)s/%(ligand_name)s' % vars()
    directories.append(work_path)

    # setup system
    jobname = '%(ligand_name)s' % vars()
    setup_system(protein_pdb_filename, ligand_basepath, ligand, common_substructure, work_path, jobname, runtype, phases)

    # update directories file
    write_file(directories_file, directories)




