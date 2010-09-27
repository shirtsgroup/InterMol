# Check ligand energies
#
# Written by Michael Zhu <mnz3v@virginia.edu> 2009-10-26
#=============================================================================================

#=============================================================================================
# IMPORTS
#=============================================================================================

from mmtools.moltools.relativefeptools import *
from mmtools.gromacstools.GromacsData import *
from mmtools.gromacstools.MdpFile import *
from math import sqrt
from numpy import *
import commands, os, re

#=============================================================================================
# PARAMETERS
#=============================================================================================

gromacs_directory = '/usr/local/gromacs/bin'
working_directory = '/home/mnz3v/mmtools/relativetools/Examples/FKBP'

# Target ligand group.
#ligand_name_list = [ ('jnk.aff-%d' % index) for index in (8, 10, 11, 15, 22, 26, 32, 52) ]
ligand_name_list = [ ('jnk.aff-%d' % index) for index in (13, 16, 53) ]

# Path to parameterized ligands.
ligand_basepath = 'ligands-parameterized'

# complex/solvent
ligand_type = 'solvent'

# Desired values to input to g_energy
xvg_values = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

#=============================================================================================
# SUBROUTINES
#=============================================================================================

def readFile(filename):
    """Write file into array of lines.

    ARGUMENTS
        filename (string) - the name of file to be read

    RETURNS
        lines (list of strings) - the lines of each file
    """

    lines = []
    f = open(filename, 'r')

    for line in f:
        line = line.strip()
        lines.append(line)
    f.close()

    return lines

def determineMatchArray(first_ligand_name, second_ligand_name, min_atoms = 4):
    """Find a list of the common substructure elements shared by both ligands.

    The atom type name strings and integer bond types are used to obtain an exact match.

    ARGUMENTS
        first_ligand_name (string) - ligand to be compared
        second_ligand_name (string) - ligand to be compared

    OPTIONAL ARGUMENTS
        min_atoms (int) - minimum number of atoms for substructure match (default: 4)

    RETURN VALUES
        match_list (list) - list of corresponding matching atoms
    """

    # First initialize ligands.
    first_ligand = OEMol(loadGAFFMolecule('%s/%s'%(working_directory,ligand_basepath), first_ligand_name))
    second_ligand = OEMol(loadGAFFMolecule('%s/%s'%(working_directory,ligand_basepath), second_ligand_name))

    # Set second ligand to be searched
    mcss = OEMCSSearch(second_ligand, OEExprOpts_StringType, OEExprOpts_IntType)

    # Set minimum amount of matched atoms
    mcss.SetMinAtoms(min_atoms)
    for match in mcss.Match(first_ligand):
        matched_atoms = list()
        for matchpair in match.GetAtoms():
            matched_atoms.append([matchpair.target.GetIdx() + 1, matchpair.pattern.GetIdx() + 1])

    return matched_atoms

def searchArray(match_list, value):
    """Searches match list for atom number and outputs matching atom number.

    ARGUMENTS
        match_list (list of integers) - the list of atoms numbers and the matching atom numbers
        value (integer) - the number of the atom to be matched

    RETURNS
        match (integer/boolean) - the atom number corresponding to value, False if there is no match
    """

    match = False

    for i in range(len(match_list)):
        if match_list[i][1] == value:
            match = match_list[i][0]

    return match

def findNumberOfAtoms(gromacs_structure):
    """Finds the number of atoms in the ligand for a gromacs structure.

    ARGUMENTS
        gromacs_structure (gromacs structure object) - the object containing atom information

    RETURNS
        length (integer) - the number of atoms in the ligand
    """

    length = 1

    while gromacs_structure.atoms[length].resname == gromacs_structure.atoms[length + 1].resname:
        length = length + 1

    return length

def findDistance(atomi, atomj):
    """Finds the cartesian distance between two atoms.

    ARGUMENTS
        atomi (list of double) - the list containing the coordinates of atom i
        atomj (list of double) - the list containing the coordinates of atom j

    RETURNS
        distance - the cartesian distance between atom i and atom j
    """

    distance = math.sqrt((atomi.x - atomj.x)**2 + (atomi.y - atomj.y)**2 + (atomi.z - atomj.z)**2)

    return distance

def calculateAtomPosition(i, position, moved_atoms, initial_GroStructure, target_GroStructure):
    """Calculates the new position of the atom in the mutated ligand based on the closest previously moved atom.

    ARGUMENTS
        i (integer) - the number of the atom being moved
        moved_atoms (list) - a list of atoms and their coordinates that have already been moved 
        initial_GroStructure (gromacs structure) - the gromacs structure containing the positions of the atoms in the correct spatial arrangement
        target_GroStructure (gromacs structure) - the gromacs structure containing the atoms to be moved

    RETURNS
        position (array) - the array specifying the new positions of the atoms to be moved
    """

    difference = array([])

    # Creates an array of differences between the moved atoms and the atom to be moved
    for j in moved_atoms:
        difference = append(difference, findDistance(target_GroStructure.atoms[i], target_GroStructure.atoms[j]))

    # Finds the closest moved atom to the atom to be moved and finds the new coordinates for the atom to be moved
    if difference.min() < 0.2:
        closest_atom = moved_atoms[difference.argmin()]
        position[i, 0] = target_GroStructure.atoms[i].x - target_GroStructure.atoms[closest_atom].x + position[closest_atom, 0]
        position[i, 1] = target_GroStructure.atoms[i].y - target_GroStructure.atoms[closest_atom].y + position[closest_atom, 1]
        position[i, 2] = target_GroStructure.atoms[i].z - target_GroStructure.atoms[closest_atom].z + position[closest_atom, 2]
        position[i, 3] = target_GroStructure.atoms[i].vx - target_GroStructure.atoms[closest_atom].vx + position[closest_atom, 3]
        position[i, 4] = target_GroStructure.atoms[i].vy - target_GroStructure.atoms[closest_atom].vy + position[closest_atom, 4]
        position[i, 5] = target_GroStructure.atoms[i].vz - target_GroStructure.atoms[closest_atom].vz + position[closest_atom, 5]

    return position

def calculatePositionArray(initial_GroStructure, target_GroStructure, location_array):
    """Calculates an array of new positions for the mutating atoms to be moved.

    ARGUMENTS
        initial_GroStructure (gromacs sturcture) - the gromacs structure containing the positions of the atoms in the correct spatial arrangement
        target_GroStructure (gromacs structure) - the gormacs structure containing the atoms to be moved
        location_array (array) - the array specifying where mutating atoms are to be moved

    RETURNS
        position (array) - the array containing the new positions of the mutating atoms
    """

    length = findNumberOfAtoms(target_GroStructure)
    unmoved_atoms = range(length + 1)
    moved_atoms = []
    position = zeros((length + 1, 6))

    # Finds the position of the atoms specified by the location array
    for i in range(length + 1):
        n = searchArray(location_array, i + 1)
        if n <> False:
            position[i, 0] = initial_GroStructure.atoms[n - 1].x
            position[i, 1] = initial_GroStructure.atoms[n - 1].y
            position[i, 2] = initial_GroStructure.atoms[n - 1].z
            position[i, 3] = initial_GroStructure.atoms[n - 1].vx
            position[i, 4] = initial_GroStructure.atoms[n - 1].vy
            position[i, 5] = initial_GroStructure.atoms[n - 1].vz
            moved_atoms.append(i)
            unmoved_atoms.remove(i)

    # Takes the atoms in the unmoved atoms list and finds the closest atom using calculateAtomPosition
    while unmoved_atoms <> []:
        for i in unmoved_atoms:
            position = calculateAtomPosition(i, position, moved_atoms, initial_GroStructure, target_GroStructure)
            if position[i, 0] <> 0:
                moved_atoms.append(i)
                unmoved_atoms.remove(i)

    return position

def mutateGroFile(initial_GroFile, target_GroFile, location_array):
    """Using the location array, the atoms in the target gromacs file are given new coordinates matching those in the initial gromacs file.

    ARGUMENTS
        initial_GroFile (string) - the name of the gromacs file containing the correct spatial arrangement
        target_GroFile (string) - the name of the gromacs file containing the atoms to be moved
        location_array (array) - the array containing the specificatiions of where to move which atom

    RETURNS
        mutated_GroStructure (gromacs structure) - the gromacs structure containing all of the information of the mutated gromacs file
    """

    initial_GroStructure = GromacsStructure()
    initial_GroStructure.load(initial_GroFile)

    target_GroStructure = GromacsStructure()
    target_GroStructure.load(target_GroFile)

    mutated_GroStructure = GromacsStructure()
    mutated_atoms_list = GromacsAtoms()

    length = findNumberOfAtoms(target_GroStructure)

    # Finds the new position for the mutating atoms
    position = calculatePositionArray(initial_GroStructure, target_GroStructure, location_array)

    # Finds the parameters for the mutated ligand and stores all parameters into mutatedGroStructure
    for i in range(length + 1):
        atom = GromacsAtom()
        n = searchArray(location_array, i + 1)
        if n == False:
            atom.resnum = target_GroStructure.atoms[i].resnum
            atom.resname = target_GroStructure.atoms[i].resname
            atom.atomname = target_GroStructure.atoms[i].atomname
            atom.atomnum = i + 1
        else:
            atom.resnum = initial_GroStructure.atoms[n - 1].resnum
            atom.resname = initial_GroStructure.atoms[n - 1].resname
            atom.atomname = initial_GroStructure.atoms[n - 1].atomname
            atom.atomnum = i + 1
        atom.x = position[i, 0]
        atom.y = position[i, 1]
        atom.z = position[i, 2]
        atom.vx = position[i, 3]
        atom.vy = position[i, 4]
        atom.vz = position[i, 5]
        mutated_atoms_list.append(atom)

    initial_GroStructure.trim(0, findNumberOfAtoms(initial_GroStructure))
    mutated_atoms_list.extend(initial_GroStructure.atoms)
    mutated_GroStructure.name = target_GroStructure.name
    mutated_GroStructure.header = target_GroStructure.header
    mutated_GroStructure.appendatoms(mutated_atoms_list)
    mutated_GroStructure.boxstring = initial_GroStructure.boxstring

    return mutated_GroStructure

def runGmx(cmd):
    """Runs a command from the command line.

    ARGUMENTS
        cmd (string) - the command to be executed

    """

    # Executes the command
    output = commands.getoutput(cmd)

    # If there is a 'Fatal error', shows which command gave the error
    if output.count('Fatal error') > 0:
        print '*** There was a fatal error in the following command: ***'
        print '*********************************************************'
        print cmd
        print '*********************************************************'
        print 'Exiting...'
        sys.exit(1)

    pass

def findEnergy(gromacs_directory, mdp_file, gro_file, top_file, ndx_file, xvg_values):
    """Runs gromacs energy calculations given an mdp, gro, and top file.

    ARGUMENTS
        gromacs_directory (string) - the location of the gromacs directory
        mdp_file (string) - the name of the mdp file
        gro_file (string) - the name of the gro file
        top_file (string) - the name of the top file
        xvg_values (list of integers) - the values to be intputted to g_energy

    RETURNS
        energy (array) - array containing xvg values
    """

    # Executes grompp_d
    grompp_d = '%s/grompp_d -f %s -c %s -p %s -o minimize.tpr -maxwarn 10000 -n %s'%(gromacs_directory, mdp_file, gro_file, top_file, ndx_file)
    runGmx(grompp_d)

    # Executes mdrun_d
    mdrun_d = '%s/mdrun_d -s minimize.tpr -x minimize.xtc -c %s -e minimize.edr -g minimize.log'%(gromacs_directory, gro_file)
    runGmx(mdrun_d)

    # Executes g_energy_d
    xvg = ''
    for i in range(len(xvg_values)):
        xvg = '%s %s'%(xvg, str(xvg_values[i]))
    g_energy_d = '/bin/echo%s 0 | %s/g_energy_d -f minimize.edr -o energy.xvg'%(xvg, gromacs_directory)
    runGmx(g_energy_d)

    # Places xvg values into a n array
    f = open('energy.xvg', 'r')
    lastline = []

    for line in f.readlines()[-1:]:
        lastline.append(line)
    f.close()
    line = line.strip()
    energy = array(re.split(' *', line)).astype(float)

    return energy

def modifySolventInTopFile(initial_top_file, modified_top_file, new_top_name):
    """Modifies the number of solvent molecules in the topology file so that it matches with the mutated atom

    ARGUMENTS
        initial_top_file (string) - the topology file whose amount of solvent molecules is to be obtained
        modified_top_file (string) - the topology file whose amount of solvent molecules is to be replaced

    """

    initial_top = open(initial_top_file, 'r')
    lastline = []

    # Obtains the number of solvent molecules and stores 
    for line in initial_top.readlines()[-1:]:
        lastline.append(line)
    initial_top.close()

    # Stores all values of the modified topology file
    modified_top = open(modified_top_file, 'r')
    modified_top_lines = modified_top.readlines()
    modified_top.close()

    # Replaces the solvent number from modified topology with the number from initial topology file
    modified_top_lines[len(modified_top_lines) - 1] = lastline[0]

    # Writes the new topology file as newTopName
    new_top = open(new_top_name, 'w')
    new_top.writelines(modified_top_lines)
    new_top.close()

    pass

def modifyTop(top_file, new_system_nmols, new_water_nmols):

    top = open(top_file, 'r')
    top_lines = top.readlines()
    top.close()

    lines = top_lines[:]
    i = 0
    nmols = []
    while i < len(lines):
        if len(lines[i]) >= len('[ molecules ]'):
            if lines[i][0:len('[ molecules ]')] == '[ molecules ]':
                lines.pop(i)
                while lines[i][0] == ';':
                    lines.pop(i)
                fields = lines.pop(i).split()
                while len(fields) >= 2:
                    nmols.append([fields[0], fields[1]])
                    if i < len(lines):
                        fields = lines.pop(i).split()
                    else:
                        fields = []
        i = i + 1

    lines = top_lines[:]
    i = 0
    systemtitle = ''
    while i < len(lines):
        if len(lines[i]) >= len('[ system ]'):
            if lines[i][0:len('[ system ]')] == '[ system ]':
                lines.pop(i)
                while lines[i][0] == ';':
                    lines.pop(i)
                systemtitle = lines.pop(i)
            else:
                i=i+1
        else:
            i=i+1

    for i in range(len(top_lines)):
        if top_lines[i].endswith('%s\n'%(nmols[1][1])) == True:
            top_lines[i] = top_lines[i].replace(nmols[1][1], str(new_water_nmols))
        if top_lines[i].startswith(systemtitle.split()[0]) == True:
            top_lines[i] = top_lines[i].replace(systemtitle.split()[0], str(new_system_nmols))

    top = open(top_file, 'w')
    top.writelines(top_lines)
    top.close()

    return

def removeIndexFromNdx(ndx_file, max_index):
    f = open(ndx_file, 'r')
    lines = f.readlines()
    f.close()
    final = []
    iset = 0

    for line in lines:
        if line == '\n':
            continue
        elif line.strip()[0] == '[':
            setname = '%s%s' %('final',iset)
            final.append(setname)
            final[iset] = []
            final[iset].append(line)
            iset += 1
        else:
            elements = line.split(' ')
            if elements[-1] == '\n':
                elements.pop(-1)
            final[iset-1].append(elements)

    top = []
    for i in range(len(final)):
        top.append(final[i][0])
        for j in range(1,len(final[i])):
            top.append(final[i][j])
            if int(final[i][j][0]) > max_index:
                top.pop()
            for k in range(len(final[i][j])-1):
                if int(final[i][j][k]) == max_index+1:
                    top.pop(j)
                    top.append(final[i][j][0:k])

    finaltop = []
    for i in range(len(top)):
        if top[i][0].strip()[0] != '[':
            top[i].append('\n')
            line = ' '.join(top[i])
        else:
            line = top[i]
        finaltop.append(line)
    f = open(ndx_file,'w')
    f.writelines(finaltop)
    f.close()

    return

def findLigandEnergy(working_directory, gromacs_directory, ligand_name_list, ligand_basepath, ligand_type, xvg_values, state):
    """Runs findEnergy on all of the ligands

    ARGUMENTS
        working_directory (string) - the current working directory
        gromacs_directory (string) - the directory containing gromacs
        ligand_name_list (array of ) - the list of ligands to be analyzed
        ligand_basepath (string) - the directory containing the files to be analyzed
        ligand_type (array of strings) - the type(s) of ligands to be analyzed
        xvg_values (array of integers) - the array containing the desired xvg_values
        state (string) - initial/mutatedInitial/mutatedFinal

    RETURNS
        ligand (array) - returns the array containing all of the energy terms for all of the comparisons for a single state
    """

    os.chdir(working_directory)
    n = len(xvg_values) + 1
    ligand = array([])

    # Initial state calculation
    if state == 'initial':

        print 'Minimizing Initial'

        # Runs findEnergy and appends values to array
        for i in range(len(ligand_name_list)):

            print 'Minimizing %s'%(ligand_name_list[i])

            new_working_directory = '%s/fep/%s/%s'%(working_directory, ligand_name_list[i], ligand_type)
            commands.getoutput('mkdir %s/temp'%(new_working_directory))
            commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'minimize.mdp', new_working_directory))
            #commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'equilibration.mdp', new_working_directory))
            #commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.g96', new_working_directory))
            commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.gro', new_working_directory))
            commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.ndx', new_working_directory))
            commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.top', new_working_directory))
            os.chdir('%s/temp'%(new_working_directory))

            # Sets initial state parameters
            mdp = MdpFile('minimize.mdp')
            mdp.setParameter('free-energy', 'no')
            # Settings for Gromacs, not Gromacs_dg
            mdp.setParameter('rlist', mdp.getParameter('rcoulomb'))
            mdp.setParameter('constraint-algorithm', 'lincs')
            # End of Gromacs settings
            mdp.write('minimize.mdp')

            line = findEnergy(gromacs_directory, 'minimize.mdp', 'system.gro', 'system.top', 'system.ndx', xvg_values)
            ligand = append(ligand.reshape(len(ligand), n), line.reshape(1, n), axis=0)
            os.chdir(new_working_directory)
            os.system('rm -rf temp')

    # Mutated Initial calculation
    elif state == 'mutatedInitial':

        print 'Minimizing Initial Mutation'

        # Runs findEnergy and appends values to array
        for i in range(len(ligand_name_list)):

            print 'Minimizing %s'%(ligand_name_list[i])

            new_working_directory = '%s/fep/%s/%s'%(working_directory, ligand_name_list[i], ligand_type)
            commands.getoutput('mkdir %s/temp'%(new_working_directory))
            commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'minimize.mdp', new_working_directory))
            #commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'equilibration.mdp', new_working_directory))
            #commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.g96', new_working_directory))
            commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.gro', new_working_directory))
            commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.ndx', new_working_directory))
            commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.top', new_working_directory))
            os.chdir('%s/temp'%(new_working_directory))

            # Sets mutated initial state parameters
            mdp = MdpFile('minimize.mdp')
            mdp.setParameter('free-energy', 'yes')
            mdp.setParameter('init-lambda', '0')
            # Settings for Gromacs, not Gromacs_dg
            mdp.setParameter('rlist', mdp.getParameter('rcoulomb'))
            mdp.setParameter('constraint-algorithm', 'lincs')
            mdp.setParameter('coulombtype','PME')
            mdp.setParameter('ewald-rtol','1e-07')
            # End of Gromacs settings
            mdp.write('minimize.mdp')

            line = findEnergy(gromacs_directory, 'minimize.mdp', 'system.gro', 'system.top', 'system.ndx', xvg_values)
            ligand = append(ligand.reshape(len(ligand), n), line.reshape(1, n), axis=0)
            os.chdir(new_working_directory)
            os.system('rm -rf temp')

    # Mutated Final calculation
    elif state == 'mutatedFinal':

        print 'Minimizing Final Mutation'

        for i in range(len(ligand_name_list)):

            new_working_directory = '%s/fep/%s/%s'%(working_directory, ligand_name_list[i], ligand_type)
            commands.getoutput('mkdir %s/temp'%(new_working_directory))
            commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'minimize.mdp', new_working_directory))
            #commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'equilibration.mdp', new_working_directory))
            #commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.g96', new_working_directory))
            commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.gro', new_working_directory))
            commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.ndx', new_working_directory))
            commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.top', new_working_directory))
            os.chdir('%s/temp'%(new_working_directory))

            # Sets mutated final state parameters
            mdp = MdpFile('minimize.mdp')
            mdp.setParameter('free-energy', 'yes')
            mdp.setParameter('init-lambda', '1')
            # Settings for Gromacs, not Gromacs_dg
            mdp.setParameter('rlist', mdp.getParameter('rcoulomb'))
            mdp.setParameter('constraint-algorithm', 'lincs')
            mdp.setParameter('coulombtype','PME')
            mdp.setParameter('ewald-rtol','1e-07')
            # End of Gromacs settings
            mdp.write('minimize.mdp')

            # Runs findEnergy and appends values to array
            line = findEnergy(gromacs_directory, 'minimize.mdp', 'system.gro', 'system.top', 'system.ndx', xvg_values)
            ligand = append(ligand.reshape(len(ligand), n), line.reshape(1, n ), axis=0)
            os.chdir(new_working_directory)
            os.system('rm -rf temp')

            for j in range(i+1, len(ligand_name_list)):

                print 'Minimizing %s  %s'%(ligand_name_list[i], ligand_name_list[j])

                new_working_directory = '%s/fep/%s/%s'%(working_directory, ligand_name_list[j], ligand_type)
                commands.getoutput('mkdir %s/temp'%(new_working_directory))
                commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'minimize.mdp', new_working_directory))
                #commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'equilibration.mdp', new_working_directory))
                commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.gro', new_working_directory))
                #commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.gro', new_working_directory))
                commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.ndx', new_working_directory))
                commands.getoutput('cp %s/%s %s/temp'%(new_working_directory, 'system.top', new_working_directory))
                os.chdir('%s/temp'%(new_working_directory))

                # Sets mutated final state parameters
                mdp = MdpFile('minimize.mdp')
                mdp.setParameter('free-energy', 'yes')
                mdp.setParameter('init-lambda', '1')
                # Settings for Gromacs, not Gromacs_dg
                mdp.setParameter('rlist', mdp.getParameter('rcoulomb'))
                mdp.setParameter('constraint-algorithm', 'lincs')
                mdp.setParameter('coulombtype','PME')
                mdp.setParameter('ewald-rtol','1e-07')
                # End of Gromacs settings
                mdp.write('minimize.mdp')

                matchArray = determineMatchArray(ligand_name_list[i], ligand_name_list[j])
                groPair1 = '%s/fep/%s/%s/%s'%(working_directory, ligand_name_list[i], ligand_type, 'system.gro')
                groPair2 = '%s/fep/%s/%s/%s'%(working_directory, ligand_name_list[j], ligand_type, 'system.gro')
                mutatedGroStructure = GromacsStructure()
                mutatedGroStructure = mutateGroFile(groPair1, groPair2, matchArray)
                mutatedGroStructure.write('system_mut')
                top1 = '%s/fep/%s/%s/%s'%(working_directory, ligand_name_list[i], ligand_type, 'system.top')
                top2 = '%s/fep/%s/%s/%s'%(working_directory, ligand_name_list[j], ligand_type, 'system.top')
                modifySolventInTopFile(top1, top2, 'system_mut.top')
                new_system_nmols = mutatedGroStructure.natoms
                new_water_nmols = (new_system_nmols-findNumberOfAtoms(mutatedGroStructure))/3
                modifyTop('system_mut.top', new_system_nmols, new_water_nmols)
                commands.getoutput('cp system.ndx system_mut.ndx')
                removeIndexFromNdx('system_mut.ndx', new_system_nmols)

                # Runs findEnergy and appends values to array
                line = findEnergy(gromacs_directory, 'minimize.mdp', 'system_mut.gro', 'system_mut.top', 'system_mut.ndx', xvg_values)
                ligand = append(ligand.reshape(len(ligand), n), line.reshape(1, n ), axis=0)
                os.chdir(new_working_directory)
                os.system('rm -rf temp')

    return ligand

def createLegend(xvg_values):
    """Creates a legend based on the chosen xvg values

    ARGUMENTS
        xvg_values (array) - the array containing the desired xvg_values

    RETURNS
        legend (array) - the array containing the desired legend
    """

    g_energy_legend = ['Time', 'Bond', 'Angle', 'Proper-Dih.', 'Ryckaert-Bell.', 'LG-14', 'Coulomb-14', 'LJ-(SR)', 'Coulomb-(SR)', 'Potential', 'Kinetic-En.', 'Total-Energy', 'Temperature', 'Pressure-(bar)', 'Vir-XX', 'Vir-XY', 'Vir-XZ', 'Vir-YX', 'Vir-YY', 'Vir-YZ', 'Vir-ZX', 'Vir-ZY', 'Vir-ZZ', 'Pres-XX-(bar)', 'Pres-XY-(bar)', 'Pres-XZ-(bar)', 'Pres-YX-(bar)', 'Pres-YY-(bar)', 'Pres-YZ-(bar)', 'Pres-ZX-(bar)', 'Pres-ZY-(bar)', 'Pres-ZZ-(bar)', '#Surf*SurfTen', 'Mu-X', 'Mu-Y', 'Mu-Z', 'T-rest']

    legend = ['Time']

    for i in xvg_values:
        legend.append(g_energy_legend[i])

    return legend

def compareInitialLigands(initial_ligand, initial_mutated_ligand, legend, ligand_name_list):
    """Prints out the ligand comparisons between the initial ligand and the initial mutated ligand that do not give the same energy terms

    ARGUMENTS
        initial_ligand (array) - the energies of the initial ligand with no free energy turned on
        initial_mutated_ligand (array) - the energies of the initial mutated ligand with free energy turned on but with no lambda values
        legend (array of strings) - the legend of the values to be compared
        ligand_name_list (array of strings) - the array containing the names of the ligands to be compared
    """

    length = len(ligand_name_list)
    unequal_atoms = []
    unequal_atoms_list = []
    extended_ligand_list = []
    extended_initial_ligand = []
    difference = []

    if len(initial_ligand) != len(initial_mutated_ligand):
        for i in range(length):
            for j in range(i+1, length):
                extended_ligand_list.append(ligand_name_list[i])
                extended_ligand_list.append(ligand_name_list[j])
                extended_initial_ligand.append(initial_ligand[i])
                extended_initial_ligand.append(initial_ligand[j])
        difference = extended_initial_ligand-initial_mutated_ligand
    else:
        extended_ligand_list = ligand_name_list
        difference = initial_ligand - initial_mutated_ligand

    # Determines which atoms do not give the same energy
    for i in range(len(difference)):
        unequal_params = [i]
        for j in range(len(legend)):
            if difference[i][j] != 0:
                unequal_params.append(j)
        if len(unequal_params) > 1:
            unequal_atoms.append(unequal_params)
    for i in range(len(unequal_atoms)):
        unequal_atoms_list = append(unequal_atoms_list, unequal_atoms[i][0])

    # Prints out the terms that do not equal along with the corresponding differences
    print '\n\nDifferences When Comparing Energies With and Without Free Energy'
    for i in range(len(difference)):
        if i in unequal_atoms_list:
            n = int(where(unequal_atoms_list == i)[0])

            print '\n%s'%(extended_ligand_list[i])
            print 'Mismatched Params.',
            for j in range(len(unequal_atoms[n]) - 1):
                print '%16s'%(legend[unequal_atoms[n][j + 1]]),
            print '\nDifference        ',
            for j in range(len(unequal_atoms[n]) - 1):
                print '%16.6f'%(difference[unequal_atoms[n][0]][unequal_atoms[n][j + 1]]),
            print ''
        else:
            print '\n%s'%(extended_ligand_list[i])
            print 'No Difference'
    pass

def compareMutatedLigands(mutated_ligand, legend, ligand_name_list):
    """Prints out the ligand comparisons between the two mutated forms that do not give the same energy terms

    ARGUMENTS
        mutated_ligand (array) - the array containing the energies of the mutated ligands
        legend (array of strings) - the legend of the values to be compared
        ligand_name_list (array of strings) - the array containing the names of the ligands to be compared
    """

    unequal_atoms = []
    unequal_atoms_list = []

    name_comparison = []
    difference = array([])
    n = len(ligand_name_list)
    index = 0

    # Creates comparison and difference arrays
    for i in range(n):
        for j in range(i+1, n):
            name_comparison.append('Comparing %s to %s'%(ligand_name_list[i], ligand_name_list[j]))
            line = mutated_ligand[i] - mutated_ligand[j]
            difference = append(difference.reshape(len(difference), len(legend)), line.reshape(1,len(legend)),axis=0)

    # Determines which atoms do not give the same energy
    for i in range(len(difference)):
        unequal_params = [i]
        for j in range(len(legend)):
            if difference[i][j] != 0:
                unequal_params.append(j)
        if len(unequal_params) > 1:
            unequal_atoms.append(unequal_params)
    for i in range(len(unequal_atoms)):
        unequal_atoms_list = append(unequal_atoms_list, unequal_atoms[i][0])

    # Prints out the terms that do not equal along with the corresponding differences
    print '\n\nDifferences When Comparing the Mutated Forms'
    for i in range(len(difference)):

        print '\n%s'%(name_comparison[i])

        if i in unequal_atoms_list:
            n = int(where(unequal_atoms_list == i)[0])

            print 'Mismatched Params.',
            for j in range(len(unequal_atoms[n]) - 1):
                print '%16s'%(legend[unequal_atoms[n][j + 1]]),
            print '\nDifference        ',
            for j in range(len(unequal_atoms[n]) - 1):
                print '%16.6f'%(difference[unequal_atoms[n][0]][unequal_atoms[n][j + 1]]),
            print ''
        else:
            print 'No Difference'
    pass

#=============================================================================================
# MAIN
#=============================================================================================

# Finds array containing ligand values
LG0 = findLigandEnergy(working_directory, gromacs_directory, ligand_name_list, ligand_basepath, ligand_type, xvg_values, 'initial')
LGM0 = findLigandEnergy(working_directory, gromacs_directory, ligand_name_list, ligand_basepath, ligand_type, xvg_values, 'mutatedInitial')
LGM1 = findLigandEnergy(working_directory, gromacs_directory, ligand_name_list, ligand_basepath, ligand_type, xvg_values, 'mutatedFinal')

# Checks to see which ligands do not give the same energy
legend = createLegend(xvg_values)
compareInitialLigands(LG0, LGM0, legend, ligand_name_list)
compareMutatedLigands(LGM1, legend, ligand_name_list)
