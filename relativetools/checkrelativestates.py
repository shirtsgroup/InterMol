#!/usr/bin/env python

# Check ligand energies

# Written by Michael Zhu <mnz3v@virginia.edu> 2009-10-26
#=============================================================================================

"""
    Requirements:   Python 2.4 or higher
                    Gromacs 4.0 or higher
"""

#=============================================================================================
# IMPORTS
#=============================================================================================

from mmtools.moltools.relativefeptools import *
from mmtools.gromacstools.GromacsData import *
from mmtools.gromacstools.MdpFile import *
from numpy import *
import commands, os, re, pdb

#=============================================================================================
# PARAMETERS
#=============================================================================================

# Desired values to input to g_energy
xvg_values = '1 2 3 4 5 6 7 8 9 10 11'

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

def readMci(mci_file):
    """Reads the mci file and returns the contents in an array.

    ARGUMENTS
        mci_file (string) - the name of file to be read

    RETURNS
        new_mci (list of integers) - the values of the mci file
    """

    new_mci = []

    # Reads in mciFile
    old_mci = readFile(mci_file)

    # Places values into a list
    for i in range(len(old_mci)):
        new_mci.append(old_mci[i].split())

    # Increases all values by 1 to match Gromacs numbering
    for i in range(len(new_mci)):
        for j in range(len(new_mci[0])):
            new_mci[i][j] = int(new_mci[i][j])+ 1

    return new_mci

def determineMatchArray(ligand_basepath, first_ligand_name, second_ligand_name, min_atoms = 4):
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
    first_ligand = OEMol(loadGAFFMolecule(ligand_basepath, first_ligand_name))
    second_ligand = OEMol(loadGAFFMolecule(ligand_basepath, second_ligand_name))

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

    distance = sqrt((atomi.x - atomj.x)**2 + (atomi.y - atomj.y)**2 + (atomi.z - atomj.z)**2)

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
        if n != False:
            position[i, 0] = initial_GroStructure.atoms[n - 1].x
            position[i, 1] = initial_GroStructure.atoms[n - 1].y
            position[i, 2] = initial_GroStructure.atoms[n - 1].z
            position[i, 3] = initial_GroStructure.atoms[n - 1].vx
            position[i, 4] = initial_GroStructure.atoms[n - 1].vy
            position[i, 5] = initial_GroStructure.atoms[n - 1].vz
            moved_atoms.append(i)
            unmoved_atoms.remove(i)

    # Takes the atoms in the unmoved atoms list and finds the closest atom using calculateAtomPosition
    while unmoved_atoms != []:
        for i in unmoved_atoms:
            position = calculateAtomPosition(i, position, moved_atoms, initial_GroStructure, target_GroStructure)
            if position[i, 0] != 0:
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
        print os.getcwd()
        print '*********************************************************'
        print 'Exiting...'
        sys.exit(1)

    pass

def findEnergy(grompp, mdrun, g_energy, mdp_file, gro_file, top_file, xvg_values):
    """Runs gromacs energy calculations given an mdp, gro, and top file.

    ARGUMENTS
        gromacs_directory (string) - the location of the gromacs directory
        mdp_file (string) - the name of the mdp file
        gro_file (string) - the name of the gro file
        top_file (string) - the name of the top file
        xvg_values (string) - the values to be intputted to g_energy

    RETURNS
        energy (array) - array containing xvg values
    """

    # Executes grompp_d
    cmd = '%s -f %s -c %s -p %s -o minimize.tpr -maxwarn 10000'%(grompp, mdp_file, gro_file, top_file)
    runGmx(cmd)

    # Executes mdrun_d
    cmd = '%s -s minimize.tpr -x minimize.xtc -c %s -e minimize.edr -g minimize.log'%(mdrun, gro_file)
    runGmx(cmd)

    # Executes g_energy_d
    cmd = '/bin/echo %s 0 | %s -f minimize.edr -o energy.xvg'%(xvg_values, g_energy)
    runGmx(cmd)

    lines = readFile('energy.xvg')
    line = lines[-1]

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

def createLegend(xvg_values):
    """Creates a legend based on the chosen xvg values

    ARGUMENTS
        xvg_values (array) - the array containing the desired xvg_values

    RETURNS
        legend (array) - the array containing the desired legend
    """

    g_energy_legend = ['Time', 'Bond', 'Angle', 'Proper-Dih.', 'Ryckaert-Bell.', 'LG-14', 'Coulomb-14', 'LJ-(SR)', 'Coulomb-(SR)', 'Potential', 'Kinetic-En.', 'Total-Energy', 'Temperature', 'Pressure-(bar)', 'Vir-XX', 'Vir-XY', 'Vir-XZ', 'Vir-YX', 'Vir-YY', 'Vir-YZ', 'Vir-ZX', 'Vir-ZY', 'Vir-ZZ', 'Pres-XX-(bar)', 'Pres-XY-(bar)', 'Pres-XZ-(bar)', 'Pres-YX-(bar)', 'Pres-YY-(bar)', 'Pres-YZ-(bar)', 'Pres-ZX-(bar)', 'Pres-ZY-(bar)', 'Pres-ZZ-(bar)', '#Surf*SurfTen', 'Mu-X', 'Mu-Y', 'Mu-Z', 'T-rest']

    legend = ['Time']

    for i in xvg_values.split():
        legend.append(g_energy_legend[int(i)])

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

class checkrelativestates:

    def __init__(self, grompp, mdrun, g_energy, ligand_name_list, mdp, gros, originaltops, relativetops, originalitps = None, relativeitps = None, matches = None, ligand_basepath = None, debug = False):
        """Base class for checking relative states

        ARGUMENTS
            grompp (string) - full binary for grompp
            mdrum (string) - full binary for mdrun
            g_energy (string) - full binary for g_energy
            ligand_name_list (array) - filenames of the ligands
            mdp (string) - path to desired mdp
            gros (array) - paths specifying gro files
            originaltops (array) - paths specifying non-mutated top files
            relaivetops (array) - paths specifying relative  top files
            originalitps (array) - paths specifying non-mutated itp (optional)
            relativeitps (array) - paths specifying relative itp files (optional)
            matches (array) - paths specifying .mci files used for matching mutating atoms
            ligand_basepath (string) - path specifying location of ligands
            debug (boolean) - setting True keeps intermediate gromacs files
        """

        self.grompp = grompp
        self.mdrun = mdrun
        self.g_energy = g_energy
        self.ligand_name_list = ligand_name_list
        self.mdp = mdp
        self.gros = gros
        self.originaltops = originaltops
        self.relativetops = relativetops
        self.originalitps = originalitps
        self.relativeitps = relativeitps
        self.matches = matches
        self.ligand_basepath = ligand_basepath
        self.xvg_values = xvg_values
        self.debug = debug
        self.calculateEnergies()
        return

    def calculateEnergies(self):
        """Creates the directories for the ligand checking and calculates the energies of each ligand The comparisons are between the original state to relative state with no lambda and between two mutated states to a common substructure
        """

        base_directory = os.getcwd()
        commands.getoutput('mkdir intermediatefiles')
        os.chdir('intermediatefiles')
        working_directory = os.getcwd()

        # Finds array containing ligand values
        standard = self.findLigandEnergy(working_directory, 'initial')
        relativelambda0 = self.findLigandEnergy(working_directory, 'mutatedInitial')
        relativelambda1 = self.findLigandEnergy(working_directory, 'mutatedFinal')

        # Checks to see which ligands do not give the same energy and outputs the differences
        legend = createLegend(self.xvg_values)
        compareInitialLigands(standard, relativelambda0, legend, self.ligand_name_list)
        compareMutatedLigands(relativelambda1, legend, self.ligand_name_list)

        if self.debug == False:
            os.chdir(base_directory)
            commands.getoutput('rm -r intermediatefiles')

        return

    def findLigandEnergy(self, working_directory, state):
        """Runs findEnergy on all of the ligands

        ARGUMENTS
            working_directory (string) - the current working directory
            gromacs_directory (string) - the directory containing gromacs
            ligand_name_list (array of ) - the list of ligands to be analyzed
            xvg_values (array of integers) - the array containing the desired xvg_values
            state (string) - initial/mutatedInitial/mutatedFinal

        RETURNS
            ligand (array) - returns the array containing all of the energy terms for all of the comparisons for a single state
        """

        os.chdir(working_directory)
        n = len(self.xvg_values.split()) + 1
        ligand = array([])

        # Initial state calculation
        if state == 'initial':

            print 'Minimizing Initial'
            working_directory = os.path.join(working_directory, 'initial')
            if not os.path.exists(working_directory):
                commands.getoutput('mkdir %s'%(working_directory))

            # Runs findEnergy and appends values to array
            for i in range(len(self.ligand_name_list)):

                print 'Minimizing %s'%(self.ligand_name_list[i])

                new_working_directory = os.path.join(working_directory, self.ligand_name_list[i])
                if not os.path.exists(new_working_directory):
                    commands.getoutput('mkdir %s'%(new_working_directory))
                commands.getoutput('cp %s %s/minimize.mdp'%(self.mdp, new_working_directory))
                commands.getoutput('cp %s %s/system.gro'%(self.gros[i], new_working_directory))
                commands.getoutput('cp %s %s/system.top'%(self.originaltops[i], new_working_directory))
                if self.originalitps != None:
                    commands.getoutput('cp %s %s'%(self.originalitps[i], new_working_directory))
                os.chdir(new_working_directory)

                # Sets initial state parameters
                mdp = MdpFile('minimize.mdp')
                mdp.setParameter('free-energy', 'no')
                mdp.setParameter('nsteps', '0')
                mdp.write('minimize.mdp')

                line = findEnergy(self.grompp, self.mdrun, self.g_energy, 'minimize.mdp', 'system.gro', 'system.top', self.xvg_values)
                ligand = append(ligand.reshape(len(ligand), n), line.reshape(1, n), axis=0)

        # Mutated Initial calculation
        elif state == 'mutatedInitial':

            print 'Minimizing Initial Mutation'
            working_directory = os.path.join(working_directory,'initialmutation')
            if not os.path.exists(working_directory):
                commands.getoutput('mkdir %s'%(working_directory))

            expanded_ligand_list = list()
            expanded_gros_list = list()

            for i in range(len(self.ligand_name_list)):
                for j in range(i+1, len(self.ligand_name_list)):
                    expanded_ligand_list.append(self.ligand_name_list[i])
                    expanded_ligand_list.append(self.ligand_name_list[j])
                    expanded_gros_list.append(self.gros[i])
                    expanded_gros_list.append(self.gros[j])

            # Runs findEnergy and appends values to array
            for i in range(len(expanded_ligand_list)):

                if i%2 == 0:
                    comparison = '%s-%s'%(expanded_ligand_list[i], expanded_ligand_list[i+1])
                    working_directory = os.path.join(working_directory, comparison)
                    commands.getoutput('mkdir %s'%(working_directory))

                print 'Minimizing %s in pair %s'%(expanded_ligand_list[i], comparison)

                new_working_directory = os.path.join(working_directory, expanded_ligand_list[i])
                if not os.path.exists(new_working_directory):
                    commands.getoutput('mkdir %s'%(new_working_directory))
                commands.getoutput('cp %s %s/minimize.mdp'%(self.mdp, new_working_directory))
                commands.getoutput('cp %s %s/system.gro'%(expanded_gros_list[i], new_working_directory))
                commands.getoutput('cp %s %s/system.top'%(self.relativetops[i], new_working_directory))
                if self.relativeitps != None:
                    commands.getoutput('cp %s %s'%(self.relativeitps[i], new_working_directory))
                os.chdir(new_working_directory)

                # Sets mutated initial state parameters
                mdp = MdpFile('minimize.mdp')
                mdp.setParameter('free-energy', 'yes')
                mdp.setParameter('init-lambda', '0')
                mdp.setParameter('nsteps','0')
                mdp.write('minimize.mdp')

                line = findEnergy(self.grompp, self.mdrun, self.g_energy, 'minimize.mdp', 'system.gro', 'system.top', self.xvg_values)
                ligand = append(ligand.reshape(len(ligand), n), line.reshape(1, n), axis=0)

        # Mutated Final calculation
        elif state == 'mutatedFinal':

            print 'Minimizing Final Mutation'
            working_directory = os.path.join(working_directory,'finalmutation')
            if not os.path.exists(working_directory):
                commands.getoutput('mkdir %s'%(working_directory))

            expanded_ligand_list = list()
            expanded_gros_list = list()

            for i in range(len(self.ligand_name_list)):
                for j in range(i+1, len(self.ligand_name_list)):
                    expanded_ligand_list.append(self.ligand_name_list[i])
                    expanded_ligand_list.append(self.ligand_name_list[j])
                    expanded_gros_list.append(self.gros[i])
                    expanded_gros_list.append(self.gros[j])

            for i in range(len(expanded_ligand_list)):

                if i%2 == 0:
                    comparison = '%s-%s'%(expanded_ligand_list[i], expanded_ligand_list[i+1])
                    commands.getoutput('mkdir %s'%(os.path.join(working_directory, comparison)))

                print 'Minimizing %s in pair %s'%(expanded_ligand_list[i], comparison)

                new_working_directory = os.path.join(working_directory,comparison, expanded_ligand_list[i])
                if not os.path.exists(new_working_directory):
                    commands.getoutput('mkdir %s'%(new_working_directory))
                commands.getoutput('cp %s %s/minimize.mdp'%(self.mdp, new_working_directory))
                commands.getoutput('cp %s %s/system.gro'%(expanded_gros_list[i], new_working_directory))
                commands.getoutput('cp %s %s/system.top'%(self.relativetops[i], new_working_directory))
                if self.relativeitps != None:
                    commands.getoutput('cp %s %s'%(self.relativeitps[i], new_working_directory))
                os.chdir(new_working_directory)

                # Sets mutated final state parameters
                mdp = MdpFile('minimize.mdp')
                mdp.setParameter('free-energy', 'yes')
                mdp.setParameter('init-lambda', '1')
                mdp.setParameter('nsteps','0')
                mdp.write('minimize.mdp')

                if i%2 == 1:
                    if self.matches != None:
                        matchArray = readMci(self.matches[(i-1)/2])
                    else:
                        matchArray = determineMatchArray(self.ligand_basepath, expanded_ligand_list[i-1], expanded_ligand_list[i])
                    mutatedGroStructure = GromacsStructure()
                    mutatedGroStructure = mutateGroFile(expanded_gros_list[i-1], expanded_gros_list[i], matchArray)
                    mutatedGroStructure.write('system.gro')
                    modifySolventInTopFile(self.relativetops[i-1], self.relativetops[i], 'system.top')

                # Runs findEnergy and appends values to array
                line = findEnergy(self.grompp, self.mdrun, self.g_energy, 'minimize.mdp', 'system.gro', 'system.top', self.xvg_values)
                ligand = append(ligand.reshape(len(ligand), n), line.reshape(1, n ), axis=0)

        return ligand


