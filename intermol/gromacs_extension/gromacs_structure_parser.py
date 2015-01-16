import numpy as np
from intermol.atom import *
from intermol.molecule import *
from intermol.system import System
import logging

logger = logging.getLogger('InterMolLog')

def readStructure(filename):
    """Read in a Gromacs structure file

    Args:
        filename (str): the file to be read in
    """

    fd = open(filename, 'r')
    lines = list(fd)
    fd.close()

    i = 2
    uniquec = set()
    # get a list of the unique components, so that we can tell how many of each of them we printed out if
    # they are split into multiple groups in the .top.

    for component in System._sys._components:
        uniquec.add(component[0])
    componentcount = dict(zip(uniquec, np.zeros(len(uniquec), int)))

    for component in System._sys._components:
        moleculetype = System._sys._molecules[component[0]]
        molecules = moleculetype.moleculeSet.list
        ncomponent = componentcount[component[0]]
        for n in range(ncomponent, ncomponent+component[1]):
            molecule = molecules[n]
            for atom in molecule._atoms:
                if lines[i]:
# NOTE: These should already be defined (see gromacs_topology_parser)
#                    atom.residue_index = int(lines[i][0:5])
#                    atom.residue_name = lines[i][5:10].strip()
#                    atom.name = lines[i][10:15].strip()
                    variables = (lines[i][20:]).split()
                    position = np.zeros([3], float) * units.nanometers
                    velocity = np.zeros([3], float) * units.nanometers / units.picoseconds
                    if len(variables) >= 3:
                        positionElements = variables[0:3]
                        for k in range(3):
                            position[k] = float(positionElements[k]) * units.nanometers
                    atom.setPosition(position[0], position[1], position[2])
                    if len(variables) >= 6:
                        velocity_elements = variables[3:6]
                        for k in range(3):
                            velocity[k] = float(velocity_elements[k]) * units.nanometers / units.picoseconds
                    atom.setVelocity(velocity[0], velocity[1], velocity[2])
                    i += 1
                else:
                    sys.exit()
        componentcount[component[0]] += component[1]

    raw_box_vector = lines[i].split()
    v = np.zeros([3, 3], float) * units.nanometers
    # diagonals
    for i in range(3):
        v[i, i] = float(raw_box_vector[i]) * units.nanometers

    if len(raw_box_vector) == 9:
        k = 3
        # Then the off-diagonals
        for i in range(3):
            for j in range(3):
                if  i != j:
                    v[i, j] = float(raw_box_vector[k]) * units.nanometers
                    k += 1

    # need to make this numpy sensitive so we can just pass the vector
    System._sys.box_vector = v


def writeStructure(filename):
    """Write the system out  in a Gromacs 4.6 format

    Args:
        filename (str): the file to write out to
    """
    lines = list()
    n = 0
    res_offset = 0
    res_offset_next = 0
    for moleculetype in System._sys._molecules.values():
        for molecule in moleculetype.moleculeSet:
            for atom in molecule._atoms:
#                if atom.name.isdigit():
#                    atom.name = "LMP_" + atom.name
                n += 1
                lines.append('%5d%-5s%5s%5d%13.8f%13.8f%13.8f%13.8f%13.8f%13.8f\n'
                             % (atom.residue_index + res_offset, atom.residue_name, atom.name, n,
                                atom._position[0].in_units_of(units.nanometers)._value,
                                atom._position[1].in_units_of(units.nanometers)._value,
                                atom._position[2].in_units_of(units.nanometers)._value,
                                atom._velocity[0].in_units_of(units.nanometers / units.picoseconds)._value,
                                atom._velocity[1].in_units_of(units.nanometers / units.picoseconds)._value,
                                atom._velocity[2].in_units_of(units.nanometers / units.picoseconds)._value))
                res_offset_next = res_offset + atom.residue_index
            res_offset = res_offset_next
    # print the box vector
    # check for rectangular; should be symmetric, so we don't have to check 6 values
    if (System._sys.box_vector[1, 0]._value == 0 and
        System._sys.box_vector[2, 0]._value == 0 and
        System._sys.box_vector[2, 1]._value == 0):
            for i in range(3):
                lines.append('%11.7f ' % System._sys.box_vector[i, i].in_units_of(units.nanometers)._value)
    else:
        for i in range(3):
            lines.append('%11.7f ' % System._sys.box_vector[i, i].in_units_of(units.nanometers)._value)
        for i in range(3):
            for j in range(3):
                if i != j:
                    lines.append('%11.7f ' % System._sys.box_vector[i, j].in_units_of(units.nanometers)._value)

    lines.append('\n')

    lines.insert(0, (System._sys._name + '\n'))
    lines.insert(1, str(n) + '\n')

    fout = open(filename, 'w')
    for line in lines:
        fout.write(line)
    fout.close()
