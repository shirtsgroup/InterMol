import logging
from re import match

from simtk.unit import nanometers
import numpy as np

logger = logging.getLogger('InterMolLog')


def _isint(word):
    """ONLY matches integers! If you have a decimal point? None shall pass!

    @param[in] word String (for instance, '123', '153.0', '2.', '-354')
    @return answer Boolean which specifies whether the string is an integer (only +/- sign followed by digits)

    """
    return match('^[-+]?[0-9]+$', word)


def _isfloat(word):
    """Matches ANY number; it can be a decimal, scientific notation, what have you
    CAUTION - this will also match an integer.

    @param[in] word String (for instance, '123', '153.0', '2.', '-354')
    @return answer Boolean which specifies whether the string is any number

    """
    return match('^[-+]?[0-9]*\.?[0-9]*([eEdD][-+]?[0-9]+)?$', word)


def _is_gro_coord(line):
    """ Determines whether a line contains GROMACS data or not

    @param[in] line The line to be tested

    """
    sline = line.split()
    if len(sline) == 6 or len(sline) == 9:
        return all([_isint(sline[2]), _isfloat(sline[3]), _isfloat(sline[4]), _isfloat(sline[5])])
    elif len(sline) == 5 or len(sline) == 8:
        return all([_isint(line[15:20]), _isfloat(sline[2]), _isfloat(sline[3]), _isfloat(sline[4])])
    else:
        return 0


def _is_gro_box(line):
    """ Determines whether a line contains a GROMACS box vector or not

    @param[in] line The line to be tested

    """
    sline = line.split()
    if len(sline) == 9 and all([_isfloat(i) for i in sline]):
        return 1
    elif len(sline) == 3 and all([_isfloat(i) for i in sline]):
        return 1
    else:
        return 0


class GromacsGroParser(object):
    """GromacsGroParser reads and writes Gromacs .gro files

    A .gro file also contains some topological information, such as elements and
    residue names, but not enough to construct a full Topology object.  This
    information is recorded and stored in the object's public fields.
    """

    def __init__(self, gro_file):
        """Load a .gro gro_file.

        The atom positions can be retrieved by calling getPositions().

        Parameters:
         - gro_file (string) the name of the gro_file to read or write
        """
        self.gro_file = gro_file

    def read(self):
        xyzs = list()
        atomname = list()
        comms = list()
        resid = list()
        resname = list()
        boxes = list()
        xyz = list()
        line_number = 0
        frame = 0
        for line in open(self.gro_file):
            if line_number == 0:
                comms.append(line.strip())
            elif line_number == 1:
                n_atoms = int(line.strip())
            elif _is_gro_coord(line):
                if frame == 0: # Create the list of residues, atom names etc. only if it's the first frame.
                    (thisresnum, thisresname, thisatomname) = [line[i*5:i*5+5].strip() for i in range(3)]
                    resname.append(thisresname)
                    resid.append(int(thisresnum))
                    atomname.append(thisatomname)
                first_decimal_pos = line.index('.', 20)
                second_decimal_pos = line.index('.', first_decimal_pos+1)
                digits = second_decimal_pos-first_decimal_pos
                pos = [float(line[20+i*digits:20+(i+1)*digits]) for i in range(3)]
                xyz.append([pos[0], pos[1], pos[2]])
            elif _is_gro_box(line) and line_number == n_atoms + 2:
                raw_box_vector = line.split()
                v = np.zeros([3, 3], float) * nanometers
                # Diagonals
                for i in range(3):
                    v[i, i] = float(raw_box_vector[i]) * nanometers
                if len(raw_box_vector) == 9:
                    k = 3
                    # Then the off-diagonals
                    for i in range(3):
                        for j in range(3):
                            if i != j:
                                v[i, j] = float(raw_box_vector[k]) * nanometers
                                k += 1
                boxes.append(v)
                xyzs.append(xyz*nanometers)
                xyz = list()
                line_number = -1
                frame += 1
            else:
                e = Exception("Unexpected line in .gro gro_file: {0}".format(line))
                logger.exception(e)
            line_number += 1

        self.positions = np.array(xyzs[0])
        self.atom_names = atomname
        self.residue_ids = resid
        self.residue_names = resname
        self.box_vector = boxes[0]

    def write(self, system):
        """Write the system out  in a Gromacs 4.6 format

        Args:
            filename (str): the file to write out to
        """
        with open(self.gro_file, 'w') as gro:
            gro.write("{0}\n".format(system.name))
            gro.write("{0}\n".format(system.n_atoms))
            for n, atom in enumerate(system.atoms):
                if atom.name.isdigit():
                    # Kluge for atoms read in from a LAMMPS data file.
                    atom.name = "LMP_{0}".format(atom.name)
                gro.write('{0:5d}{1:<4s}{2:6s}{3:5d}'.format(
                        atom.residue_index, atom.residue_name, atom.name, n))
                for pos in atom.position:
                    gro.write('{0:17.12f}'.format(pos.value_in_unit(nanometers)))
                gro.write('\n')

            # Check for rectangular; should be symmetric, so we don't have to
            # check 6 values
            if (system.box_vector[1, 0]._value == 0 and
                system.box_vector[2, 0]._value == 0 and
                system.box_vector[2, 1]._value == 0):
                    for i in range(3):
                        gro.write('{0:11.7f}'.format(system.box_vector[i, i].value_in_unit(nanometers)))
            else:
                for i in range(3):
                    gro.write('{0:11.7f}'.format(system.box_vector[i, i].value_in_unit(nanometers)))
                for i in range(3):
                    for j in range(3):
                        if i != j:
                            gro.write('{0:11.7f}'.format(system.box_vector[i, j].value_in_unit(nanometers)))

            gro.write('\n')


if __name__ == "__main__":
    import pdb
    gro = GromacsGroParser('../../tests/gromacs/unit_tests/bond1/bond1.gro')

    pdb.set_trace()
