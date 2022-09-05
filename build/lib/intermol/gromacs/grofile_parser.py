import logging

from parmed.unit import nanometers, picoseconds
import numpy as np

logger = logging.getLogger('InterMolLog')


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
        atomname = list()
        resid = list()
        resname = list()
        boxes = list()
        xyzs = list()
        vels = list()
        with open(self.gro_file) as gro:
            next(gro)
            n_atoms = int(next(gro).strip())

            for _ in range(n_atoms):
                line = next(gro)
                (thisresnum, thisresname, thisatomname) = [line[i*5:i*5+5].strip() for i in range(3)]
                resname.append(thisresname)
                resid.append(int(thisresnum))
                atomname.append(thisatomname)

                entries = line[20:].split()
                # If there aren't 6, then fixed column, presumably 8 digit
                if len(entries) not in [3, 6]:
                    data = line[20:]
                    entries = []
                    spacing = 8
                    for j in range(0, len(data), spacing):
                        entry = data[j:j+spacing].strip()
                        if len(entry) > 0:
                            entries.append(entry)
                entries = [float(x) for x in entries]
                xyz = [x * nanometers for x in entries[:3]]
                xyzs.append(xyz)
                if len(entries) == 6:
                    vel = [v * nanometers / picoseconds for v in entries[3:6]]
                else:
                    vel = [v * nanometers / picoseconds for v in [0., 0., 0.]]
                vels.append(vel)

            line = next(gro)
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

        self.positions = np.array(xyzs)
        self.velocities = np.array(vels)
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
                # .gro wraps at 100,0000, which is why the field is 5 width.    
                gro.write('{0:5d}{1:<5s}{2:5s}{3:5d}'.format(
                        atom.residue_index%100000, atom.residue_name, atom.name, (n + 1)%100000))
                for pos in atom.position:
                    gro.write('{0:17.12f}'.format(pos.value_in_unit(nanometers)))
                if np.any(atom.velocity):
                    for vel in atom.velocity:
                        gro.write('{0:17.12f}'.format(vel.value_in_unit(nanometers / picoseconds)))
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

