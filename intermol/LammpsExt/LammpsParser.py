import os
import warnings
import re
import pdb

import intermol.unit as units
from intermol.Atom import Atom
from intermol.Molecule import Molecule
from intermol.System import System
from intermol.Types import *
from intermol.Force import *
from intermol.HashMap import *

def writeData(filename, unit_set='real'):
    """Write system to LAMMPS data file

    Args:
        filename (str): name of output file
        unit_set (str): LAMMPS unit set for output file
    """
    # TODO: allow user to specify different unit sets
    if unit_set == 'real':
        DIST = units.angstroms
        VEL = units.angstroms / units.femtosecond
        ENERGY = units.kilocalorie
        MASS = units.grams / units.mole
        CHARGE = units.elementary_charge
        MOLE = units.mole

    # prepare a bunch of lists that we want to append to
    # these loosely correspond to GROMACS [ directives ]
    mass_list = list()
    mass_list.append('\n')
    mass_list.append('Masses\n')
    mass_list.append('\n')

    pair_coeff_list = list()
    pair_coeff_list.append('\n')
    pair_coeff_list.append('Pair Coeffs\n')
    pair_coeff_list.append('\n')

    bond_coeff_list = list()
    bond_coeff_list.append('\n')
    bond_coeff_list.append('Bond Coeffs\n')
    bond_coeff_list.append('\n')

    angle_coeff_list = list()
    angle_coeff_list.append('\n')
    angle_coeff_list.append('Angle Coeffs\n')
    angle_coeff_list.append('\n')

    dihedral_coeff_list = list()
    dihedral_coeff_list.append('\n')
    dihedral_coeff_list.append('Dihedral Coeffs\n')
    dihedral_coeff_list.append('\n')

    improper_coeff_list = list()
    improper_coeff_list.append('\n')
    improper_coeff_list.append('Improper Coeffs\n')
    improper_coeff_list.append('\n')

    atom_list = list()
    atom_list.append('\n')
    atom_list.append('Atoms\n')
    atom_list.append('\n')

    vel_list = list()
    vel_list.append('\n')
    vel_list.append('Velocities\n')
    vel_list.append('\n')

    bond_list = list()
    bond_list.append('\n')
    bond_list.append('Bonds\n')
    bond_list.append('\n')

    angle_list = list()
    angle_list.append('\n')
    angle_list.append('Angles\n')
    angle_list.append('\n')

    dihedral_list = list()
    dihedral_list.append('\n')
    dihedral_list.append('Dihedrals\n')
    dihedral_list.append('\n')

    improper_list = list()
    improper_list.append('\n')
    improper_list.append('Impropers\n')
    improper_list.append('\n')

    # read all atom specific information
    atom_type_dict = dict()  # str_type:int_type
    a_i = 1  # counter for atom types
    for moleculeType in System._sys._molecules.itervalues():
        for molecule in moleculeType.moleculeSet:
            for atom in molecule._atoms:
                # type, mass and pair coeffs
                if atom._atomtype[0] not in atom_type_dict:
                    atom_type_dict[atom._atomtype[0]] = a_i
                    mass_list.append('%d %8.4f\n'
                                % (a_i,
                                   atom._mass[0].in_units_of(MASS)._value))
                    pair_coeff_list.append('%d %8.4f %8.4f\n'
                                % (a_i,
                                   atom._sigma[0].in_units_of(DIST)._value,
                                   atom._epsilon[0].in_units_of(ENERGY/MOLE)._value))
                    a_i += 1
                # atom
                atom_list.append('%-6d %-6d %-6d %5.8f %8.3f %8.3f %8.3f\n'
                    % (atom.atomIndex,
                       atom.residueIndex,
                       atom_type_dict[atom._atomtype[0]],
                       atom._charge[0].in_units_of(CHARGE)._value,
                       atom._position[0].in_units_of(DIST)._value,
                       atom._position[1].in_units_of(DIST)._value,
                       atom._position[2].in_units_of(DIST)._value))
                # velocity
                vel_list.append('%-6d %8.4f %8.4f %8.4f\n'
                    % (atom.atomIndex,
                       atom._velocity[0].in_units_of(VEL)._value,
                       atom._velocity[1].in_units_of(VEL)._value,
                       atom._velocity[2].in_units_of(VEL)._value))

    # read all connectivity information
    bond_style = 0
    bond_type_dict = dict()
    bond_numbering = dict()
    b_type_i = 1
    b_mol_i = 1
    b_i = 1

    for moleculeType in System._sys._molecules.itervalues():
        # atom index offsets from 1 for each molecule
        offsets = list()
        for molecule in moleculeType.moleculeSet:
            offsets.append(molecule._atoms[0].atomIndex - 1)

        # bond types
        if moleculeType.bondForceSet:
            for bond in moleculeType.bondForceSet.itervalues():
                if isinstance(bond, Bond):
                    if bond_style == 0:
                        bond_style = 'harmonic'
                    elif bond_style != 'harmonic':
                        warnings.warn("More than one bond style found!")
                        # this may need to be an error

                    atomtype1 = moleculeType.moleculeSet[0]._atoms[bond.atom1 - 1].getAtomType()[0]
                    atomtype2 = moleculeType.moleculeSet[0]._atoms[bond.atom2 - 1].getAtomType()[0]

                    temp = BondType(atomtype1,
                            atomtype2,
                            1,
                            bond.length,
                            bond.k)
                    if temp not in bond_type_dict.itervalues():
                        bond_coeff_list.append('%d %18.8e %18.8e\n'
                                    % (b_type_i,
                                       bond.length.in_units_of(DIST)._value,
                                       bond.k.in_units_of(ENERGY/(DIST*DIST * MOLE))._value))

                else:
                    warnings.warn("Found unsupported bond type for LAMMPS!")

                if temp not in bond_type_dict.itervalues():
                    bond_type_dict[b_type_i] = temp
                    b_type_i += 1
                    bond_numbering[b_mol_i] = b_type_i
                    b_mol_i += 1
                else:
                    bond_numbering[]

            pdb.set_trace()
            for offset in offsets:
                for per_mol, bond in enumerate(moleculeType.bondForceSet.itervalues()):
                    bond_list.append('%6d %6d %6d %6d\n'
                                % (b_i,
                                   bond_numbering[per_mol + 1],
                                   bond.atom1 + offset,
                                   bond.atom2 + offset))
                    b_i += 1


    pdb.set_trace()

    # actual data file writing
    with open(filename, 'w') as f:
        f.write(System._sys._name + '\n')
        f.write('\n')

        n_atoms = len(atom_list) - 3
        n_bonds = len(bond_list) - 3
        n_angles = len(angle_list) - 3
        n_dihedrals = len(dihedral_list) - 3
        n_impropers = len(improper_list) - 3

        n_atom_types = len(pair_coeff_list) - 3
        n_bond_types = len(bond_coeff_list) - 3
        n_angle_types = len(angle_coeff_list) - 3
        n_dihedral_types = len(dihedral_coeff_list) - 3
        n_improper_types = len(improper_coeff_list) - 3

        f.write('{0} atoms\n'.format(n_atoms))
        f.write('{0} bonds\n'.format(n_bonds))
        f.write('{0} angles\n'.format(n_angles))
        f.write('{0} dihedrals\n'.format(n_dihedrals))
        f.write('{0} impropers\n'.format(n_impropers))
        f.write('\n')

        f.write('{0} atom types\n'.format(n_atom_types))
        if n_bonds > 0:
            f.write('{0} bond types\n'.format(n_bond_types))
        if n_angles > 0:
            f.write('{0} angle types\n'.format(n_angle_types))
        if n_dihedrals > 0:
            f.write('{0} dihedral types\n'.format(n_dihedral_types))
        if n_impropers > 0:
            f.write('{0} improper types\n'.format(n_improper_types))
        f.write('\n')

        f.write('0.0 %10.6f xlo xhi\n' % (System._sys._v1x.in_units_of(DIST)._value))
        f.write('0.0 %10.6f ylo yhi\n' % (System._sys._v2y.in_units_of(DIST)._value))
        f.write('0.0 %10.6f zlo zhi\n' % (System._sys._v3z.in_units_of(DIST)._value))

        # write masses
        for mass in mass_list:
            f.write(mass)

        # write coefficients
        for pair in pair_coeff_list:
            f.write(pair)
        if len(bond_coeff_list) > 3:
            for bond in bond_coeff_list:
                f.write(bond)
        if len(angle_coeff_list) > 3:
            for angle in angle_coeff_list:
                f.write(pair)
        if len(dihedral_coeff_list) > 3:
            for dihedral in dihedral_coeff_list:
                f.write(dihedral)
        if len(improper_coeff_list) > 3:
            for improper in improper_coeff_list:
                f.write(improper)

        # write atoms and velocities
        for atom in atom_list:
            f.write(atom)
        for vel in vel_list:
            f.write(vel)

        # write connectivity
        if len(bond_list) > 3:
            for bond in bond_list:
                f.write(bond)
        if len(angle_list) > 3:
            for angle in angle_list:
                f.write(angle)
        if len(dihedral_list) > 3:
            for dihedral in dihedral_list:
                f.write(dihedral)
        if len(improper_list) > 3:
            for improper in improper_list:
                f.write(improper)

    # write the corresponding input file
    basename = os.path.splitext(filename)[0]
    input_filename = basename + '.input'
    with open(input_filename, 'w') as f:
        f.write('units {0}\n'.format(unit_set))
        f.write('atom_style full\n')  # TODO
        f.write('\n')

        f.write('dimension 3\n')  # TODO
        f.write('boundary p p p\n')  # TODO
        f.write('\n')

        # non-bonded
        f.write('pair_style lj/cut/coul/long 10.0\n')  # TODO: match mdp
        f.write('pair_modify mix geometric\n')  # TODO: match defaults
        f.write('kspace_style pppm 1.0e-4\n')  # TODO: match mdp
        f.write('\n')

        # bonded
        if len(bond_coeff_list) > 3:
            f.write('bond_style {0}\n'.format(bond_style))
        if len(angle_coeff_list) > 3:
            f.write('angle_style {0}\n'.format(angle_style))
        if len(dihedral_coeff_list) > 3:
            f.write('dihedral_style {0}\n'.format(dihedral_style))
        if len(improper_coeff_list) > 3:
            f.write('improper_style {0}\n'.format(improper_style))
        f.write('\n')

        # read data
        f.write('read_data {0}\n'.format(os.path.basename(filename)))
        f.write('\n')

        # output energies
        energy_terms = " ".join(['ebond',
                                 'eangle',
                                 'edihed',
                                 'epair',
                                 'evdwl',
                                 'ecoul',
                                 #'epoten',
                                 #'ekinet',
                                 'etotal',
                                 'temp'])

        f.write('thermo_style custom {0}\n'.format(energy_terms))
        f.write('\n')

        f.write('run 0\n')

        # special bonds? CHARMM?
        # multiple types in one category?
