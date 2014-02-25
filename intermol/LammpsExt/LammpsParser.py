import os
import warnings
import re
import numpy as np
import sys
import pdb

import intermol.unit as units
from intermol.Atom import Atom
from intermol.Molecule import Molecule
from intermol.System import System
from intermol.Types import *
from intermol.Force import *
from intermol.HashMap import *

def readData(data_file, input_file, verbose=False):
    """Reads a LAMMPS data file

    *** Only works for directives delimited by blank lines ***

    Currently supports the following directives:
        Masses
        Pair Coeffs (must be mix geometry)
        Bond Coeffs (must be harmonic)
        Angle Coeffs (must be harmonic)
        Dihedral Coeffs (must be OPLS)
        Atoms
        Bonds
        Angles
        Dihedrals

    TODO:
        -handling for comments
        -handling for directives not delimited by blank lines
        -allow specification of forcefield styles

    Args:
        data_file (str): name of LAMMPS data file to read in
        input_file (str): name of LAMMPS input file to read in
    """
    # read input file...
    with open(input_file, 'r') as f:
        input_lines = f.readlines()

    # NOTE: there are likely several ways to make this regex more robust and
    #       actually match whatever parsing the actual LAMMPS engine does
    directives = re.compile(r"""
        ((?P<units>^\s*units\s+\w+)
        |
        (?P<atom_style>^\s*atom_style\s+\w+)
        |
        (?P<dimension>^\s*dimension\s+\d)
        |
        (?P<boundary>^\s*boundary\s+\w{1,2}\s+\w{1,2}\s+\w{1,2})
        |
        (?P<pair_style>^\s*pair_style\s+\w+)
        |
        (?P<kspace_style>^\s*kspace_style\s+\w+)
        |
        (?P<pair_modify>^\s*pair_modify\s+\w+)
        |
        (?P<bond_style>^\s*bond_style\s+\w+)
        |
        (?P<angle_style>^\s*angle_style\s+\w+)
        |
        (?P<dihedral_style>^\s*dihedral_style\s+\w+)
        |
        (?P<improper_style>^\s*improper_style\s+\w+\s+)
        |
        (?P<special_bonds>^\s*special_bonds\s+\w+)
        |
        (?P<Impropers>\s*Impropers))
        """, re.VERBOSE)

    i = 0
    while i < len(input_lines):
        match = directives.match(input_lines[i])
        if match:
            if match.group('units'):
                fields = input_lines.pop(i).split()
                unit_set = fields[1]

            elif match.group('atom_style'):
                # NOTE: assuming 'full' as default for everything else
                # TODO: support more options
                fields = input_lines.pop(i).split()
                if len(fields) == 2:
                    atom_style = fields[1]
                else:
                    warnings.warn("Found unsupported atom_style in input file!")

            elif match.group('dimension'):
                fields = input_lines.pop(i).split()
                dimension = int(fields[1])
                if dimension not in [2, 3]:
                    raise Exception("Invalid dimension specified in input file (must be 2 or 3).")

            elif match.group('boundary'):
                fields = input_lines.pop(i).split()
                boundaries = [fields[1], fields[2], fields[3]]

            elif match.group('pair_style'):
                # TODO: support more options
                fields = input_lines.pop(i).split()
                pair_style = []
                if fields[1] == 'hybrid':
                    pass
                    # TODO: parse multiple pair styles
                elif fields[1] == 'lj/cut/coul/long':
                    pair_style.append(fields[1])
                    System._sys._nbFunc = 1

            elif match.group('kspace_style'):
                # TODO: currently ignored
                fields = input_lines.pop(i).split()
                if fields[1] == 'pppm':
                    pass

            elif match.group('pair_modify'):
                # TODO: support more options
                fields = input_lines.pop(i).split()
                if fields[1] == 'mix':
                    if fields[2] == 'geometric':
                        System._sys._combinationRule = 3
                    elif fields[2] == 'arithmetic':
                        System._sys._combinationRule = 2
                    else:
                        warnings.warn("Found unsupported `pair_modify mix` argument in input file!")
                else:
                    warnings.warn("Found unsupported `pair_modify` style in input file!")

            elif match.group('bond_style'):
                fields = input_lines.pop(i).split()
                bond_style = []
                if len(fields) == 2:
                    bond_style.append(fields[1])
                elif fields[1] == 'hybrid':
                    bond_style = []
                    for style in fields[2:]:
                        bond_style.append(style)
                else:
                    raise Exception("Invalid bond_style in input file!")

            elif match.group('angle_style'):
                fields = input_lines.pop(i).split()
                angle_style = []
                if len(fields) == 2:
                    angle_style.append(fields[1])
                elif fields[1] == 'hybrid':
                    angle_style = []
                    for style in fields[2:]:
                        angle_style.append(style)
                else:
                    raise Exception("Invalid angle_style in input file!")

            elif match.group('dihedral_style'):
                fields = input_lines.pop(i).split()
                dihedral_style = []
                if len(fields) == 2:
                    dihedral_style.append(fields[1])
                    # TODO: correctly determine gen-pairs state
                    if dihedral_style == 'opls':
                        System._sys._genpairs = 'yes'
                elif fields[1] == 'hybrid':
                    dihedral_style = []
                    for style in fields[2:]:
                        dihedral_style.append(style)
                else:
                    raise Exception("Invalid dihedral_style in input file!")

            elif match.group('improper_style'):
                fields = input_lines.pop(i).split()
                improper_style = []
                if len(fields) == 2:
                    improper_style.append(fields[1])
                elif fields[1] == 'hybrid':
                    improper_style = []
                    for style in fields[2:]:
                        improper_style.append(style)
                else:
                    raise Exception("Invalid improper_style in input file!")

            elif match.group('special_bonds'):
                # TODO: support more options
                fields = input_lines.pop(i).split()
                if 'lj/coul' in fields:
                    System._sys._ljCorrection = float(fields[fields.index('lj/coul') + 3])
                    System._sys._coulombCorrection = float(fields[fields.index('lj/coul') + 3])
                elif 'lj' in fields:
                    System._sys._ljCorrection = float(fields[fields.index('lj') + 3])
                elif 'coul' in fields:
                    System._sys._coulombCorrection = float(fields[fields.index('coul') + 3])
            else:
                warnings.warn("Found unused regex match in input file!")
                i += 1
        else:
            i += 1

    # define units
    RAD = units.radians
    DEGREE = units.degrees
    if unit_set == 'real':
        DIST = units.angstroms
        VEL = units.angstroms / units.femtosecond
        ENERGY = units.kilocalorie / units.mole
        MASS = units.grams / units.mole
        CHARGE = units.elementary_charge
        MOLE = units.mole
    else:
        raise Exception("Unsupported unit set specified in input file: {0}".format(unit_set))

    # read type info from data file...
    with open(data_file, 'r') as f:
        data_lines = f.readlines()

    # TODO: improve robustness of box regex
    directives = re.compile(r"""
        ((?P<box>.+xlo)
        |
        (?P<Masses>\s*Masses)
        |
        (?P<PairCoeffs>\s*Pair\sCoeffs)
        |
        (?P<BondCoeffs>\s*Bond\sCoeffs)
        |
        (?P<AngleCoeffs>\s*Angle\sCoeffs)
        |
        (?P<DihedralCoeffs>\s*Dihedral\sCoeffs)
        |
        (?P<ImproperCoeffs>\s*Improper\sCoeffs))
        """, re.VERBOSE)

    i = 0
    while i < len(data_lines):
        match = directives.match(data_lines[i])
        if match:
            if match.group('box'):
                # TODO: support for non-orthogonal boxes
                v = np.zeros([3,3], float)
                for j in range(3):
                    fields = [float(field) for field in data_lines.pop(i).split()[:2]]
                    box_length = fields[1] - fields[0]
                    if box_length > 0:
                        v[j, j] = box_length
                    else:
                        raise Exception("Negative box length specified in data file.")
                System._sys.setBoxVector(v * DIST)

            elif match.group('Masses'):
                if verbose:
                    print 'Parsing Masses...'
                data_lines.pop(i)
                data_lines.pop(i)

                mass_dict = dict()  # type:mass
                #     not end of file         not blank line
                while i < len(data_lines) and data_lines[i].strip():
                    fields = data_lines.pop(i).split()
                    mass_dict[int(fields[0])] = float(fields[1]) * MASS

            elif match.group('PairCoeffs'):
                nb_types = dict()
                if verbose:
                    print 'Parsing Pair Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = data_lines.pop(i).split()
                    # TODO: support more options
                    if len(pair_style) == 1:
                        # TODO: may have to do lookup of pairstyle to determine format
                        if System._sys._nbFunc == 1:
                            nb_types[int(fields[0])] = [float(fields[1]) * ENERGY,
                                                        float(fields[2]) * DIST]
                        else:
                            warnings.warn("Found unsupported pair coeff formatting in data file!")
                    else:
                        warnings.warn("Found unsupported pair coeff formatting in data file!")

            elif match.group('BondCoeffs'):
                bond_types = dict()
                if verbose:
                    print 'Parsing Bond Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = [float(field) for field in data_lines.pop(i).split()]
                    # TODO: support more options
                    if len(bond_style) == 1:
                        if 'harmonic' in bond_style:
                            bond_types[int(fields[0])] = [2 * fields[1] * ENERGY / (DIST * DIST),
                                                          fields[2] * DIST]
                        else:
                            warnings.warn("Found unsupported bond coeff formatting in data file!")
                    else:
                        warnings.warn("Found unsupported bond coeff formatting in data file!")

            elif match.group('AngleCoeffs'):
                angle_types = dict()
                if verbose:
                    print 'Parsing Angle Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = [float(field) for field in data_lines.pop(i).split()]
                    # TODO: support more options
                    if len(angle_style) == 1:
                        if 'harmonic' in angle_style:
                            angle_types[int(fields[0])] = [2 * fields[1] * ENERGY / (RAD * RAD),
                                                           fields[2] * DEGREE]
                        else:
                            warnings.warn("Found unsupported angle coeff formatting in data file!")
                    else:
                        warnings.warn("Found unsupported angle coeff formatting in data file!")

            elif match.group('DihedralCoeffs'):
                dihedral_types = dict()
                if verbose:
                    print 'Parsing Dihedral Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = [float(field) for field in data_lines.pop(i).split()]
                    if len(dihedral_style) == 1:
                        if 'opls' in dihedral_style:
                            dihedral_types[int(fields[0])] = [fields[1] * ENERGY,
                                                              fields[2] * ENERGY,
                                                              fields[3] * ENERGY,
                                                              fields[4] * ENERGY]
                        else:
                            warnings.warn("Found unsupported dihedral coeff formatting in data file!")
                    else:
                        warnings.warn("Found unsupported dihedral coeff formatting in data file!")

            elif match.group('ImproperCoeffs'):
                improper_types = dict()
                if verbose:
                    print 'Parsing Improper Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = [float(field) for field in data_lines.pop(i).split()]
                    # TODO: not implemented...
                    warnings.warn("Found unsupported dihedral coeff formatting in data file!")
            else:
                warnings.warn("Found unused regex match in data file!")
                i += 1
        else:
            i += 1

    # NOTE: current implementation treats the entire LAMMPS system as one molecule
    moleculeName = data_lines[0].strip()
    currentMolecule = Molecule(moleculeName)
    System._sys.addMolecule(currentMolecule)
    currentMoleculeType = System._sys._molecules[moleculeName]
    currentMoleculeType.nrexcl = 3 # TODO: figure this number out based on lj/coul fudge

    directives = re.compile(r"""
        ((?P<n_atoms>\s*\d+\s+atoms)
        |
        (?P<n_bonds>\s*\d+\s+bonds)
        |
        (?P<n_angles>\s*\d+\s+angles)
        |
        (?P<n_dihedrals>\s*\d+\s+dihedrals)
        |
        (?P<Atoms>\s*Atoms)
        |
        (?P<Bonds>\s*Bonds)
        |
        (?P<Angles>\s*Angles)
        |
        (?P<Dihedrals>\s*Dihedrals)
        |
        (?P<Impropers>\s*Impropers))
        """, re.VERBOSE)

    i = 0
    while i < len(data_lines):
        match = directives.match(data_lines[i])
        if match:
            if match.group('n_atoms'):
                fields = data_lines.pop(i).split()
                n_atoms = int(fields[0])
                xyz = np.empty(shape=(n_atoms, 3))
                types = np.empty(shape=(n_atoms), dtype='int')
                masses = np.empty(shape=(n_atoms))
                charges = np.empty(shape=(n_atoms))

            elif match.group('n_bonds'):
                fields = data_lines.pop(i).split()
                bonds = np.empty(shape=(float(fields[0]), 3), dtype='int')

            elif match.group('n_angles'):
                fields = data_lines.pop(i).split()
                angles = np.empty(shape=(float(fields[0]), 4), dtype='int')

            elif match.group('n_dihedrals'):
                fields = data_lines.pop(i).split()
                dihedrals = np.empty(shape=(float(fields[0]), 5), dtype='int')

            elif match.group('Atoms'):
                if verbose:
                    print 'Parsing Atoms...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = data_lines.pop(i).split()
                    if len(fields) in [7, 10]:
                        if len(fields) == 10:
                            # TODO: store image flags?
                            pass

                        newAtomType = None
                        if System._sys._combinationRule == 1:
                            pass

                        elif System._sys._combinationRule in [2, 3]:
                            # NOTE: assuming 'A' as default ptype
                            newAtomType = AtomCR23Type(fields[2],  # atomtype or name
                                    fields[2],                     # bondtype
                                    -1,                            # Z
                                    mass_dict[int(fields[2])],     # mass
                                    float(fields[3]) * CHARGE,     # charge
                                    'A',                           # ptype
                                    nb_types[int(fields[2])][1],   # sigma
                                    nb_types[int(fields[2])][0])   # epsilon

                        System._sys._atomtypes.add(newAtomType)

                        atom = Atom(int(fields[0]),  # AtomNum
                                fields[2],           # atomName
                                int(fields[1]),      # resNum (molNum for LAMMPS)
                                fields[1])           # resName (molNum for LAMMPS)
                        atom.setAtomType(0, fields[2])  # atomNum for LAMMPS
                        atom.cgnr = 0  # TODO: look into alternatives because this could cause performance issues
                        atom.setCharge(0, float(fields[3]) * CHARGE)
                        atom.setMass(0, mass_dict[int(fields[2])])
                        atom.setPosition(float(fields[4]) * DIST,
                                float(fields[5]) * DIST,
                                float(fields[6]) * DIST)

                        for ab_state, atomType in enumerate(atom._atomtype):
                            # Searching for a matching atomType to pull values from
                            tempType = AbstractAtomType(atom._atomtype[ab_state])
                            atomType = System._sys._atomtypes.get(tempType)
                            if atomType:
                                atom.setSigma(ab_state, atomType.sigma)
                                atom.setEpsilon(ab_state, atomType.epsilon)
                                atom.bondtype =  atomType.bondtype
                            else:
                                warnings.warn("A corresponding AtomType was not found. Insert missing values yourself.\n")
                        currentMolecule.addAtom(atom)

            elif match.group('Bonds'):
                if verbose:
                    print 'Parsing Bonds...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = [int(field) for field in data_lines.pop(i).split()]
                    #bonds[fields[0] - 1] = fields[1:]

                    newBondForce = None
                    atomtype1 = currentMolecule._atoms[fields[2]-1].bondtype
                    atomtype2 = currentMolecule._atoms[fields[3]-1].bondtype

                    # TODO: implement bond type creation
                    """
                    if len(bond_style) == 1:
                        if bond_style[0] == 'harmonic':
                            tempBondStyle = 1
                    tempType = AbstractBondType(atomtype1, atomtype2, tempBondStyle)
                    bondType = bond_types.get(tempType)
                    if not bondType:
                        # we only have the reversed bond order stored, flip the atoms
                        tempType = AbstractBondType(atomtype2, atomtype1, tempBondStyle)
                        bondType = self.bondtypes.get(tempType)
                    if not bondType:
                        raise Exception("Bondtype lookup failed for '{0}'".format(" ".join(split)))

                    if isinstance(bondType, BondType):
                        split.append(bondType.length)
                        split.append(bondType.k)

                    else:
                        warnings.warn("Bondtype '{0}' is unsupported or something more complicated went wrong".format(bondType.type))
                    """
                    if len(bond_style) == 1:
                        if bond_style[0] == 'harmonic':
                            fields.append(bond_types[int(fields[1])][1])
                            fields.append(bond_types[int(fields[1])][0])
                            newBondForce = Bond(fields[2],
                                    fields[3],
                                    fields[4],
                                    fields[5])
                        currentMoleculeType.bondForceSet.add(newBondForce)
                    System._sys._forces.add(newBondForce)


            elif match.group('Angles'):
                if verbose:
                    print 'Parsing Angles...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = [int(field) for field in data_lines.pop(i).split()]

                    newAngleForce = None
                    atomtype1 = currentMolecule._atoms[fields[2]-1].bondtype
                    atomtype2 = currentMolecule._atoms[fields[3]-1].bondtype
                    atomtype3 = currentMolecule._atoms[fields[4]-1].bondtype
                    if len(angle_style) == 1:
                        if angle_style[0] == 'harmonic':
                            fields.append(angle_types[int(fields[1])][1])
                            fields.append(angle_types[int(fields[1])][0])
                            newAngleForce = Angle(fields[2],
                                    fields[3],
                                    fields[4],
                                    fields[5],
                                    fields[6])
                        currentMoleculeType.angleForceSet.add(newAngleForce)
                    System._sys._forces.add(newAngleForce)

            elif match.group('Dihedrals'):
                if verbose:
                    print 'Parsing Dihedrals...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = [int(field) for field in data_lines.pop(i).split()]

                    newDihedralForce = None
                    atomtype1 = currentMolecule._atoms[fields[2]-1].bondtype
                    atomtype2 = currentMolecule._atoms[fields[3]-1].bondtype
                    atomtype3 = currentMolecule._atoms[fields[4]-1].bondtype
                    atomtype4 = currentMolecule._atoms[fields[5]-1].bondtype
                    if len(dihedral_style) == 1:
                        if dihedral_style[0] == 'opls':
                            Fs = [dihedral_types[int(fields[1])][0],
                                  dihedral_types[int(fields[1])][1],
                                  dihedral_types[int(fields[1])][2],
                                  dihedral_types[int(fields[1])][3]]
                            Cs_temp = ConvertFromOPLSToRBDihedral(Fs[0]._value,
                                    Fs[1]._value,
                                    Fs[2]._value,
                                    Fs[3]._value)
                            Cs = [param * Fs[0].unit for param in Cs_temp]
                            newDihedralForce = RBDihedral(fields[2],
                                    fields[3],
                                    fields[4],
                                    fields[5],
                                    Cs[0],
                                    Cs[1],
                                    Cs[2],
                                    Cs[3],
                                    Cs[4],
                                    Cs[5],
                                    Cs[6])
                        currentMoleculeType.dihedralForceSet.add(newDihedralForce)
                    System._sys._forces.add(newDihedralForce)


            elif match.group('Impropers'):
                if verbose:
                    print 'Parsing Impropers...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = [int(field) for field in data_lines.pop(i).split()]

            else:
                i += 1
        else:
            i += 1


def writeData(filename, unit_set='real'):
    """Write system to LAMMPS data file

    Args:
        filename (str): name of output file
        unit_set (str): LAMMPS unit set for output file
    """

    # TODO: allow user to specify different unit sets
    RAD = units.radians
    DEGREE = units.degrees
    if unit_set == 'real':
        DIST = units.angstroms
        VEL = units.angstroms / units.femtosecond
        ENERGY = units.kilocalorie / units.mole
        MASS = units.grams / units.mole
        CHARGE = units.elementary_charge
        MOLE = units.mole
    else:
        raise Exception("Unsupported unit set specified: {0}".format(unit_set))

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


    # dicts for type information
    atom_type_dict = dict()  # str_type:int_type
    a_type_i = 1  # counter for atom types

    bond_style = []
    bond_type_dict = dict()  # typeObject:int_type
    b_type_i = 1  # counter for bond types

    angle_style = []
    angle_type_dict = dict()
    ang_type_i = 1

    dihedral_style = []
    dihedral_type_dict = dict()
    dih_type_i = 1

    # read all atom specific and FF information
    for moleculeType in System._sys._molecules.itervalues():
        # bond types
        if moleculeType.bondForceSet:
            for bond in moleculeType.bondForceSet.itervalues():
                if isinstance(bond, Bond):
                    if 'harmonic' not in bond_style:
                        bond_style.append('harmonic')
                    if len(bond_style) > 1:
                        # this may need to be an error
                        # or require some form of conversion if possible
                        warnings.warn("More than one bond style found!")

                    atom1 = moleculeType.moleculeSet[0]._atoms[bond.atom1 - 1]
                    atomtype1 = atom1.bondtype
                    atom2 = moleculeType.moleculeSet[0]._atoms[bond.atom2 - 1]
                    atomtype2 = atom2.bondtype

                    temp = BondType(atomtype1,
                            atomtype2,
                            1,
                            bond.length,
                            bond.k)
                    # NOTE: k includes the factor of 0.5 for harmonics in LAMMPS
                    #       For now, I will assume that we always store it without internally
                    #       since that's the way GROMACS does it.
                    if temp not in bond_type_dict:
                        bond_type_dict[temp] = b_type_i
                        bond_coeff_list.append('{0:d} {1:18.8f} {2:18.8f}\n'.format(
                                b_type_i,
                                0.5 * bond.k.in_units_of(ENERGY/(DIST*DIST))._value,
                                bond.length.in_units_of(DIST)._value))
                        b_type_i += 1
                else:
                    warnings.warn("Found unsupported bond type for LAMMPS!")
        # angle types
        if moleculeType.angleForceSet:
            for angle in moleculeType.angleForceSet.itervalues():
                if isinstance(angle, Angle):
                    if 'harmonic' not in angle_style:
                        angle_style.append('harmonic')
                    if len(angle_style) > 1:
                        # this may need to be an error
                        # or require some form of conversion if possible
                        warnings.warn("More than one angle style found!")

                    atom1 = moleculeType.moleculeSet[0]._atoms[angle.atom1 - 1]
                    atomtype1 = atom1.bondtype
                    atom2 = moleculeType.moleculeSet[0]._atoms[angle.atom2 - 1]
                    atomtype2 = atom2.bondtype
                    atom3 = moleculeType.moleculeSet[0]._atoms[angle.atom3 - 1]
                    atomtype3 = atom3.bondtype

                    temp = AngleType(atomtype1,
                            atomtype2,
                            atomtype3,
                            1,
                            angle.theta,
                            angle.k)
                    # NOTE: k includes the factor of 0.5 for harmonics in LAMMPS
                    #       For now, I will assume that we always store it without internally
                    #       since that's the way GROMACS does it.
                    if temp not in angle_type_dict:
                        angle_type_dict[temp] = ang_type_i
                        angle_coeff_list.append('{0:d} {1:18.8f} {2:18.8f}\n'.format(
                                ang_type_i,
                                0.5 * angle.k.in_units_of(ENERGY/(RAD*RAD))._value,
                                angle.theta.in_units_of(DEGREE)._value))
                        ang_type_i += 1
                else:
                    warnings.warn("Found unsupported angle type for LAMMPS!")

        # dihedral types
        if moleculeType.dihedralForceSet:
            for dihedral in moleculeType.dihedralForceSet.itervalues():
                if isinstance(dihedral, ProperDihedral1):
                    continue
                if isinstance(dihedral, RBDihedral):
                    if 'opls' not in dihedral_style:
                        dihedral_style.append('opls')
                    if len(dihedral_style) > 1:
                        # this may need to be an error
                        # or require some form of conversion if possible
                        warnings.warn("More than one dihedral style found!")

                    atom1 = moleculeType.moleculeSet[0]._atoms[dihedral.atom1 - 1]
                    atomtype1 = atom1.bondtype
                    atom2 = moleculeType.moleculeSet[0]._atoms[dihedral.atom2 - 1]
                    atomtype2 = atom2.bondtype
                    atom3 = moleculeType.moleculeSet[0]._atoms[dihedral.atom3 - 1]
                    atomtype3 = atom3.bondtype
                    atom4 = moleculeType.moleculeSet[0]._atoms[dihedral.atom4 - 1]
                    atomtype4 = atom4.bondtype

                    temp = RBDihedralType(atomtype1,
                            atomtype2,
                            atomtype3,
                            atomtype4,
                            3,
                            dihedral.C0,
                            dihedral.C1,
                            dihedral.C2,
                            dihedral.C3,
                            dihedral.C4,
                            dihedral.C5,
                            dihedral.C6)
                    if temp not in dihedral_type_dict:
                        Fs_temp = ConvertFromRBToOPLSDihedral(dihedral.C0._value,
                                dihedral.C1._value,
                                dihedral.C2._value,
                                dihedral.C3._value,
                                dihedral.C4._value,
                                dihedral.C5._value,
                                dihedral.C6._value)
                        Fs = [param * dihedral.C0.unit for param in Fs_temp]
                        dihedral_type_dict[temp] = dih_type_i
                        dihedral_coeff_list.append('{0:d} {1:18.8f} {2:18.8f} {3:18.8f} {4:18.8f}\n'.format(
                                dih_type_i,
                                Fs[0].in_units_of(ENERGY)._value,
                                Fs[1].in_units_of(ENERGY)._value,
                                Fs[2].in_units_of(ENERGY)._value,
                                Fs[3].in_units_of(ENERGY)._value))
                        dih_type_i += 1
                else:
                    warnings.warn("Found unsupported dihedral type for LAMMPS!")

        # atom specific information
        x_min = y_min = z_min = np.inf
        for molecule in moleculeType.moleculeSet:
            for atom in molecule._atoms:
                # type, mass and pair coeffs
                if atom._atomtype[0] not in atom_type_dict:
                    atom_type_dict[atom._atomtype[0]] = a_type_i
                    mass_list.append('%d %8.4f\n'
                                % (a_type_i,
                                   atom._mass[0].in_units_of(MASS)._value))
                    pair_coeff_list.append('{0:d} {1:8.4f} {2:8.4f}\n'.format(
                                   a_type_i,
                                   atom._epsilon[0].in_units_of(ENERGY)._value,
                                   atom._sigma[0].in_units_of(DIST)._value))
                    a_type_i += 1

                # box minima
                x = atom._position[0].in_units_of(DIST)._value
                y = atom._position[1].in_units_of(DIST)._value
                z = atom._position[2].in_units_of(DIST)._value
                if x < x_min:
                    x_min = x
                if y < y_min:
                    y_min = y
                if z < z_min:
                    z_min = z

                # atom
                atom_list.append('{0:-6d} {1:-6d} {2:-6d} {3:5.8f} {4:8.3f} {5:8.3f} {6:8.3f}\n'.format(
                        atom.atomIndex,
                        atom.residueIndex,
                        atom_type_dict[atom._atomtype[0]],
                        atom._charge[0].in_units_of(CHARGE)._value,
                        x,
                        y,
                        z))
                # velocity
                vel_list.append('{0:-6d} {1:8.4f} {2:8.4f} {3:8.4f}\n'.format(
                        atom.atomIndex,
                        atom._velocity[0].in_units_of(VEL)._value,
                        atom._velocity[1].in_units_of(VEL)._value,
                        atom._velocity[2].in_units_of(VEL)._value))

    # read all connectivity information
    for moleculeType in System._sys._molecules.itervalues():
        # atom index offsets from 1 for each molecule
        offsets = list()
        for molecule in moleculeType.moleculeSet:
            offsets.append(molecule._atoms[0].atomIndex - 1)

        for i, offset in enumerate(offsets):
            for j, bond in enumerate(moleculeType.bondForceSet.itervalues()):
                atom1 = moleculeType.moleculeSet[0]._atoms[bond.atom1 - 1]
                atomtype1 = atom1.bondtype
                atom2 = moleculeType.moleculeSet[0]._atoms[bond.atom2 - 1]
                atomtype2 = atom2.bondtype

                temp = BondType(atomtype1,
                        atomtype2,
                        1,
                        bond.length,
                        bond.k)

                bond_list.append('{0:-6d} {1:6d} {2:6d} {3:6d}\n'.format(
                        i + j + 1,
                        bond_type_dict[temp],
                        bond.atom1 + offset,
                        bond.atom2 + offset))

            for j, angle in enumerate(moleculeType.angleForceSet.itervalues()):
                atom1 = moleculeType.moleculeSet[0]._atoms[angle.atom1 - 1]
                atomtype1 = atom1.bondtype
                atom2 = moleculeType.moleculeSet[0]._atoms[angle.atom2 - 1]
                atomtype2 = atom2.bondtype
                atom3 = moleculeType.moleculeSet[0]._atoms[angle.atom3 - 1]
                atomtype3 = atom3.bondtype

                temp = AngleType(atomtype1,
                        atomtype2,
                        atomtype3,
                        1,
                        angle.theta,
                        angle.k)

                angle_list.append('{0:-6d} {1:6d} {2:6d} {3:6d} {4:6d}\n'.format(
                        i + j + 1,
                        angle_type_dict[temp],
                        angle.atom1 + offset,
                        angle.atom2 + offset,
                        angle.atom3 + offset))

            for j, dihedral in enumerate(moleculeType.dihedralForceSet.itervalues()):
                atom1 = moleculeType.moleculeSet[0]._atoms[dihedral.atom1 - 1]
                atomtype1 = atom1.bondtype
                atom2 = moleculeType.moleculeSet[0]._atoms[dihedral.atom2 - 1]
                atomtype2 = atom2.bondtype
                atom3 = moleculeType.moleculeSet[0]._atoms[dihedral.atom3 - 1]
                atomtype3 = atom3.bondtype
                atom4 = moleculeType.moleculeSet[0]._atoms[dihedral.atom4 - 1]
                atomtype4 = atom4.bondtype

                if isinstance(dihedral, ProperDihedral1):
                    continue
                temp = RBDihedralType(atomtype1,
                            atomtype2,
                            atomtype3,
                            atomtype4,
                            3,
                            dihedral.C0,
                            dihedral.C1,
                            dihedral.C2,
                            dihedral.C3,
                            dihedral.C4,
                            dihedral.C5,
                            dihedral.C6)

                dihedral_list.append('{0:-6d} {1:6d} {2:6d} {3:6d} {4:6d} {5:6d}\n'.format(
                        i + j + 1,
                        dihedral_type_dict[temp],
                        dihedral.atom1 + offset,
                        dihedral.atom2 + offset,
                        dihedral.atom3 + offset,
                        dihedral.atom4 + offset))

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
        if n_bond_types > 0:
            f.write('{0} bond types\n'.format(n_bond_types))
        if n_angle_types > 0:
            f.write('{0} angle types\n'.format(n_angle_types))
        if n_dihedral_types > 0:
            f.write('{0} dihedral types\n'.format(n_dihedral_types))
        if n_improper_types > 0:
            f.write('{0} improper types\n'.format(n_improper_types))
        f.write('\n')
        # shifting of box dimensions
        f.write('{0:10.6f} {1:10.6f} xlo xhi\n'.format(
                x_min,
                x_min + System._sys._boxVector[0][0].in_units_of(DIST)._value))
        f.write('{0:10.6f} {1:10.6f} ylo yhi\n'.format(
                y_min,
                y_min + System._sys._boxVector[1][1].in_units_of(DIST)._value))
        f.write('{0:10.6f} {1:10.6f} zlo zhi\n'.format(
                z_min,
                z_min + System._sys._boxVector[2][2].in_units_of(DIST)._value))

        # write masses
        for mass in mass_list:
            f.write(mass)

        # write coefficients
        if len(pair_coeff_list) > 3:
            for pair in pair_coeff_list:
                f.write(pair)
        if len(bond_coeff_list) > 3:
            for bond in bond_coeff_list:
                f.write(bond)
        if len(angle_coeff_list) > 3:
            for angle in angle_coeff_list:
                f.write(angle)
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
        f.write('pair_style lj/cut/coul/long 9.0 10.0\n')  # TODO: match mdp
        f.write('pair_modify mix geometric\n')  # TODO: match defaults
        f.write('kspace_style pppm 1.0e-4\n')  # TODO: match mdp
        f.write('\n')

        # bonded
        if len(bond_coeff_list) > 3:
            if len(bond_style) == 1:
                f.write('bond_style {0}\n'.format(" ".join(bond_style)))
            else:
                f.write('bond_style hybrid {0}\n'.format(" ".join(bond_style)))
        if len(angle_coeff_list) > 3:
            if len(angle_style) == 1:
                f.write('angle_style {0}\n'.format(" ".join(angle_style)))
            else:
                f.write('angle_style hybrid {0}\n'.format(" ".join(angle_style)))
        if len(dihedral_coeff_list) > 3:
            if len(dihedral_style) == 1:
                f.write('dihedral_style {0}\n'.format(" ".join(dihedral_style)))
            else:
                f.write('dihedral_style hybrid {0}\n'.format(" ".join(dihedral_style)))
        if len(improper_coeff_list) > 3:
            if len(improper_style) == 1:
                f.write('improper_style {0}\n'.format(" ".join(improper_style)))
            else:
                f.write('improper_style hybrid {0}\n'.format(" ".join(improper_style)))

        f.write('special_bonds lj {0} {1} {2} coul {3} {4} {5}\n'.format(0.0,
                0.0,
                System._sys._ljCorrection,
                0.0,
                0.0,
                System._sys._coulombCorrection))
        f.write('\n')

        # read data
        f.write('read_data {0}\n'.format(os.path.basename(filename)))
        f.write('\n')

        # output energies
        energy_terms = " ".join(['ebond',
                                 'eangle',
                                 'edihed',
                                 'eimp',
                                 'epair',
                                 'evdwl',
                                 'ecoul',
                                 'pe'])

        f.write('thermo_style custom {0}\n'.format(energy_terms))
        f.write('\n')

        f.write('run 0\n')

        # special bonds? CHARMM?
        # multiple types in one category?
