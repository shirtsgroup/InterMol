from __future__ import print_function

from intermol.forces.forcedata import *
from intermol.forces.forcefunctions import *

"""To add another bonded force interacton, add to forcedata.py in one of the
existing categories:

1. add it to the appropriate forcelist
2  define the parameters
3. define the units for the parameters
4. write the docstrings

Adding new categories will take extra effort, adding an
abstract_force_type, and additional logic, but should be clear.

Currently, the following files are printed, so should not be edited
directly. (take the output of 'python make_forces.py' to update this
section)

force_types/lj_c_nonbonded_type.py
force_types/lj_sigeps_nonbonded_type.py
force_types/buckingham_nonbonded_type.py
force_types/lj_c_pair_type.py
force_types/lj_sigeps_pair_type.py
force_types/ljq_c_pair_type.py
force_types/ljq_sigeps_pair_type.py
force_types/lj_default_pair_type.py
force_types/ljq_default_pair_type.py
force_types/harmonic_bond_type.py
force_types/fene_bond_type.py
force_types/connection_bond_type.py
force_types/morse_bond_type.py
force_types/cubic_bond_type.py
force_types/harmonic_potential_bond_type.py
force_types/g96_bond_type.py
force_types/harmonic_angle_type.py
force_types/urey_bradley_angle_type.py
force_types/cross_bond_angle_angle_type.py
force_types/cross_bond_bond_angle_type.py
force_types/g96_angle_type.py
force_types/quartic_angle_type.py
force_types/improper_harmonic_dihedral_type.py
force_types/trig_dihedral_type.py
force_types/fourier_dihedral_type.py
force_types/proper_periodic_dihedral_type.py
force_types/rb_dihedral_type.py
force_types/two_virtual_type.py
force_types/three_linear_virtual_type.py
force_types/three_fd_virtual_type.py
force_types/three_fad_virtual_type.py
force_types/three_out_virtual_type.py
force_types/four_fdn_virtual_type.py
"""

all_unitlist = dict()
for name, uset in master_unitlist.items():
    # Shouldn't really matter which unitset we use for writing the function
    # types, as they are just for testing functional compatibility.
    unitset = specify(ProgramUnitSets['gromacs'], uset, None, False)
    all_unitlist[name] = unitset

print("Generating the following files")
for forcelist in forcelists:
    # Define the things that are different for each force type.
    if forcelist is nonbondedlist:
        ftype = "nonbonded"
        natoms = 2
    if forcelist is pairlist:
        ftype = "pair"
        natoms = 2
    if forcelist is bondlist:
        ftype = "bond"
        natoms = 2
    elif forcelist is anglelist:
        ftype = "angle"
        natoms = 3
    elif forcelist is dihedrallist:
        ftype = "dihedral"
        natoms = 4
    elif forcelist is virtualsitelist:
        ftype = "virtual"

    abstractparams = AbstractOptParams[ftype]
    abstractdefaults = AbstractOptParamsDefaults[ftype]

    # Construct the abstract params, abstract params default, and abstract slots
    # string for printing from the lists.
    absParamStr = ''
    absParamStrDefaults = ''
    absParamStrSelf = ''
    absSlots = ''
    for i, p in enumerate(abstractparams):
        absParamStr += p
        absParamStrDefaults += (p + '=' + abstractdefaults[i])
        absParamStrSelf += (p + '=' + p)
        absSlots += '\'' + p + '\''
        if p != abstractparams[-1]:
            absParamStr += ', '
            absParamStrDefaults += ', '
            absParamStrSelf += ', '
            absSlots += ', '

    for force in forcelist:
        # These depend on the number of atoms but for virtual sites, the number
        # of atoms is inside, so we put this inside.
        if forcelist == virtualsitelist:
            if 'two' in force:
                natoms = 3
            if 'three' in force:
                natoms = 4
            if 'four' in force:
                natoms = 5
            if 'n_' in force:
                break   # We don't support these yet.

        # The abstract type module.
        if forcelist == virtualsitelist:
            abstract_type_file = 'abstract_%d_' % (natoms-1) + ftype + '_type'
        else:
            abstract_type_file = 'abstract_' + ftype + '_type'

        # The name of the abstract class used
        abstract_type = capifyname(abstract_type_file)

        bondingtypes = ''
        # For when bonding types are optional, in the force definition.
        optbondingtypes = ''
        atoms = ''
        for i in range(natoms):
            si = str(i+1)
            bondingtypes += 'bondingtype' + si + ', '
            optbondingtypes += 'bondingtype' + si + '=None, '
            atoms += 'atom' + si + ', '

        forcename = force + '_' + ftype
        filename = forcename + '_type.py'

        with open(filename, 'w') as f:
            # Imports
            f.write('import parmed.unit as units\n\n')
            f.write('from intermol.decorators import accepts_compatible_units\n')
            f.write('from intermol.forces.{0} import {1}\n\n\n'.format(abstract_type_file, abstract_type))

            # Name of the class in camelCase.
            capname = capifyname(forcename)

            # Name of the class type in camelCase.
            capnametype = capname + 'Type'

            # Define the class type.
            f.write('class ' + capnametype + '(' + abstract_type + '):\n')

            # List the slots for the type .
            f.write('    __slots__ = [')
            for param in master_paramlist[forcename]:
                f.write('\'{0}\', '.format(param))
            f.write('{0}]\n\n'.format(absSlots))

            # Now write the accepts_compatible_units section.  Depends on the
            # num of atoms, number of parameters, and abstract class pameters.
            spaces = ' ' * 30
            f.write('    @accepts_compatible_units(')
            # Units for atoms (None)
            for i in range(natoms):
                f.write('None, ')
            f.write('\n')
            # Write out the explicit units.
            for i, p in enumerate(master_paramlist[forcename]):
                f.write('{0}{1}={2},\n'.format(
                        spaces, p, all_unitlist[forcename][i]))
            # Write none the number of times of the abstract default parameters.
            for p in abstractparams:
                f.write('{0}{1}=None'.format(spaces, p))
                if p == abstractparams[-1]:
                    f.write(')\n')
                else:
                    f.write(',\n')

            # Now write the init of the force type.
            f.write('    def __init__(self, {0}\n'.format(bondingtypes))
            n_params = len(master_paramlist[forcename])
            spaces = ' ' * 17
            for i, param in enumerate(master_paramlist[forcename]):
                f.write('{0}{1}=0.0 * {2}'.format(
                        spaces, param, all_unitlist[forcename][i]))
                if i == n_params:
                    f.write(')\n\n')
                else:
                    f.write(',\n')
            f.write('{0}{1}):\n'.format(spaces, absParamStrDefaults))

            # Initialize the abstract type of this force.
            spaces = ' ' * 8
            f.write('{0}{1}.__init__(self, {2}'.format(
                    spaces, abstract_type, bondingtypes))
            f.write(absParamStr + ')\n')

            # Initialize the specific parameters for this force
            for param in master_paramlist[forcename]:
                f.write('{0}self.{1} = {1}\n'.format(spaces, param))
            f.write('\n\n')

            # Write the class for the force.
            f.write('class {0}({1}):\n'.format(capname, capnametype))
            # Include the docstring here.
            f.write('    """\n')
            f.write('    ' + doclist[forcename])
            f.write('    """\n')

            # Define the init for the class.
            f.write('    def __init__(self, {0}{1}\n'.format(
                    atoms, optbondingtypes))
            spaces = ' ' * 17
            for i, param in enumerate(master_paramlist[forcename]):
                f.write('{0}{1}=0.0 * {2}'.format(spaces, param, all_unitlist[forcename][i]))
                if i == n_params:
                    f.write(')\n\n')
                else:
                    f.write(',\n')
            f.write('{0}{1}):\n'.format(spaces, absParamStrDefaults))

            # Init the atoms.
            for atom_index in range(1, natoms+1):
                f.write('        self.atom{0} = atom{0}\n'.format(atom_index))

            # Init the type of this force.
            f.write('        {0}.__init__(self, {1}\n'.format(capnametype, bondingtypes))
            spaces = ' ' * 16
            for i, param in enumerate(master_paramlist[forcename]):
                f.write('{0}{1}={1},\n'.format(spaces, param))
            f.write('{0}{1})'.format(spaces, absParamStrSelf))

            # Print which files should not be edited directly.
            print(filename)
