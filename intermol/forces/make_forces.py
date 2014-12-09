import os
from forcedata import *
from forcefunctions import *

"""
To add another bonded force interacton, add to forcedata.py
 in one of the existing categories:

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
for name, uset in forcedata.master_unitlist.iteritems():
    # shouldn't really matter which unitset we use for writing the function types, as they 
    # are just for testing functional compatibility
    unitset = specify(forcedata.ProgramUnitSets['gromacs'], uset, None, False)  
    all_unitlist[name] = unitset

# loop over the force lists
print "Generating the following files"

for forcelist in forcelists:

    # define the things that are different for each force type
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
        # natoms will be set inside the loop

    abstractparams = AbstractOptParams[ftype]
    abstractdefaults = AbstractOptParamsDefaults[ftype]

    # construct the abstract params, abstract params default, and abstract slots string
    # for printing from the lists
    absParamStr = ''
    absParamStrDefaults = ''
    absParamStrSelf = ''
    absSlots = ''
    for i, p in enumerate(abstractparams):
        absParamStr += p
        absParamStrDefaults += (p + ' = ' + abstractdefaults[i])
        absParamStrSelf += (p + ' = ' + p)
        absSlots += '\'' + p + '\''
        if p != abstractparams[-1]:
            absParamStr += ', '
            absParamStrDefaults += ', '
            absParamStrSelf += ', '
            absSlots += ', '

    # loop over the individual members of the force list

    for force in forcelist:

        # these depend on the number of atoms but for virtual sites, the number of atoms
        # is inside, so we put this inside

        if forcelist == virtualsitelist:
            if 'two' in force:
                natoms = 3
            if 'three' in force:
                natoms = 4
            if 'four' in force:
                natoms = 5
            if 'n_' in force:
                break   # we don't support these yet

        # The abstract type module
        if forcelist == virtualsitelist:
            abstract_type_file = 'abstract_%d_' % (natoms-1) + ftype + '_type'
        else:
            abstract_type_file = 'abstract_' + ftype + '_type'

        # The name of the abstract class used
        abstract_type = capifyname(abstract_type_file)

        bondingtypes = ''
        optbondingtypes = ''   # for when bonding types are optional, in the force definition
        atoms = ''
        for i in range(natoms):
            si = str(i+1)
            bondingtypes += 'bondingtype' + si + ', '
            optbondingtypes += 'bondingtype' + si + ' = None, '
            atoms += 'atom' + si + ', '

        # name of the force
        forcename = force + '_' + ftype

        # actual name of the file for the module
        filename = forcename + '_type.py'

        # open the file up
        f = open(filename,'w')

        # imports
        f.write('from intermol.decorators import *\n')
        f.write('from ' + abstract_type_file + ' import *\n\n')

        # name of the class in camelCase
        capname = capifyname(forcename)

        # name of the class type in camelCase
        capnametype = capname + 'Type'

        # define the class type
        f.write('class ' + capnametype + '(' + abstract_type + '):\n')

        # list the slots for the type -- depends on parameters and abstract class parameters
        f.write('    __slots__ = [')
        for param in forcedata.master_paramlist[forcename]:
            f.write('\'' + param + '\', ')
        f.write(absSlots + ']\n\n')

        # now write the accepts_compatible_units section.  Depends on the num of atoms,
        # number of parameters, and abstract class pameters

        spaces = ' ' * 30
        f.write('    @accepts_compatible_units(None,\n')
        # units for atoms (None)
        for i in range(natoms-1):
            f.write(spaces + 'None,\n')
        # write out the explicit units
        for i, p in enumerate(forcedata.master_paramlist[forcename]):
            f.write(spaces + p + ' = ' + all_unitlist[forcename][i] + ',\n')
        # write none the number of times of the abstract default parameters
        for p in abstractparams:
            f.write(spaces + p + ' = None')
            if p == abstractparams[-1]:
                f.write(')\n\n')
            else:
                f.write(',\n')

        # now write the init of the force type
        f.write('    def __init__(self, ' + bondingtypes)
        for i, param in enumerate(forcedata.master_paramlist[forcename]):
            f.write(param + ' = 0.0 * ' + all_unitlist[forcename][i] + ', ')
        f.write(absParamStrDefaults + '):\n')

        # initialize the abstract type of this force
        spaces = ' ' * 8
        f.write(spaces + abstract_type + '.__init__(self, ' + bondingtypes)
        f.write(absParamStr + ')\n')

        #initialize the specific parameters for this force
        for param in forcedata.master_paramlist[forcename]:
            f.write(spaces + 'self.' + param + ' = ' + param + '\n')
        f.write('\n')

        # write the class for the force
        f.write('class ' + capname + '(' + capnametype + '):\n')
        # include the docstring here
        f.write('    """\n')
        f.write('    ' + doclist[forcename])
        f.write('    """\n')

        # define the init for the class
        f.write('    def __init__(self, '+ atoms + optbondingtypes)
        for i, param in enumerate(forcedata.master_paramlist[forcename]):
            f.write(param + ' = 0.0 * ' + all_unitlist[forcename][i] + ', '),
        f.write(absParamStrDefaults + '):\n')

        # init the atoms
        for i in range(natoms):
            iatom = str(i+1)
            f.write('        self.atom' + iatom + ' = atom' + iatom + '\n')

        # init the type of this force
        f.write('        ' + capnametype + '.__init__(self, ' + bondingtypes)
        for param in forcedata.master_paramlist[forcename]:
            f.write(param + ' = ' + param + ', ')
        f.write(absParamStrSelf + ')\n\n')

        # define equality for the class
        # f.write('    def __eq__(self, ' + ftype + '):\n')
        #
        # #Check whether atoms are equal to each other
        # f.write('\n        equality = (\n')
        # for i in range(natoms):
        #     f.write('            self.atom%d == ' % (i+1) + ftype + '.atom%d' % (i+1))
        #     if i==natoms-1:
        #         f.write('\n            ) or (\n')
        #     else:
        #         f.write(' and\n')

        # this might be bad for virtuals - not symmetric.  Leave for now.  Can't have matching virtuals, as a
        # virtual can't appear in a virtual list.
        # for i in range(natoms):
        #     f.write('            self.atom%d == ' % (i+1) + ftype + '.atom%d' % (natoms-i))
        #     if i==natoms-1:
        #         f.write(')\n\n')
        #     else:
        #         f.write(' and\n')
        # if ftype == 'nonbonded':
        #     f.write('        equality = equality and self.type == nonbonded.type\n')
        # f.write('        return equality\n')
        # f.write('\n')

        # define the hash
        # f.write('    def __hash__(self):\n')
        # f.write('        return hash(tuple([')
        # for i in range(natoms):
        #     f.write('self.atom' + str(i+1))
        #     if (i == natoms-1):
        #         f.write(']))\n\n')
        #     else:
        #         f.write(', ')

        f.close()
        # print the filename so we know which files should not be edited directly.
        print filename
