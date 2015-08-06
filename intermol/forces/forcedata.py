from __future__ import print_function

"""
To add another bonded force interacton in one of the existing categories:

1. Add the force name to the appropriate forcelist

2. Define the parameters.  Order should reflect the most common order.
   Different orderers are specified by providing a parameterlist with
   the same terms in a different order.

3. Define the units for the parameters.  Ordering should match the
   master parameter list.  It is automatically reshuffled to match the
   program specific parameter list, so only one master unitlist is
   required.

4. write the docstrings

5. add the entry in __init__.py

Adding new categories will take extra effort, adding an
abstract_force_type, and additional logic, but should be clear.
"""

programs = ['desmond', 'gromacs', 'lammps']
ProgramUnitSets = dict()
ProgramUnitSets['desmond'] = {'angleD': 'units.degrees',
                              'angleR': 'units.radians',
                              'amount': 'units.mole',
                              'charge': 'units.elementary_charge',
                              'energy': 'units.kilocalories_per_mole',
                              'length': 'units.angstroms',
                              'mass': 'units.amu',
                              'time': 'units.picoseconds',
                              'temperature': 'units.kelvin'
                              }
ProgramUnitSets['gromacs'] = {'angleD': 'units.degrees',
                              'angleR': 'units.radians',
                              'amount': 'units.mole',
                              'charge': 'units.elementary_charge',
                              'energy': 'units.kilojoules_per_mole',
                              'length': 'units.nanometers',
                              'mass': 'units.amu',
                              'time': 'units.picoseconds',
                              'temperature': 'units.kelvin'
                              }
ProgramUnitSets['lammps'] = {'angleD': 'dumself.DEGREE',
                             'angleR': 'dumself.RAD',
                             'amount': 'dumself.MOLE',
                             'charge': 'dumself.CHARGE',
                             'energy': 'dumself.ENERGY',
                             'length': 'dumself.DIST',
                             'mass': 'dumself.MASS',
                             'time': 'dumself.DIST / dumself.VEL',
                             'temperature': 'dumself.TEMP'
                             }

AbstractOptParams = dict()
AbstractOptParamsDefaults = dict()

# defined by abstract_nonbonded_type
nonbondedlist = ['lj_c',
                 'lj_sigeps',
                 'buckingham'
                 ]
AbstractOptParams['nonbonded'] = ['type']
AbstractOptParamsDefaults['nonbonded'] = ['False']

# defined by abstract_pair_type
pairlist = ['lj_c',
            'lj_sigeps',
            'ljq_c',
            'ljq_sigeps',
            'lj_default',
            'ljq_default'
            ]
AbstractOptParams['pair'] = ['scaleLJ', 'scaleQQ', 'long']
AbstractOptParamsDefaults['pair'] = ['None', 'None', 'False']

# defined by abstract_bond_type
bondlist = ['harmonic',
            'fene',
            'fene_expandable',
            'connection',
            'morse',
            'nonlinear',
            'quartic',
            'quartic_breakable',
            'cubic',
            'harmonic_potential',
            'g96'
            ]
AbstractOptParams['bond'] = ['order', 'c']
AbstractOptParamsDefaults['bond'] = ['1', 'False']

# defined by abstract_angle_type
anglelist = ['harmonic',
             'urey_bradley',
             'urey_bradley_noharm',
             'cross_bond_angle',
             'cross_bond_bond',
             'cosine',
             'cosine_squared',
             'quartic',
             'restricted_bending'
             ]
AbstractOptParams['angle'] = ['c']  # constrained
AbstractOptParamsDefaults['angle'] = ['False']

# defined by abstract_dihedral_type
dihedrallist = ['improper_harmonic',
                'trig',
                'fourier',
                'proper_periodic',
                'rb',
                'restricted_bending',
                'bending_torsion'
                ]
AbstractOptParams['dihedral'] = ['improper']
AbstractOptParamsDefaults['dihedral'] = ['False']

virtualsitelist = ['two',
                   'three_linear',
                   'three_fd',
                   'three_fad',
                   'three_out',
                   'four_fdn',
                   'n_cog',
                   'n_com',
                   'n_cow'
                   ]
AbstractOptParams['virtual'] = ['placeholder']
AbstractOptParamsDefaults['virtual'] = ['False']

forcelists = [nonbondedlist, pairlist, bondlist, anglelist, dihedrallist,
              virtualsitelist]

master_unitlist = dict()
master_paramlist = dict()

# Include the differences in ordering from the master paramlist.
# Gromacs is the default ordering for convenience when it also
# supports.

desmond_paramlist = dict()
lammps_paramlist = dict()
gromacs_paramlist = dict()

# not using these yet
desmond_unitlist = dict()
lammps_unitlist = dict()
gromacs_unitlist = dict()

doclist = dict()

# =========
# nonbonded
# =========
doclist['lj_c_nonbonded'] = 'stub documentation\n'
master_paramlist['lj_c_nonbonded'] = ['C6', 'C12']
master_unitlist['lj_c_nonbonded'] = ['energy * length ** (6)',
                                     'energy * length ** (12)'
]

doclist['lj_sigeps_nonbonded'] = 'stub documentation\n'
master_paramlist['lj_sigeps_nonbonded'] = ['sigma', 'epsilon']
master_unitlist['lj_sigeps_nonbonded'] = ['length',
                                          'energy'
]

doclist['buckingham_nonbonded'] = 'stub documentation\n'
master_paramlist['buckingham_nonbonded'] = ['a', 'b', 'C6']
master_unitlist['buckingham_nonbonded'] = ['energy',
                                           'length ** (-1)',
                                           'energy * length ** (6)'
]

# =====
# pairs
# =====
doclist['lj_c_pair'] = 'stub documentation\n'
master_paramlist['lj_c_pair'] = ['C6', 'C12']
master_unitlist['lj_c_pair'] = ['energy * length ** (6)',
                                'energy * length ** (12)'
]

doclist['lj_sigeps_pair'] = 'stub documentation\n'
master_paramlist['lj_sigeps_pair'] = ['sigma', 'epsilon']
master_unitlist['lj_sigeps_pair'] = ['length',
                                     'energy'
]

doclist['ljq_c_pair'] = 'stub documentation\n'
master_paramlist['ljq_c_pair'] = ['qi', 'qj', 'C6', 'C12']
master_unitlist['ljq_c_pair'] = ['charge',
                                 'charge',
                                 'energy * length ** (6)',
                                 'energy * length ** (12)'
]

doclist['ljq_sigeps_pair'] = 'stub documentation\n'
master_paramlist['ljq_sigeps_pair'] = ['qi', 'qj', 'sigma', 'epsilon']
master_unitlist['ljq_sigeps_pair'] = ['charge',
                                      'charge',
                                      'length',
                                      'energy'
]

doclist['lj_default_pair'] = 'stub documentation\n'
master_paramlist['lj_default_pair'] = []
master_unitlist['lj_default_pair'] = []

doclist['ljq_default_pair'] = 'stub documentation\n'
master_paramlist['ljq_default_pair'] = []
master_unitlist['ljq_default_pair'] = []

# ====
# bonds
# ====
doclist['connection_bond'] = 'stub documentation\n'
master_paramlist['connection_bond'] = []
master_unitlist['connection_bond'] = []

doclist['cubic_bond'] = 'stub documentation\n'
master_paramlist['cubic_bond'] = ['length', 'C2', 'C3']
master_unitlist['cubic_bond'] = ['length',
                                 'energy * length ** (-2)',
                                 'energy * length ** (-3)'
]

doclist['quartic_bond'] = 'stub documentation\n'
master_paramlist['quartic_bond'] = ['length', 'C2', 'C3', 'C4']
master_unitlist['quartic_bond'] = ['length',
                                   'energy * length ** (-2)',
                                   'energy * length ** (-3)',
                                   'energy * length ** (-4)'
]

doclist[
    'quartic_breakable_bond'] = 'http://lammps.sandia.gov/doc/bond_quartic.html\n'
master_paramlist['quartic_breakable_bond'] = ['k', 'B1', 'B2', 'Rc', 'U0']
master_unitlist['quartic_breakable_bond'] = ['energy * length ** (-4)',
                                             'length',
                                             'length',
                                             'length',
                                             'energy'
]

doclist['nonlinear_bond'] = 'http://lammps.sandia.gov/doc/bond_nonlinear.html\n'
master_paramlist['nonlinear_bond'] = ['epsilon', 'r0',
                                      'lamda']  # mispell lambda since 'lambda' is a python reserved word.
master_unitlist['nonlinear_bond'] = ['energy',
                                     'length',
                                     'length'
]

doclist['fene_bond'] = 'stub documentation\n'
master_paramlist['fene_bond'] = ['length', 'kb']
master_unitlist['fene_bond'] = ['length',
                                'energy * length ** (-2)'
]

doclist['fene_expandable_bond'] = 'stub documentation\n'
master_paramlist['fene_expandable_bond'] = ['k', 'length', 'epsilon', 'sigma',
                                            'delta']
master_unitlist['fene_expandable_bond'] = ['energy * length ** (-2)', 'length',
                                           'energy', 'length', 'length']

doclist['g96_bond'] = 'stub documentation\n'
master_paramlist['g96_bond'] = ['length', 'k']
master_unitlist['g96_bond'] = ['length',
                               'energy * length ** (-4)'
]

doclist['harmonic_bond'] = 'stub documentation\n'
master_paramlist['harmonic_bond'] = ['length', 'k']
master_unitlist['harmonic_bond'] = ['length',
                                    'energy * length ** (-2)'
]
lammps_paramlist['harmonic_bond'] = ['k', 'length']

doclist['harmonic_potential_bond'] = 'stub documentation\n'
master_paramlist['harmonic_potential_bond'] = ['length', 'k']
master_unitlist['harmonic_potential_bond'] = ['length',
                                              'energy * length ** (-2)'
]

doclist['morse_bond'] = 'stub documentation\n'
master_paramlist['morse_bond'] = ['length', 'D', 'beta']
master_unitlist['morse_bond'] = ['length',
                                 'energy',
                                 'length ** (-1)']
lammps_paramlist['morse_bond'] = ['D', 'beta', 'length']

# =====
# angles
# =====
doclist['urey_bradley_angle'] = 'stub documentation\n'
master_paramlist['urey_bradley_angle'] = ['theta', 'k', 'r', 'kUB']
master_unitlist['urey_bradley_angle'] = ['angleD',
                                         'energy * angleR ** (-2)',
                                         'length',
                                         'energy * length ** (-2)'
]
lammps_paramlist['urey_bradley_angle'] = ['k', 'theta', 'kUB', 'r']

# urey_bradley with no harmonic terms (in Desmond)
doclist['urey_bradley_noharm_angle'] = 'stub documentation\n'
master_paramlist['urey_bradley_noharm_angle'] = ['r', 'kUB']
master_unitlist['urey_bradley_noharm_angle'] = ['length',
                                                'energy * length ** (-2)'
]

doclist['harmonic_angle'] = 'stub documentation\n'
master_paramlist['harmonic_angle'] = ['theta', 'k']
master_unitlist['harmonic_angle'] = ['angleD',
                                     'energy * angleR **(-2)'
]
lammps_paramlist['harmonic_angle'] = ['k', 'theta']  # reversed

doclist['cross_bond_bond_angle'] = 'stub documentation\n'
master_paramlist['cross_bond_bond_angle'] = ['r1', 'r2', 'k']
master_unitlist['cross_bond_bond_angle'] = ['length',
                                            'length',
                                            'energy * length ** (-2)'
]

doclist['cross_bond_angle_angle'] = 'stub documentation\n'
master_paramlist['cross_bond_angle_angle'] = ['r1', 'r2', 'r3', 'k']
master_unitlist['cross_bond_angle_angle'] = ['length',
                                             'length',
                                             'length',
                                             'energy * length ** (-2)'
]

doclist['cosine_angle'] = 'http://lammps.sandia.gov/doc/angle_cosine.html\n'
master_paramlist['cosine_angle'] = ['k']
master_unitlist['cosine_angle'] = ['energy'
]

doclist['cosine_squared_angle'] = 'stub documentation\n'
master_paramlist['cosine_squared_angle'] = ['theta', 'k']
master_unitlist['cosine_squared_angle'] = ['angleD',
                                           'energy'
]
lammps_paramlist['cosine_squared_angle'] = ['k', 'theta'] # reveresd

doclist['quartic_angle'] = 'stub documentation\n'
master_paramlist['quartic_angle'] = ['theta', 'C0', 'C1', 'C2', 'C3', 'C4']
master_unitlist['quartic_angle'] = ['angleD',
                                    'energy',
                                    'energy * angleR ** (-1)',
                                    'energy * angleR ** (-2)',
                                    'energy * angleR ** (-3)',
                                    'energy * angleR ** (-4)'
]

doclist['restricted_bending_angle'] = 'stub documentation\n'
master_paramlist['restricted_bending_angle'] = ['theta', 'k']
master_unitlist['restricted_bending_angle'] = ['angleD','energy']

# ========
# dihedrals
# ========
doclist['improper_harmonic_dihedral'] = 'stub documentation\n'
master_paramlist['improper_harmonic_dihedral'] = ['xi', 'k']
master_unitlist['improper_harmonic_dihedral'] = ['angleD',
                                                 'energy * angleR **(-2)'
]
lammps_paramlist['improper_harmonic_dihedral'] = ['k', 'xi']

doclist['fourier_dihedral'] = 'stub documentation\n'
master_paramlist['fourier_dihedral'] = ['c1', 'c2', 'c3', 'c4', 'c5']
master_unitlist['fourier_dihedral'] = ['energy',
                                       'energy',
                                       'energy',
                                       'energy',
                                       'energy'
]

doclist['trig_dihedral'] = 'stub documentation\n'
master_paramlist['trig_dihedral'] = ['phi', 'fc0', 'fc1', 'fc2', 'fc3', 'fc4',
                                     'fc5', 'fc6']
master_unitlist['trig_dihedral'] = ['angleD',
                                    'energy',
                                    'energy',
                                    'energy',
                                    'energy',
                                    'energy',
                                    'energy',
                                    'energy'
]

doclist['proper_periodic_dihedral'] = 'stub documentation\n'
master_paramlist['proper_periodic_dihedral'] = ['phi', 'k', 'multiplicity',
                                                'weight']
master_unitlist['proper_periodic_dihedral'] = ['angleD',
                                               'energy',
                                               'units.dimensionless',
                                               'units.dimensionless'
]
lammps_paramlist['proper_periodic_dihedral'] = ['k', 'multiplicity', 'phi',
                                                'weight']
gromacs_paramlist['proper_periodic_dihedral'] = ['phi', 'k', 'multiplicity']
desmond_paramlist['proper_periodic_dihedral'] = ['phi', 'k', 'multiplicity']

doclist['rb_dihedral'] = 'stub documentation\n'
master_paramlist['rb_dihedral'] = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6']
# Gromacs has one less.
gromacs_paramlist['rb_dihedral'] = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
# Lammps has two less, called multi/harmonic.
lammps_paramlist['rb_dihedral'] = ['C0', 'C1', 'C2', 'C3', 'C4']

master_unitlist['rb_dihedral'] = ['energy',
                                  'energy',
                                  'energy',
                                  'energy',
                                  'energy',
                                  'energy',
                                  'energy']


doclist['restricted_bending_dihedral'] = 'stub documentation\n'
master_paramlist['restricted_bending_dihedral'] = ['theta', 'k']
master_unitlist['restricted_bending_dihedral'] = ['angleD',
                                                  'energy']


doclist['bending_torsion_dihedral'] = 'stub documentation\n'
master_paramlist['bending_torsion_dihedral'] = ['a0', 'a1','a2','a3','a4']
master_unitlist['bending_torsion_dihedral'] = ['energy',
                                               'energy',
                                               'energy',
                                               'energy',
                                               'energy']

doclist['two_virtual'] = 'stub documentation\n'
master_paramlist['two_virtual'] = ['a']
master_unitlist['two_virtual'] = ['units.dimensionless']

doclist['three_linear_virtual'] = 'stub documentation\n'
master_paramlist['three_linear_virtual'] = ['a', 'b']
master_unitlist['three_linear_virtual'] = ['units.dimensionless',
                                           'units.dimensionless']

doclist['three_fd_virtual'] = 'stub documentation\n'
master_paramlist['three_fd_virtual'] = ['a', 'd']
master_unitlist['three_fd_virtual'] = ['units.dimensionless', 'length']

doclist['three_fad_virtual'] = 'stub documentation\n'
master_paramlist['three_fad_virtual'] = ['theta', 'd']
master_unitlist['three_fad_virtual'] = ['angleD', 'length']

doclist['three_out_virtual'] = 'stub documentation\n'
master_paramlist['three_out_virtual'] = ['a', 'b', 'c']
master_unitlist['three_out_virtual'] = ['units.dimensionless',
                                        'units.dimensionless', 'length ** (-1)']

doclist['four_fdn_virtual'] = 'stub documentation\n'
master_paramlist['four_fdn_virtual'] = ['a', 'b', 'c']
master_unitlist['four_fdn_virtual'] = ['units.dimensionless',
                                       'units.dimensionless', 'length']

doclist['n_cog_virtual'] = 'stub documentation\n'
master_paramlist['n_cog_virtual'] = []
master_unitlist['n_cog_virtual'] = []

doclist['n_cow_virtual'] = 'stub documentation\n'
master_paramlist['n_cow_virtual'] = []
master_unitlist['n_cow_virtual'] = []

doclist['n_cow_virtual'] = 'stub documentation\n'
master_paramlist['n_cow_virtual'] = []
master_unitlist['n_cow_virtual'] = []

ProgramUnitLists = {'lammps': lammps_unitlist,
                    'gromacs': gromacs_unitlist,
                    'desmond': desmond_unitlist}
