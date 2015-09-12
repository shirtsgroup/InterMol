from collections import OrderedDict
import logging
import os
import subprocess

import simtk.unit as units

import intermol.tests
from intermol.tests.testing_tools import which
import shutil

logger = logging.getLogger('InterMolLog')


# --------- energy evaluation methods ---------- #
key_dict = {'totE': 'Potential',
            'bond': 'Bond',
            'angle': 'Angle',
            'dihed': 'All dihedrals',
            '1-4 vdW': 'LJ-14',
            '1-4 ee': 'Coulomb-14',
            'ee': 'Coulomb', 
            'vdW': 'LJ'
            }

def standardize_key(in_key):
    if in_key in key_dict:
        out_key = key_dict[in_key]
    else:
        out_key = in_key
    return out_key

def amber_energies(ligandloc):

    """
    Evalutes energies of AMBER ligands (given the dictionary)

    #repetitive to keep opening the same file, but better than muddling convert.py
    """

#    f = open('/Users/mrshirts/work/SAMPL5/new_simulation_inputs/simulation_inputs/Amber_energies.csv','r')
    f = open('/Users/mrshirts/work/SAMPL5/simulation_inputs/Amber_energies.csv','r')
    allline = f.readlines()
    if len(allline) == 1:
        lines = allline[0].split('\r')
    ligands = dict()
    
    ligand = ligandloc.split('/')[-2]
    for line in lines:
        if line[0:4] == 'name':
            things = line.split(',')
            enames = things[1:]
        else:
            things = line.split(',')
            if ligand == things[0]:
                energies = OrderedDict()
                for i,ename in enumerate(enames):
                    if ename in key_dict.keys():
                        gname = key_dict[ename]
                        energies[gname] = float(things[i+1])*units.kilocalories_per_mole

    logger.info('Extracting energy of {0}'.format(ligand))

    return energies
