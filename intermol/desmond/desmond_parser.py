import logging
from warnings import warn
import math
import numpy as np
import shlex

import parmed.unit as units
from intermol.atom import Atom
from intermol.forces import *
import intermol.forces.forcefunctions as ff
from intermol.exceptions import (UnimplementedFunctional, UnsupportedFunctional,
                                 UnimplementedSetting, UnsupportedSetting,
                                 DesmondError, InterMolError)
from intermol.molecule import Molecule
from intermol.moleculetype import MoleculeType
from intermol.system import System
from intermol.desmond import cmap_parameters


#MRS for old desmond functionality
import re
import copy 

logger = logging.getLogger('InterMolLog')

ENGINE = 'desmond'


# driver helper functions
def load(cms_file):
    """Load a DESMOND input file into a 'System'

    Args:
        cms_file:
        include_dir:
    Returns:
        system:
    """
    parser = DesmondParser(cms_file)
    return parser.read()


def save(cms_file, system):
    """Unpacks a 'System' into a DESMOND input file

    Args:
        cms_file:
        system:
    """
    parser = DesmondParser(cms_file, system)
    return parser.write()


# parser helper functions

def end_header_section(blank_section, header, header_lines):

    if blank_section:
        header_lines = list()
        header_lines.append(header)
        header_lines.append('      :::\n')
    else:
        header_lines[0] = header
    return header_lines


def split_with_quotes(line):

    if '"' in line:
        elements = shlex.split(line)
        for e in elements:
            e.replace(' ','_')
    else:
        elements = line.split()
    return elements

def create_lookup(forward_dict):
    return dict((v, k) for k, v in forward_dict.items())

def create_type(forward_dict):
    return dict((k, eval(v.__name__ + 'Type')) for k, v in forward_dict.items())

class DesmondParser(object):
    """
    A class containing methods required to read in a Desmond CMS File
    """

    # 'lookup_*' is the inverse dictionary typically used for writing
    desmond_combination_rules = {'1': 'Multiply-C6C12',
                                 '2': 'Lorentz-Berthelot',
                                 '3': 'Multiply-Sigeps'
                                 }
    lookup_desmond_combination_rules = create_lookup(desmond_combination_rules)

    desmond_pairs = {'LJ12_6_SIG_EPSILON': LjSigepsPair,
                     'LJ': LjDefaultPair,
                     'COULOMB': LjDefaultPair
                     }
    lookup_desmond_pairs = create_lookup(desmond_pairs)  # not unique
    desmond_pair_types = create_type(desmond_pairs)
     
    desmond_bonds = {'HARM_CONSTRAINED': HarmonicBond,
                     'HARM': HarmonicBond
                     }

    lookup_desmond_bonds = create_lookup(desmond_bonds)  # not unique - revisit.
    desmond_bond_types = create_type(desmond_bonds)

    def canonical_bond(self, bond, params, direction='into', name=None):

        if direction == 'into':
            canonical_force_scale = self.canonical_force_scale_into
            phase = 'Read'
        else:
            try:
                name = self.lookup_desmond_bonds[bond.__class__]  # check to make sure this OK given the c
            except:
                raise UnsupportedFunctional(bond, ENGINE)

            canonical_force_scale = self.canonical_force_scale_from
            phase = 'Write'

        names = []
        paramlists = []

        if bond.__class__ in [HarmonicBond, HarmonicPotentialBond]:

            if direction == 'into':
                bond.k *= canonical_force_scale
                if name == 'HARM_CONSTRAINED':
                    bond.c = True
                elif name == 'HARM':
                    bond.c = False
                else:
                    warn("ReadError: Found unsupported bond in Desmond {:s}".format(name))
                return bond

            else:
                params['k'] *= canonical_force_scale
                # harmonic potentials in Gromacs should be constrained (??: check what this means)
                name = 'HARM'
                if hasattr(bond,'c'):
                    if getattr(bond,'c') and not isinstance(bond, HarmonicPotentialBond):
                        name = 'HARM_CONSTRAINED'

                names.append(name)
                paramlists.append(params)
            return names, paramlists

    desmond_angles = {'HARM_CONSTRAINED': HarmonicAngle,
                      'HARM': HarmonicAngle,
                      'UB': UreyBradleyNoharmAngle
                      }

    lookup_desmond_angles = create_lookup(desmond_angles)
    desmond_angle_types = create_type(desmond_angles)

    def canonical_angle(self, angle, params, direction='into', name=None,
                        molecule_type=None):
        """
        Args:
            name:
            kwds:
            angle:
            direction: 'into' means into the canonical form, 'from' means from the
                        canonical form into Desmond
            current molecule type (would like to be able to get rid of this, but need it to search angles for now

        Returns:
            modified list of keywords and names
        """

        if direction == 'into':
            canonical_force_scale = self.canonical_force_scale_into
        else:
            # we'd like to automate this, but currently have to state explicitly.
            if angle.__class__ not in [HarmonicAngle, UreyBradleyAngle]:
               raise UnsupportedFunctional(angle, ENGINE)

            canonical_force_scale = self.canonical_force_scale_from
            phase = 'Write'

        names = []
        paramlists = []

        if angle.__class__ in [HarmonicAngle, UreyBradleyAngle, UreyBradleyNoharmAngle]:

            if direction == 'into':
                if angle.__class__ in [UreyBradleyAngle, UreyBradleyNoharmAngle]:
                    angle.kUB *=  canonical_force_scale
                if angle.__class__ in [UreyBradleyAngle, HarmonicAngle]:
                    angle.k *=  canonical_force_scale

                if name == 'HARM_CONSTRAINED':  # this needs to go first because HARM is a substring
                    angle.c = True
                elif name == 'HARM':
                    angle.c = False
            else:
                params['k'] = canonical_force_scale * params['k']
                name = 'HARM'
                if hasattr(angle,'c'):
                    if getattr(angle,'c'):
                        name = 'HARM_CONSTRAINED'

            if direction == 'into' and angle.__class__ in [UreyBradleyNoharmAngle,HarmonicAngle]:
                if angle.__class__ == UreyBradleyNoharmAngle:
                    # Urey-Bradley is implemented in DESMOND differently, with the
                    # terms implemented in a new angle term independent of the harmonic term.
                    # Instead, we will add everything together afterwards into a single term
                    angle = self.create_forcetype(UreyBradleyAngle,[angle.atom1,angle.atom2,angle.atom3],
                                                  [0,0,angle.r._value,angle.kUB._value]) # this seems kludgy
                    # next, find if we already have this angle somewhere
                    matched_angle = molecule_type.match_angles(angle)
                    if matched_angle:   # we found one, if false, we haven't seen it yet, we'll add later
                        if matched_angle.__class__ == HarmonicAngle:
                            angle.k = matched_angle.k
                            angle.theta = matched_angle.theta
                            molecule_type.angle_forces.remove(matched_angle)
                elif angle.__class__ == HarmonicAngle:
                    matched_angle = molecule_type.match_angles(angle)
                    if matched_angle and matched_angle.__class__ == UreyBradleyAngle:
                        # just copy over the information into the old angle.
                        matched_angle.k = angle.k
                        matched_angle.theta = angle.theta
                        angle = None

            elif direction == 'from' and angle.__class__ in [UreyBradleyAngle]:
                params_harmpart = {k:v for (k,v) in params.items() if k in ['theta','k','c'] }
                names.append(name)
                paramlists.append(params_harmpart)
                name = 'UB'
                params['kUB'] *= canonical_force_scale
                params_ubpart = {k:v for (k,v) in params.items() if k in ['r','kUB'] }
                names.append(name)
                paramlists.append(params_ubpart)
            else:
                if direction == 'from':
                    names.append(name)
                    paramlists.append(params)

            if direction == 'into':
                return angle
            elif direction == 'from':
                return names, paramlists

        else:
            raise UnsupportedFunctional(angle, ENGINE)

    desmond_dihedrals = {'IMPROPER_HARM': ImproperHarmonicDihedral,
                         'PROPER_TRIG': TrigDihedral,
                         'IMPROPER_TRIG': TrigDihedral,
                         'OPLS_PROPER': TrigDihedral,
                         'OPLS_IMPROPER': TrigDihedral
                         }
    lookup_desmond_dihedrals = {TrigDihedral: 'PROPER_TRIG',
                                ImproperHarmonicDihedral: 'IMPROPER_HARM'
                                }

    lookup_desmond_dihedral = create_lookup(desmond_dihedrals)
    desmond_dihedral_types = create_type(desmond_dihedrals)

    def canonical_dihedral(self, dihedral, params, direction = 'into', name = None, molecule_type = None):

        if direction == 'into':
            canonical_force_scale = self.canonical_force_scale_into
            phase = 'Read'
        else:
            try:
                name = self.lookup_desmond_dihedrals[dihedral.__class__]
            except:
                raise UnsupportedFunctional(dihedral, ENGINE)

            canonical_force_scale = self.canonical_force_scale_from
            phase = 'Write'

        if dihedral.__class__ in [ImproperHarmonicDihedral, TrigDihedral]:

            if direction == 'into':
            #Improper Diehdral 2 ---NOT SURE ABOUT MULTIPLICITY

                if name == "IMPROPER_HARM":
                    dihedral.improper = True
                elif name == "PROPER_TRIG" or name == "IMPROPER_TRIG":
                    if name == "IMPROPER_TRIG":
                        dihedral.improper = True
                    else:
                        dihedral.improper = False
                elif name == "OPLS_PROPER" or name == "OPLS_IMPROPER":
                # OPLS_IMPROPER actually isn't any different from OPLS_PROPER
                    dihedral.improper = False
                try:
                    # we can have multiple parameters with DESMOND, and append if we do
                    dihedralmatch = molecule_type.match_dihedrals(dihedral)
                    # this will fail if it's the wrong type of dihedral
                    if dihedralmatch:
                        dihedralmatch.sum_parameters(dihedral)
                except Exception as e:
                    logger.exception(e)
                return dihedral

            else:
                names = []
                paramlists = []

                if dihedral.__class__ in [ImproperHarmonicDihedral]:
                    params['k'] = params['k'] * canonical_force_scale
                    name = 'IMPROPER_HARM'

                elif dihedral.__class__ in [TrigDihedral]:
                    name = 'PROPER_TRIG'
                    if hasattr(dihedral,'improper'):
                        if getattr(dihedral,'improper'):
                            name = 'IMPROPER_TRIG'

                names.append(name)
                paramlists.append(params)

            return names, paramlists

    def __init__(self, cms_file, system=None):
        """
        Initializes a DesmondParse object which serves to read in a CMS file
        into the abstract representation.

        Args:

        """

        self.cms_file = cms_file
        if not system:
            system = System()
        self.system = system

        self.vdwtypes = []
        self.vdwtypeskeys = []

        self.viparr = 1
        self.fmct_blockpos = []
        self.atom_blockpos = []
        self.bond_blockpos = []
        self.ffio_blockpos = []

        self.paramlist = ff.build_paramlist('desmond')
        self.unitvars = ff.build_unitvars('desmond', self.paramlist)

        self.canonical_force_scale_into = 2.0
        self.canonical_force_scale_from = 0.5

        self.atom_col_vars = ['i_m_mmod_type',
                              'r_m_x_coord',
                              'r_m_y_coord',
                              'r_m_z_coord',
                              'i_m_residue_number',
                              's_m_pdb_residue_name',
                              'i_m_atomic_number',
                              's_m_pdb_atom_name',
                              's_m_atom_name',
                              'r_ffio_x_vel',
                              'r_ffio_y_vel',
                              'r_ffio_z_vel'
                              ]

        self.atom_box_vars = ['r_chorus_box_ax',
                              'r_chorus_box_ay',
                              'r_chorus_box_az',
                              'r_chorus_box_bx',
                              'r_chorus_box_by',
                              'r_chorus_box_bz',
                              'r_chorus_box_cx',
                              'r_chorus_box_cy',
                              'r_chorus_box_cz'
                              ]

    def get_parameter_list_from_kwds(self, force, kwds):
        return ff.get_parameter_list_from_kwds(force, kwds, self.paramlist)

    def get_parameter_list_from_force(self, force):
        return ff.get_parameter_list_from_force(force, self.paramlist)

    def get_parameter_kwds_from_force(self, force):
        return ff.get_parameter_kwds_from_force(force, self.get_parameter_list_from_force, self.paramlist)

    def create_kwd_dict(self, forcetype_object, values, optvalues = None):
        kwd = ff.create_kwd_dict(self.unitvars, self.paramlist, forcetype_object, values, optvalues = optvalues)
        return kwd

    def create_forcetype(self, forcetype_object, paramlist, values, optvalues = None):
        return forcetype_object(*paramlist, **self.create_kwd_dict(forcetype_object, values, optvalues))

    def parse_ffio_block(self,start,end):
    #LOAD FFIO BLOCKS IN FIRST (CONTAINS TOPOLOGY)

            # read in a ffio_block that isn't ffio_ff and split it into the
            # commands and the values.
            # lots of room for additional error checking here, such as whether
            # each entry has the correct number of data values, whether they are the correct type, etc.

            # scroll to the next ffio entry
        while not 'ffio_' in self.lines[start]:
            # this is not an ffio block! or, we have reached the end of the file
            if 'ffio_' not in self.lines[start]:
                start += 1
            if start >= end:
                return 'Done with ffio', 0, None, None, None, start

        components = re.split('\W', self.lines[start].split()[0]) # get rid of whitespace, split on nonword
        ff_type = components[0]
        ff_number = int(components[1])
        i = start+1
        entry_data = []
        while not ':::' in self.lines[i]:
            entry_data.append(self.lines[i].split()[0])
            i+=1
        i+=1 # skip the separator we just found
        entry_values = []
        while not ':::' in self.lines[i]:
            if self.lines[i].strip():  # skip the blank spaces.
                entry_values.append(self.lines[i])
            i+=1
        while '}' not in self.lines[i]:  # wait until we hit an end to the block
            i+=1
        i+=1 # step past the end of the block

        entry_dict = dict()
        for j, d in enumerate(entry_data):
            entry_dict[d] = j+1  # the first one is the entry number

        return ff_type, ff_number, entry_data, entry_values, entry_dict, i

    def store_ffio_data(self, ff_type, ff_number, entry_data, entry_values, entry_dict):

        self.stored_ffio_data[ff_type] = dict()
        self.stored_ffio_data[ff_type]['ff_type'] = ff_type
        self.stored_ffio_data[ff_type]['ff_number'] = ff_number
        self.stored_ffio_data[ff_type]['entry_data'] = entry_data
        self.stored_ffio_data[ff_type]['entry_values'] = entry_values
        self.stored_ffio_data[ff_type]['entry_dict'] = entry_dict

    def retrieve_ffio_data(self, ff_type):
        return [self.stored_ffio_data[ff_type]['ff_number'],
                self.stored_ffio_data[ff_type]['entry_data'],
                self.stored_ffio_data[ff_type]['entry_values'],
                self.stored_ffio_data[ff_type]['entry_dict']
                ]

    def parse_vdwtypes(self, type, current_molecule_type):
        ff_number, entry_data, entry_values, entry_dict = self.retrieve_ffio_data(type)
        # molecule name is at sites, but vdwtypes come
        # before sites. So we store info in vdwtypes and
        # edit it later at sites. Eventually, we should
        # probably move to a model where we store sections
        # we can't use yet, and then process them in the
        # order we want.
        logger.debug("Parsing [ vdwtypes ] ...")

        for j in range(ff_number):
            self.vdwtypes.append(entry_values[j].split()[3:]) #THIS IS ASSUMING ALL VDWTYPES ARE STORED AS LJ12_6_SIG_EPSILON
            self.vdwtypeskeys.append(entry_values[j].split()[1])

    def parse_sites(self, type, molname, i, start):
        ff_number, entry_data, entry_values, entry_dict = self.retrieve_ffio_data(type)

        #correlate with atomtypes and atoms in GROMACS
        logger.debug("Parsing [ sites ] ...")

        #set indices to avoid continually calling list functions.
        ivdwtype = entry_data.index('s_ffio_vdwtype')+1
        icharge = entry_data.index('r_ffio_charge')+1
        imass = entry_data.index('r_ffio_mass')+1
        stemp = None
        etemp = None
        if 'i_ffio_resnr' in entry_data:
            iresnum = entry_data.index('i_ffio_resnr')+1
            iresidue = entry_data.index('s_ffio_residue')+1

        cgnr = 0
        # create the atom type container for the datax
        current_molecule_type = MoleculeType(name=molname)
        current_molecule_type.nrexcl = 0   #PLACEHOLDER FOR NREXCL...WE NEED TO FIND OUT WHERE IT IS
                                           #MRS: basically, we have to figure out the furthest number of bonds out
                                           # to exclude OR explicitly set gromacs exclusions. Either should work.
                                           # for now, we'll go with the latter

        self.system.add_molecule_type(current_molecule_type)

        current_molecule = Molecule(name=molname)  # should this be the same molname several as lines up?
        for j in range(ff_number):
            values = split_with_quotes(entry_values[j])
            if values[1] == "atom":
                if ('i_ffio_resnr' in entry_data):
                    atom = Atom(int(values[0]), values[ivdwtype],
                                int(values[iresnum]),
                                values[iresidue])
                else:
                    # No residuenr, means we will have identical atoms sharing this.
                    atom = Atom(int(values[0]), values[ivdwtype])

                atom.atomtype = (0, values[ivdwtype])
                atom.charge = (0, float(values[icharge])*units.elementary_charge)
                atom.mass = (0, float(values[imass]) * units.amu)
                stemp = float(self.vdwtypes[self.vdwtypeskeys.index(values[ivdwtype])][0]) * units.angstroms #was in angstroms
                etemp = float(self.vdwtypes[self.vdwtypeskeys.index(values[ivdwtype])][1]) * units.kilocalorie_per_mole #was in kilocal per mol
                atom.sigma = (0, stemp)
                atom.epsilon = (0, etemp)
                atom.cgnr = cgnr
                cgnr+=1

                newAtomType = None
                current_molecule.add_atom(atom)
                if not self.system._atomtypes.get(AbstractAtomType(atom.atomtype.get(0))): #if atomtype not in self.system, add it
                    if self.system.combination_rule == 'Multiply-C6C12':
                        sigma = (etemp/stemp)**(1/6)
                        epsilon = (stemp)/(4*sigma**6)
                        newAtomType = AtomCType(values[ivdwtypes],             #atomtype/name
                                      values[ivdwtype],                             #bondtype
                                      -1,                               #atomic_number
                                      float(values[imass]) * units.amu,      #mass
                                      float(values[icharge]) * units.elementary_charge,  #charge--NEED TO CONVERT TO ACTUAL UNIT
                                      'A',                             #pcharge...saw this in top--NEED TO CONVERT TO ACTUAL UNITS
                                      sigma * units.kilocalorie_per_mole * angstroms**(6),
                                      epsilon * units.kilocalorie_per_mole * unit.angstro,s**(12))
                    elif (self.system.combination_rule == 'Lorentz-Berthelot') or (self.system.combination_rule == 'Multiply-Sigeps'):
                        newAtomType = AtomSigepsType(values[ivdwtype], #atomtype/name
                                      values[ivdwtype],                 #bondtype
                                      -1,                   #atomic_number
                                      float(values[imass]) * units.amu,  #mass--NEED TO CONVERT TO ACTUAL UNITS
                                      float(values[icharge]) * units.elementary_charge,  #charge--NEED TO CONVERT TO ACTUAL UNIT
                                      'A',                  #pcharge...saw this in top--NEED TO CONVERT TO ACTUAL UNITS
                                      stemp,
                                      etemp)
                    self.system.add_atomtype(newAtomType)

        if len(self.atom_blockpos) > 1:  #LOADING M_ATOMS
            if self.atom_blockpos[0] < start:
                # generate the new molecules for this block; the number of molecules depends on
                # The number of molecules depends on the number of entries in ffio_sites (ff_number)
                new_molecules = self.loadMAtoms(self.lines, self.atom_blockpos[0], i, current_molecule, ff_number)
                self.atom_blockpos.pop(0)

        index = 0
        for molecule in new_molecules:
            self.system.add_molecule(molecule)
            # now construct an atomlist with all the atoms

            for atom in molecule.atoms:
                # does this need to be a deep copy?
                # tmpatom = copy.deepcopy(atom)
                #    tmpatom.index = index
                self.atomlist.append(atom)
                index +=1

        return self.system._molecule_types[molname]

    def parse_bonds(self, type, current_molecule_type, i, start):

        ff_number, entry_data, ev, ed = self.retrieve_ffio_data(type)

        if len(self.bond_blockpos) > 1:  #LOADING M_BONDS
            if self.bond_blockpos[0] < start:
                for molecule in iter(current_molecule_type.molecules):
                    npermol = len(molecule.atoms)
                    break
                # of the parsers, this is the only one that uses 'lines'. Can we remove?
                current_molecule_type.bond_forces = self.loadMBonds(self.lines, self.bond_blockpos[0], i, npermol)
                self.bond_blockpos.pop(0)

        logger.debug("Parsing [ bonds ]...")

        for j in range(ff_number):
            values = split_with_quotes(ev[j])
            key = values[ed['s_ffio_funct']].upper()
            atomnames = ['i_ffio_ai','i_ffio_aj']
            atoms = [int(values[ed[a]]) for a in atomnames]
            bondingtypes = [self.atomlist[atom-1].name for atom in atoms]
            atoms.extend(bondingtypes)

            cnames = ['r_ffio_c1','r_ffio_c2']
            params = [float(values[ed[x]]) for x in cnames]
            new_bond = self.create_forcetype(self.desmond_bonds[key], atoms, params)
            kwds = self.get_parameter_kwds_from_force(new_bond)
            new_bond = self.canonical_bond(new_bond, kwds, direction = 'into', name = key)

            # removing the placeholder from matoms (should be a better way to do this?)
            if new_bond:
                old_bond = current_molecule_type.match_bonds(new_bond)
                if old_bond:
                    new_bond.order = old_bond.order
                    current_molecule_type.bond_forces.remove(old_bond)
                current_molecule_type.bond_forces.add(new_bond)

    def parse_pairs(self, type, current_molecule_type):

        ff_number, entry_data, ev, ed = self.retrieve_ffio_data(type)
        logger.debug("Parsing [ pairs ] ...")

        for j in range(ff_number):
            ljcorr = False
            coulcorr = False
            new_pair = None
            values = split_with_quotes(ev[j])
            atomnames = ['i_ffio_ai','i_ffio_aj']
            atoms = [int(values[ed[a]]) for a in atomnames]
            bondingtypes = [self.atomlist[atom-1].name for atom in atoms]
            params = atoms + bondingtypes
            key = values[ed['s_ffio_funct']].upper()
            if key == "LJ12_6_SIG_EPSILON":
                new_pair = self.create_forcetype(LjSigepsPair, params, 
                                                 [float(values[ed['r_ffio_c1']]),float(values[ed['r_ffio_c2']])])
            elif key == "LJ" or key == "COULOMB":
                # I think we just need LjSigepsPair, not LjPair?
                new_pair = self.create_forcetype(LjDefaultPair, params, [0, 0])
                if key == "LJ":
                    ljcorr = float(values[ed['r_ffio_c1']])
                    new_pair.scaleLJ = ljcorr
                elif key == "COULOMB":
                    coulcorr = float(values[ed['r_ffio_c1']])
                    new_pair.scaleQQ = coulcorr
            else:
                warn("ReadError: didn't recognize type {:s} in line {:s}".format(key, ev[j]))

            # now, we catch the matches and read them into a single potential
            pair_match = current_molecule_type.match_pairs(new_pair)
            if pair_match:  # we found a pair with the same atoms; let's insert or delete information as needed.
                remove_old = False
                remove_new = False
                if isinstance(new_pair, LjSigepsPair) and isinstance(pair_match, LjDefaultPair) and pair_match.scaleQQ:
                    #Need to add old scaleQQ to this new pair
                    new_pair.scaleQQ = pair_match.scaleQQ
                    remove_old = True
                elif isinstance(pair_match, LjSigepsPair) and isinstance(new_pair, LjDefaultPair) and new_pair.scaleQQ:
                    #Need to add the scaleQQ to the old pair
                    pair_match.scaleQQ = new_pair.scaleQQ
                    remove_new = True
                elif isinstance(new_pair,LjDefaultPair) and isinstance(pair_match,LjDefaultPair):
                    if pair_match.scaleQQ and not new_pair.scaleQQ:
                        new_pair.scaleQQ = pair_match.scaleQQ
                        remove_old = True
                    elif not pair_match.scaleQQ and new_pair.scaleQQ:
                        pair_match.scaleQQ = new_pair.scaleQQ
                        remove_new = True
                    if pair_match.scaleLJ and not new_pair.scaleLJ:
                        new_pair.scaleLJ = pair_match.scaleLJ
                        remove_new = True
                    elif not pair_match.scaleLJ and new_pair.scaleLJ:
                        pair_match.scaleLJ = new_pair.scaleLJ
                        remove_old = True

                if remove_old:
                    current_molecule_type.pair_forces.remove(pair_match)
                if remove_new:
                    new_pair = None

            if coulcorr:
                self.system.coulomb_correction = coulcorr  # need this for gromacs to have the global declared
                #If we have difference between global and local, catch in gromacs.

            if ljcorr:
                self.system.lj_correction = ljcorr  # need this for gromacs to have the global declared
                #If we have difference between global and local, catch in gromacs.

            if new_pair:
                current_molecule_type.pair_forces.add(new_pair)
                # IMPORTANT: we are going to assume that all pairs are both LJ and COUL.
                # if COUL is not included, then it is because the charges are zero, and they will give the
                # same energy.  This could eventually be improved by checking versus the sites.

    def parse_angles(self, type, current_molecule_type):

        ff_number, entry_data, ev, ed = self.retrieve_ffio_data(type)
        logger.debug("Parsing [ angles ] ...")

        for j in range(ff_number):
            values = split_with_quotes(ev[j])
            key = values[ed['s_ffio_funct']].upper()
            atomnames = ['i_ffio_ai','i_ffio_aj','i_ffio_ak']
            atoms = [int(values[ed[a]]) for a in atomnames]
            bondingtypes = [self.atomlist[atom-1].name for atom in atoms]
            atoms.extend(bondingtypes)
            kwds = [float(values[ed['r_ffio_c1']]), float(values[ed['r_ffio_c2']])]
            new_angle = self.create_forcetype(self.desmond_angles[key], atoms, kwds)
            kwds = self.get_parameter_kwds_from_force(new_angle)
            new_angle = self.canonical_angle(new_angle, kwds, direction = 'into', name = key,
                                             molecule_type = current_molecule_type)
            if new_angle:
                current_molecule_type.angle_forces.add(new_angle)

    def parse_dihedrals(self, type, current_molecule_type):

        ff_number, entry_data, ev, ed = self.retrieve_ffio_data(type)
        logger.debug("Parsing [ dihedrals ] ...")

        for j in range(ff_number):
            values = split_with_quotes(ev[j])
            new_dihedral = None
            dihedral_type = None
            key = values[ed['s_ffio_funct']].upper()
            atomnames = ['i_ffio_ai','i_ffio_aj','i_ffio_ak','i_ffio_al']
            atoms = [int(values[ed[a]]) for a in atomnames]
            bondingtypes = [self.atomlist[atom-1].name for atom in atoms]
            atoms.extend(bondingtypes)
            # not sure how to put the following lines in canonical, since it expects keywords,
            # not strings of variable length.  will have to fix later.
            cnames = []
            for e in entry_data:
                if 'r_ffio_c' in e:
                    cnames.append(e)   # append all of the constants, in order

            if key == "IMPROPER_HARM":
                kwds = [float(values[ed[cnames[0]]]), 2*float(values[ed[cnames[1]]])] # harmonic, multiple x2 for desmond convention.
            elif key == "PROPER_TRIG" or key == "IMPROPER_TRIG":
                kwds = [float(values[ed[x]]) for x in cnames]
            elif key == "OPLS_PROPER" or key == "OPLS_IMPROPER":
                opls_vals = [float(values[ed[x]]) * units.kilocalorie_per_mole for x in cnames[1:5]]
                opls_kwds = dict(zip([x[-2:] for x in cnames[1:5]],opls_vals))
                opls_kwds = convert_dihedral_from_fourier_to_trig(opls_kwds)
                kwds = np.zeros(8) # will fill this in later.
            new_dihedral = self.create_forcetype(self.desmond_dihedrals[key], atoms, kwds)
            # really should be some way to get rid of this code below
            if key == "OPLS_PROPER" or key == "OPLS_IMPROPER":
                for key in opls_kwds.keys():
                    setattr(new_dihedral,key,opls_kwds[key])
            # really should be some way to get rid of this code above
            kwds = self.get_parameter_kwds_from_force(new_dihedral)
            new_dihedral = self.canonical_dihedral(new_dihedral, kwds, direction = 'into', name = key,
                                                   molecule_type = current_molecule_type)

            if new_dihedral:
                current_molecule_type.dihedral_forces.add(new_dihedral)

    def parse_torsion_torsion(self, type, current_molecule_type):

        ff_number, entry_data, ev, ed = self.retrieve_ffio_data(type)
        logger.debug("Parsing [ torsion-torsion ] ...")

        for j in range(ff_number):
            values = split_with_quotes(ev[j])
            new_torsiontorsion = None
            key = values[ed['s_ffio_funct']].upper()
            if key == "CMAP":
                # we shouldn't need to try/accept because there are no units.
                new_torsiontorsion = TorsionTorsionCMAP(int(values[ed['i_ffio_ai']]),
                                                        int(values[ed['i_ffio_aj']]),
                                                        int(values[ed['i_ffio_ak']]),
                                                        int(values[ed['i_ffio_al']]),
                                                        int(values[ed['i_ffio_am']]),
                                                        int(values[ed['i_ffio_an']]),
                                                        int(values[ed['i_ffio_ao']]),
                                                        int(values[ed['i_ffio_ap']]),
                                                        'cmap',
                                                        int(values[ed['i_ffio_c1']]))
            else:
                warn("ReadError: found unsupported torsion-torsion type in: {:s}".format(str(line[i])))
            if new_torsiontorsion:
                current_molecule_type.torsiontorsion_forces.add(new_torsiontorsion)

    def parse_exclusions(self, type, current_molecule_type):

        ff_number, entry_data, ev, entry_dict = self.retrieve_ffio_data(type)
        logger.debug("Parsing [ exclusions ] ...")

        # currently, assumes no comments, could be dangerous?
        for j in range(ff_number):
            temp = split_with_quotes(ev[j])
            temp.remove(temp[0])
            current_molecule_type.exclusions.add(tuple([int(x) for x in temp]))

    def parse_restraints(self, type, current_molecule_type):

        ff_number, entry_data, ev, ed = self.retrieve_ffio_data(type)
        logger.debug("Warning: Parsing [ restraints] not yet implemented")

    def parse_constraints(self, type, current_molecule_type):

        ff_number, entry_data, ev, ed = self.retrieve_ffio_data(type)
        logger.debug("Parsing [ constraints ] ...")

        ctype = 1
        funct_pos = 0
        atompos = [] #position of atoms in constraints; spread all over the place
        lenpos = [] #position of atom length; spread all over the place
        tempatom = []
        templength = []
        templen = 0
        for j in range(len(entry_data)):
            if entry_data[j] == 's_ffio_funct':
                funct_pos = ctype
            elif 'i_ffio' in entry_data[j]:
                atompos.append(ctype)
            elif 'r_ffio' in entry_data[j]:
                lenpos.append(ctype)
            ctype+=1

        for j in range(ff_number):
            # water constraints actually get written to rigidwater (i.e. settles) constraints.

            if 'HOH' in ev[j] or 'AH' in ev[j]:
                values = split_with_quotes(ev[j])
                tempatom = []
                templength = []
                for a in atompos:
                    if not '<>' in values[a]:
                        tempatom.append(int(values[a]))
                    else:
                        tempatom.append(None)
                for l in lenpos:
                    if not '<>' in values[l]:
                        if 'AH' in ev[j]:
                            templength.append(float(values[l])*units.angstroms) # Check units?
                    else:
                        templength.append(None*units.angstroms)

                constr_type = values[funct_pos]
                if 'HOH' in constr_type:
                    dOH = float(values[lenpos[1]])
                    if dOH != float(values[lenpos[2]]):
                        logger.debug("Warning: second length in a rigid water specification {:s} is not the same as the first {:s}".format(values[lenpos[1]],values[lenpos[2]]))
                    angle = float(values[lenpos[0]])/(180/math.pi)
                    dHH = 2*dOH*math.sin(angle/2)
                    params = [atompos[0], atompos[1], atompos[2], dOH*units.angstroms, dHH*units.angstroms]
                    new_rigidwater = RigidWater(*params)
                    if new_rigidwater:
                        current_molecule_type.rigidwaters.add(new_rigidwater)

                elif 'AH' in constr_type:
                    templen = int(list(constr_type)[-1])
                    params = [tempatom[0], tempatom[1], templength[0], constr_type]
                    for t in range(2,templen+1):
                        params.extend([tempatom[t],templength[t-1]])
                    new_constraint = Constraint(*params)
                    if new_constraint:
                        current_molecule_type.constraints.add(new_constraint)
            else:
                warn("ReadError: found unsupported constraint type {:s}".format(ev[j]))

    def load_ffio_block(self, molname, start, end):

#        Loading in ffio blocks from Desmond format
#        Args:
#            molname: name of current molecule
#            start: beginning of where ffio_ff starts for each molecule
#            end: ending of where ffio_ff ends for each molecule

        i = start
        j = start

        self.stored_ffio_types = []  # a list of stored ffio_type to keep track
                              # of the ordering later
        self.stored_ffio_data = {}  # dictionary of stored ffio_entries

        values = []
        constraints = []
        temp = []

        current_molecule_type = None

        #There are several sections which require sites information to
        #process.  We keep a flag for this so that we are aware when
        #we have seen the sites
        bPreambleRead = False

        namecol = 0
        combrcol = 0
        vdwtypercol = 0

        #DEFAULT VALUES WHEN CONVERTING TO GROMACS
        self.system.nonbonded_function = 1
        self.system.genpairs = 'yes'

        logger.debug('Parsing [ molecule {:s} ]'.format(molname))
        logger.debug('Parsing [ ffio ]')

        while i < end:
            if not bPreambleRead:
                # read the first section for the forces field info
                while not (':::' in self.lines[i]):
                    if 's_ffio_name' in self.lines[i]:
                        namecol = i-start-1
                    elif 's_ffio_comb_rule' in self.lines[i]:
                        combrcol = i-start-1
                    elif 's_ffio_vdw_func' in self.lines[i]:
                        vdwtypercol = i-start-1
                    i+=1
                i+=1 # skip the ':::'
                # figure out combination rule
                combrule = self.lines[i+combrcol].upper()
                if "ARITHMETIC/GEOMETRIC" in combrule:
                    self.system.combination_rule = 'Lorentz-Berthelot'
                elif "GEOMETRIC" in combrule:
                    self.system.combination_rule = 'Multiply-Sigeps'
                elif "LJ12_6_C6C12" in combrule:
                    self.system.combination_rule = 'Multiply-C6C12'
                if (vdwtypercol > 0):
                    vdwrule = self.lines[i+vdwtypercol]
                # MISSING: need to identify vdw rule here -- currently assuming LJ12_6_sig_epsilon!

                # skip to the next ffio entry
                while not ('ffio' in self.lines[i]):
                    i+=1
                bPreambleRead = True

            ff_type, ff_number, entry_data, entry_values, entry_dict, i = self.parse_ffio_block(i, end)

            self.stored_ffio_types.append(ff_type)
            self.store_ffio_data(ff_type,ff_number, entry_data, entry_values, entry_dict)

        # Reorder so 'vdwtypes' is first, then 'sites'.  Could eventually get some simplification
        # by putting sites first, but too much rewriting for now.

        self.stored_ffio_types.insert(0, self.stored_ffio_types.pop(self.stored_ffio_types.index('ffio_sites')))
        self.stored_ffio_types.insert(0, self.stored_ffio_types.pop(self.stored_ffio_types.index('ffio_vdwtypes')))

        # now process all the data
        for type in self.stored_ffio_types:
            if type in self.sysDirective:
                params = [type]
                if type == 'ffio_sites':
                    params += [molname, i, start]
                else:
                    params += [current_molecule_type]
                    if type == 'ffio_bonds':
                        params += [i, start]
                if type == 'ffio_sites':
                    current_molecule_type = self.sysDirective[type](*params)
                else:
                    self.sysDirective[type](*params)
            elif type == 'Done with ffio':
                continue
            else:
                while '}' not in self.lines[i]:
                    i+=1   # not the most robust if there is nesting in a particular pattern

    def loadMBonds(self, lines, start, end, npermol): #adds new bonds for each molecule in System

#        Loading in m_bonds in Desmond format
#        Args:
#            lines: list of all data in CMS format
#           start: beginning of where m_bonds starts for each molecule
#           end: ending of where m_bondsends for each molecule

        logger.debug("Parsing [ m_bonds ] ...")

        bg = False
        newbond_force = None
        values = []
        i = start
        bonds = set()
        # right now, hard coded with header:
        ## First column is bond index #
        # i_m_from
        # i_m_to
        # i_m_order
        # which might need to be read at some point.

        while i < end:
            if ':::' in lines[i]:
                if bg:
                    break
                else:
                    bg = True
                    i+=1
            if bg:
                values = split_with_quotes(lines[i])
                atomi = int(values[1])
                atomj = int(values[2])
                bondingtypei = self.atomlist[atomi-1].name
                bondingtypej = self.atomlist[atomj-1].name
                params = [atomi, atomj, bondingtypei, bondingtypej]
                if atomi > npermol:  # we've collected the number of atoms per molecule.  Exit.
                    break
                order = int(values[3])
                kwd = [0, 0]
                optkwd = {'order': order, 'c': False}
                new_bond = self.create_forcetype(HarmonicBond, params, kwd, optkwd)
                bonds.add(new_bond)
            i+=1

        return bonds

    def loadMAtoms(self, lines, start, end, currentMolecule, slength): #adds positions and such to atoms in each molecule in System

#        Loading in m_atoms from Desmond format
#        Args:
#            lines: list of all data in CMS format
#           start: beginning of where m_atoms starts for each molecule
#           end: ending of where m_atoms ends for each molecule
#           currentMolecule
#           slength: number of unique atoms in m_atoms, used to calculate repetitions

        logger.debug("Parsing [ m_atom ] ...")

        i = start
        bg = False
        pdbaname = ""
        aname = ""

        mult = int(re.split('\W',lines[start].split()[0])[1])/slength

        cols = dict()
        while i < end:
            if ':::' in lines[i]:
                i+=1
                break
            else:
                if 'First column' in lines[i]:
                    start += 1
                for c in self.atom_col_vars:
                    if c in lines[i]:
                        logger.debug("   Parsing [ {:s} ] ...".format(c))
                        cols[c] = i - start
                        break
            i+=1

        logger.debug("   Parsing atoms...")

        atom = None
        newMoleculeAtoms = []
        molecules = []
        j = 0
        while j < mult:
            newMolecule = copy.deepcopy(currentMolecule)
            for atom in newMolecule.atoms:
                if ':::' in lines[i]:
                    break
                else:
                    aline = split_with_quotes(lines[i])
                    atom.residue_index = int(aline[cols['i_m_residue_number']])
                    atom.residue_name = aline[cols['s_m_pdb_residue_name']].strip()
                    try:
                        atom.atomic_number = int(aline[cols['i_m_atomic_number']])
                    except Exception as e:
                        logger.exception(e) # EDZ: just pass statement before, now exception is recorded, but supressed

                    atom.position = [float(aline[cols['r_m_x_coord']]) * units.angstroms,
                                     float(aline[cols['r_m_y_coord']]) * units.angstroms,
                                     float(aline[cols['r_m_z_coord']]) * units.angstroms]
                    atom.velocity = [0.0 * units.angstroms * units.picoseconds**(-1),
                                     0.0 * units.angstroms * units.picoseconds**(-1),
                                     0.0 * units.angstroms * units.picoseconds**(-1)]
                    if 'r_ffio_x_vel' in cols:
                        atom.velocity[0] = float(aline[cols['r_ffio_x_vel']]) * units.angstroms * units.picoseconds**(-1)
                    if 'r_ffio_y_vel' in cols:
                        atom.velocity[1] = float(aline[cols['r_ffio_y_vel']]) * units.angstroms * units.picoseconds**(-1)
                    if 'r_ffio_z_vel' in cols:
                        atom.velocity[2] = float(aline[cols['r_ffio_z_vel']]) * units.angstroms * units.picoseconds**(-1)

                    if 's_m_pdb_atom_name' in cols:
                        pdbaname = aline[cols['s_m_pdb_atom_name']].strip()
                    if 's_m_atom_name' in cols:
                        aname = aline[cols['s_m_atom_name']].strip()
                    if re.match('$^',pdbaname) and not re.match('$^',aname):
                        atom.name = aname
                    elif re.match('$^',aname) and not re.match('$^',pdbaname):
                        atom.name = pdbaname
                    elif re.search("\d+",pdbaname) and not re.search("\d+",aname):
                        if re.search("\D+",pdbaname) and re.search("\w+",pdbaname):
                            atom.name = pdbaname
                        else:
                            atom.name = aname
                    elif re.search("\d+",aname) and not re.search("\d+",pdbaname):
                        if re.search("\D+",aname) and re.search("\w+",aname):
                            atom.name = aname
                        else:
                            atom.name = pdbaname
                    elif re.match('$^',pdbaname) and re.match('$^',aname):
                        atom.name = "None"
                    else:
                        atom.name = aname  #doesn't matter which we choose, so we'll go with atom name instead of pdb
                    i+=1

            molecules.append(newMolecule)
            j+=1

        return molecules


    def load_box_vector(self, lines, start, end):

#       Loading Box Vector
#       Create a Box Vector to load into the System
#        Args:
#            lines: all the lines of the file stored in an array
#            start: starting position
#            end: ending position

        v = np.zeros([3, 3]) * units.angstroms
        for i, line in enumerate(lines[start:end]):
            if self.atom_box_vars[0] in line:
                startboxlabel = i
            if ':::' in line:
                endlabel = i + start
                break
        startbox = startboxlabel + endlabel

        for nvec, line in enumerate(lines[startbox:startbox + 9]):
            j = nvec // 3
            k = nvec % 3
            v[j, k] = float(line.strip()) * units.angstrom

        self.system.box_vector = v

    def read(self):

#        Load in data from file
#       Read data in Desmond format
#        Args:

        molnames = []
        with open(self.cms_file, 'r') as fl:
            self.lines = list(fl)
        i=0
        j=0

        self.atomtypes = dict()

        self.atomlist = []

        # figure out on which lines the different blocks begin and end.
        for line in self.lines:
            if 'f_m_ct' in line:
                if j > 0:
                    self.fmct_blockpos.append(i)
                j+=1
            if 'm_atom' in line and not (('i_m' in line) or ('s_m' in line)):
                if j > 1:
                    self.atom_blockpos.append(i)
                j+=1
            if 'm_bond' in line:
                if j > 2:
                    self.bond_blockpos.append(i)
                j+=1
            if 'ffio_ff' in line:
                if j > 2:
                    self.ffio_blockpos.append(i)
                j+=1
            i+=1
        i-=1

        self.fmct_blockpos.append(i)
        self.atom_blockpos.append(i)
        self.bond_blockpos.append(i)
        self.ffio_blockpos.append(i)

        self.sysDirective = {'ffio_vdwtypes': self.parse_vdwtypes,
                             'ffio_sites': self.parse_sites,
                             'ffio_bonds': self.parse_bonds,
                             'ffio_pairs': self.parse_pairs,
                             'ffio_angles': self.parse_angles,
                             'ffio_dihedrals': self.parse_dihedrals,
                             'ffio_torsion_torsion': self.parse_torsion_torsion,
                             'ffio_constraints': self.parse_constraints,
                             'ffio_exclusions': self.parse_exclusions,
                             'ffio_restraints': self.parse_restraints
                             }

        #LOADING Ffio blocks
        logger.debug("Reading ffio block...")
        #MRS: warning -- currently no check to avoid duplicated molecule names. Investigate.
        i = 0
        j = 0
        while i < (len(self.ffio_blockpos)-1):
            j = self.fmct_blockpos[i]
            while ':::' not in self.lines[j]:
                j+=1
            # make sure we have reasonable molecular names.
            molname = self.lines[j+1].strip()
            molname = molname.replace("\"","")  # get rid of quotation marks so we can find unique names
            if molname == "":
                molname = "Molecule_"+str(len(molnames)+1)
            molnames.append(molname)
            self.load_ffio_block(molname, self.ffio_blockpos[i], self.fmct_blockpos[i+1]-1)
            i+=1
        i = 0

        #LOAD RAW BOX VECTOR-Same throughout cms

        logger.debug("Reading Box Vector...")
        self.load_box_vector(self.lines, self.fmct_blockpos[0], self.atom_blockpos[0])

        return self.system

    def write_vdwtypes_and_sites(self, molecule):

        logger.debug("   -Writing vdwtypes...")

        i = 0
        sites = []
        vdwtypes = []
        sig = None
        ep = None
        stemp = None
        etemp = None
        combrule = self.system.combination_rule
        for atom in molecule.atoms:
            i+=1
            if atom.residue_index:
                sites.append(' {:3d} {:5s} {:9.8f} {:9.8f} {:2s} {:1d} {:4s}\n'.format(
                        i, 'atom',
                        atom._charge[0].value_in_unit(units.elementary_charge),
                        atom._mass[0].value_in_unit(units.atomic_mass_unit),
                        atom.atomtype[0], atom.residue_index, atom.residue_name))
            else:
                sites.append(' {:3d} {:5s} {:9.8f} {:9.8f} {:2s}\n'.format(
                        i, 'atom',
                        atom._charge[0].value_in_unit(units.elementary_charge),
                        atom._mass[0].value_in_unit(units.atomic_mass_unit),
                        atom.atomtype[0]))

            sig = float(atom.sigma[0].value_in_unit(units.angstroms))
            ep = float(atom.epsilon[0].value_in_unit(units.kilocalorie_per_mole))
            if combrule == 'Multiply-C6C12':   #MRS: seems like this should be automated more?
                stemp = ep * (4 * (sig**6))
                etemp = stemp * (sig**6)
            elif combrule in ['Lorentz-Berthelot','Multiply-Sigeps']:
                stemp = sig
                etemp = ep
            vdwstring = ' {:2s} {:18s} {:8.8f} {:8.8f}\n'.format(atom.atomtype[0], 
                                                                 "LJ12_6_sig_epsilon", 
                                                                 float(stemp), float(etemp))
            if vdwstring not in vdwtypes:
                vdwtypes.append(vdwstring)

        lines = []
        lines.append("    ffio_vdwtypes[{:d}] {{\n".format(len(vdwtypes)))
        lines.append("      s_ffio_name\n")
        lines.append("      s_ffio_funct\n")
        lines.append("      r_ffio_c1\n")
        lines.append("      r_ffio_c2\n")
        lines.append("      :::\n")
        i = 0
        for v in vdwtypes:
            i+=1
            lines.append('      {:d}{:2s}'.format(i,v))
        lines.append("      :::\n")
        lines.append("    }\n")

        logger.debug("   -Writing sites...")

        lines.append("    ffio_sites[{:d}] {{\n".format(len(sites)))
        lines.append("      s_ffio_type\n")
        lines.append("      r_ffio_charge\n")
        lines.append("      r_ffio_mass\n")
        lines.append("      s_ffio_vdwtype\n")
        if len(sites[0].split()) > 5:   # fix this to explicitly ask if resnr is in here rather than length
            lines.append("      i_ffio_resnr\n")
            lines.append("      s_ffio_residue\n")
        lines.append("      :::\n")
        for s in sites:
            lines.append('   {:s}'.format(s))

        lines.append("      :::\n")
        lines.append("    }\n")

        return lines

    def write_bonds(self, moleculetype):

        logger.debug("   -Writing bonds...")

        dlines = list()
        hlines = list()
        hlines.append('ffio_bonds_placeholder\n')
        hlines.append("      i_ffio_ai\n")
        hlines.append("      i_ffio_aj\n")
        hlines.append("      s_ffio_funct\n")
        hlines.append("      r_ffio_c1\n")
        hlines.append("      r_ffio_c2\n")
        hlines.append("      :::\n")

        i = 0
        bondlist = sorted(list(moleculetype.bond_forces), key=lambda x: (x.atom1, x.atom2))
        for bond in bondlist:
            atoms = [bond.atom1,bond.atom2]
            kwds = self.get_parameter_kwds_from_force(bond)
            names, paramlists = self.canonical_bond(bond, kwds, direction = 'from')
            # could in general return multiple types and paramlists
            for nbond, name in enumerate(names):
                i += 1
                converted_bond = self.desmond_bonds[name](*atoms, **paramlists[nbond])
                line = '      {:d} {:d} {:d} {:s}'.format(i, atoms[0], atoms[1], name)
                bond_params = self.get_parameter_list_from_force(converted_bond)
                param_units = self.unitvars[converted_bond.__class__.__name__]
                for param, param_unit in zip(bond_params, param_units):
                    line += " {:15.8f}".format(param.value_in_unit(param_unit))
                line += '\n'
                dlines.append(line)
        header = "    ffio_bonds[{:d}] {{\n".format(i)
        hlines = end_header_section(i==0,header,hlines)

        dlines.append("      :::\n")
        dlines.append("    }\n")
        hlines.extend(dlines)
        return hlines

    def write_angles(self, moleculetype):

        logger.debug("   -Writing angles...")

        dlines = list()
        hlines = list()
        hlines.append("    ffio_angles_placeholder\n")
        hlines.append("      i_ffio_ai\n")
        hlines.append("      i_ffio_aj\n")
        hlines.append("      i_ffio_ak\n")
        hlines.append("      s_ffio_funct\n")
        hlines.append("      r_ffio_c1\n")
        hlines.append("      r_ffio_c2\n")
        hlines.append("      :::\n")

        i = 0
        anglelist = sorted(list(moleculetype.angle_forces), key=lambda x: (x.atom1,x.atom2,x.atom3))
        for angle in anglelist:
            atoms = [angle.atom1,angle.atom2,angle.atom3]
            kwds = self.get_parameter_kwds_from_force(angle)
            names, paramlists = self.canonical_angle(angle, kwds, direction = 'from')
            # could return multiple names and kwd lists
            for nangle, name in enumerate(names):
                i+=1
                converted_angle = self.desmond_angles[name](*atoms, **paramlists[nangle])
                line = '      {:d} {:d} {:d} {:d} {:s}'.format(i, atoms[0], atoms[1], atoms[2], name)
                angle_params = self.get_parameter_list_from_force(converted_angle)
                param_units = self.unitvars[converted_angle.__class__.__name__]
                for param, param_unit in zip(angle_params, param_units):
                    line += " {:15.8f}".format(param.value_in_unit(param_unit))
                line += '\n'
                dlines.append(line)

        header = "    ffio_angles[{:d}] {{\n".format(i)
        hlines = end_header_section(i==0,header,hlines)

        dlines.append("      :::\n")
        dlines.append("    }\n")
        hlines.extend(dlines)
        return hlines

    def write_dihedrals(self, moleculetype):

        logger.debug("   -Writing dihedrals...")

        dlines = list()
        hlines = list()
        hlines.append("    ffio_dihedrals_placeholder\n")
        hlines.append("      i_ffio_ai\n")
        hlines.append("      i_ffio_aj\n")
        hlines.append("      i_ffio_ak\n")
        hlines.append("      i_ffio_al\n")
        hlines.append("      s_ffio_funct\n")
        # we assume the maximum number of dihedral terms

        hmax = 8
        # assume the maximum number of dihedral terms (8) to simplify things for now
        for ih in range(hmax):
            hlines.append("      r_ffio_c{:d}\n".format(ih))
        hlines.append("      :::\n")

        i = 0

        # sorting dihedrals by first index
        dihedrallist = sorted(list(moleculetype.dihedral_forces), key=lambda x: (x.atom1, x.atom2, x.atom3, x.atom4))
        for dihedral in dihedrallist:
            atoms = [dihedral.atom1,dihedral.atom2,dihedral.atom3,dihedral.atom4]
            kwds = self.get_parameter_kwds_from_force(dihedral)
            names, paramlists = self.canonical_dihedral(dihedral, kwds, direction = 'from')
            for ndihedrals, name in enumerate(names):
                i+=1
                line = '      {:d} {:d} {:d} {:d} {:d} {:s}'.format(i, atoms[0], atoms[1], atoms[2], atoms[3], name)
                converted_dihedral= self.desmond_dihedrals[name](*atoms,**paramlists[ndihedrals])
                dihedral_params = self.get_parameter_list_from_force(converted_dihedral)
                param_units = self.unitvars[converted_dihedral.__class__.__name__]
                for param, param_unit in zip(dihedral_params, param_units):
                    line += " {:15.8}".format(param.value_in_unit(param_unit))
                for j in range(8-len(dihedral_params)):
                    line += " {:6.3f}".format(0.0)
                line += '\n'
                dlines.append(line)

        header = "    ffio_dihedrals[{:d}] {{\n".format(i)
        hlines = end_header_section(i==0,header,hlines)

        dlines.append("      :::\n")
        dlines.append("    }\n")
        hlines.extend(dlines)
        return hlines

    def write_torsion_torsion(self, moleculetype):

        # currently just cmap
        logger.debug("   -Writing torsion-torsions...")

        hlines = list()
        dlines = list()
        hlines.append("    ffio_torsion_torsion_placeholder\n")
        hlines.append("      i_ffio_ai\n")
        hlines.append("      i_ffio_aj\n")
        hlines.append("      i_ffio_ak\n")
        hlines.append("      i_ffio_al\n")
        hlines.append("      i_ffio_am\n")
        hlines.append("      i_ffio_an\n")
        hlines.append("      i_ffio_ao\n")
        hlines.append("      i_ffio_ap\n")
        hlines.append("      s_ffio_func\n")
        hlines.append("      i_ffio_c1\n")
        hlines.append("      :::\n")
        i = 0

        for torsiontorsion in moleculetype.torsiontorsion_forces:
            i+=1
            # only type of torsion/torsion is CMAP currently
            dlines.append('      {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:s} {:d}\n'.format(
                    i,
                    int(torsiontorsion.atom1), int(torsiontorsion.atom2),
                    int(torsiontorsion.atom3), int(torsiontorsion.atom4),
                    int(torsiontorsion.atom5), int(torsiontorsion.atom6),
                    int(torsiontorsion.atom7), int(torsiontorsion.atom8),
                    'cmap', torsiontorsion.chart))
        header = "    ffio_torsion_torsion[{:d}] {{\n".format(i)

        hlines = end_header_section(i==0,header,hlines)
        dlines.append("      :::\n")
        dlines.append("    }\n")
        hlines.extend(dlines)
        dlines = list()

        # write out the cmap terms: for now, write out all the
        # charts.  Later, we can scan through and only print out the ones we use
        # and only include the relevant charts

        if (i > 0):  # only include cmap_charts if we need to
            cmap_charts = cmap_parameters.get_cmap_charts()
            for chart in cmap_charts:
                chartlines = chart.split('\n')
                for line in chartlines:
                    dlines.append(line + '\n')

        hlines.extend(dlines)
        return hlines

    def write_exclusions(self, moleculetype):

        logger.debug("   -Writing exclusions...")

        hlines = list()
        dlines = list()
        hlines.append("    ffio_exclusions_placeholder\n")
        hlines.append("      i_ffio_ai\n")
        hlines.append("      i_ffio_aj\n")
        hlines.append("      :::\n")

        i = 0
        if moleculetype.nrexcl == 0:

            # Should probably be determined entirely by the bonds,
            # since settles now adds bonds.  For now, leave this in
            # for Desmond to Desmond conversion, where nrexcl is not
            # determined.  Probably should switch eventually.

            exclusionlist = sorted(list(moleculetype.exclusions), key=lambda x: (x[0], x[1]))
            for exclusion in moleculetype.exclusions:
                i+=1
                dlines.append('      {:d} {:d} {:d}\n'.format(i, int(exclusion[0]), int(exclusion[1])))

        else:
            if moleculetype.nrexcl > 4:
                warn("Can't handle more than excluding 1-4 interactions right now!")

            fullbondlist = []
            fullbondlist = sorted(list(moleculetype.bond_forces), key=lambda x: (x.atom1, x.atom2))
            # exclude HarmonicPotential types, which do not have exclusions.
            bondlist = [bond for bond in fullbondlist if (not isinstance(bond, HarmonicPotentialBond))]

            # first, figure out the first appearance of each atom in the bondlist
            currentatom = 0
            atompos = []
            bondindex = 0
            for molecule in moleculetype.molecules:
                nsize = len(molecule.atoms)+1
                break  # only need the first
            atombonds = np.zeros([nsize,8],int)  # assume max of 8 for now
            natombonds = np.zeros(nsize,int)
            for bond in bondlist:
                atombonds[bond.atom1,natombonds[bond.atom1]] = bond.atom2
                natombonds[bond.atom1] += 1
                atombonds[bond.atom2,natombonds[bond.atom2]] = bond.atom1
                natombonds[bond.atom2] += 1

            for atom in range(1,nsize):
                atomexclude = set()  # will be a unique set

                # need to make this recursive! And there must be a better algorithm
                for j1 in range(natombonds[atom]):
                    toatom1 = atombonds[atom,j1];
                    atomexclude.add(toatom1)
                    if moleculetype.nrexcl > 1:
                        for j2 in range(natombonds[toatom1]):
                            toatom2 = atombonds[toatom1,j2]
                            atomexclude.add(toatom2)
                            if moleculetype.nrexcl > 2:
                                for j3 in range(natombonds[toatom2]):
                                    toatom3 = atombonds[toatom2,j3]
                                    atomexclude.add(toatom3)
                                    if moleculetype.nrexcl > 3:
                                        for j4 in range(natombonds[toatom3]):
                                            toatom4 = atombonds[toatom1,j4]
                                            atomexclude.add(toatom4)

                uniqueexclude = set(atomexclude)
                for a in atomexclude:
                    if (a > atom):
                        i+=1
                        dlines.append('      {:d} {:d} {:d}\n'.format(i, atom, a))


        header = "    ffio_exclusions[{:d}] {{\n".format(i)
        hlines = end_header_section(i==0,header,hlines)

        dlines.append("      :::\n")
        dlines.append("    }\n")
        hlines.extend(dlines)
        return hlines

    def write_pairs(self, moleculetype):

        logger.debug("   -Writing pairs...")

        dlines = list()
        hlines = list()
        hlines.append("ffio_pairs_placeholder\n")
        hlines.append("      i_ffio_ai\n")
        hlines.append("      i_ffio_aj\n")
        hlines.append("      s_ffio_funct\n")
        hlines.append("      r_ffio_c1\n")
        hlines.append("      r_ffio_c2\n")
        hlines.append("      :::\n")

        i = 0
        for pair in sorted(list(moleculetype.pair_forces), key=lambda x: (x.atom1, x.atom2)):
            atoms = ' {:d} {:d} '.format(pair.atom1, pair.atom2)
            # first, the COUL part.
            if pair.__class__ in (LjDefaultPair, LjqDefaultPair, LjSigepsPair, LjCPair):
                                 # the first two appear to be duplicates: consider merging.
                if pair.scaleQQ:
                    scaleQQ = pair.scaleQQ
                else:
                    scaleQQ = self.system.coulomb_correction
                i += 1
                dlines += '       {:d} {:s} Coulomb {:10.8f} <>\n'.format(i, atoms, scaleQQ)
            elif pair._class in (LjqSigepsPair, LjqCPair):
                warn("Desmond does not support pairtype {:s}!".format(pair.__class__.__name__ ))  # may not be true?
            else:
                warn("Unknown pair type {:s}!".format(pair.__class__.__name__))

            # now the LJ part.
            if pair.__class__ in (LjDefaultPair,LjqDefaultPair):
                if pair.scaleLJ:
                    scaleLJ = pair.scaleLJ
                else:
                    scaleLJ = self.system.lj_correction
                i += 1
                dlines += '       {:d} {:s} LJ {:10.8f} <>\n'.format(i, atoms, scaleLJ)
            elif pair.__class__ in (LjSigepsPair, LjqSigepsPair, LjCPair, LjqCPair):
                # Check logic here -- not clear that we can correctly determine which type it is.
                # Basically, I think it's whether scaleLJ is defined or not.
                if pair.__class__ in (LjCPair, LjqCPair):
                    epsilon = 0.25 * (pair.C6**2) / pair.C12    # (16*eps^2*sig^12 / 4 eps*sig^12) = 4 eps
                    sigma = (0.25 * pair.C6 / epsilon)**(1.0/6.0)  # (0.25 * 4 eps sig^6 / eps)^(1/6)
                elif pair.__class__ in (LjSigepsPair, LjqCPair):
                    epsilon = pair.epsilon
                    sigma = pair.sigma
                i += 1
                dlines += '       {:d} {:s} LJ12_6_sig_epsilon {:10.8f} {:10.8f}\n'.format(i, atoms,
                                                                               sigma.value_in_unit(units.angstroms),
                                                                               epsilon.value_in_unit(units.kilocalorie_per_mole))
            else:
                warn("Unknown pair type {:s}!".format(pair.__class__.__name__))

        header = "    ffio_pairs[{:d}] {{\n".format(i)
        hlines = end_header_section(i==0,header,hlines)

        dlines.append("      :::\n")
        dlines.append("    }\n")
        hlines.extend(dlines)
        return hlines

    def write_constraints(self, moleculetype):

        logger.debug("   -Writing constraints...")
        isHOH = False

        if len(moleculetype.rigidwaters) > 0:
            alen = 3
            clen = 3
        else:
            alen = 0
            clen = 0

        alen_max = alen
        clen_max = clen

        for constraint in moleculetype.constraints:
            if constraint.type[0:2] == 'AH':
                alen = constraint.n+1
                clen = alen-1
            if alen_max < alen:
                alen_max = alen
                clen_max = clen
        # we now know the maximum length of all constraint types


        # not sure we need to sort these, but makes it easier to debug
        constraintlist = sorted(list(moleculetype.constraints),key=lambda x: x.atom1)
        dlines = list()
        hlines = list()

        i = 0
        for constraint in constraintlist: #calculate the max number of atoms in constraint
            i+=1
            if constraint.type == 'HOH':
                cline = '      {:d} {:d} {:d} {:d} '.format(i,int(constraint.atom1),int(constraint.atom2),int(constraint.atom3))
                for j in range(alen_max-3):
                    cline += '0 '
                cline += constraint.type
                cline += ' {:10.8f}'.format(float(constraint.length1.value_in_unit(units.degrees)))
                cline += ' {:10.8f}'.format(float(constraint.length2.value_in_unit(units.angstroms)))
                cline += ' {:10.8f}'.format(float(constraint.length2.value_in_unit(units.angstroms)))
                for j in range(clen_max-3):
                    cline += ' <>'
            elif constraint.type[0:2] == 'AH':
                alen = constraint.n+1
                clen = alen-1
                catoms = [constraint.atom1]
                clengths = []
                for j in range(1,alen+1):
                    atomname = 'atom'+str(j+1)
                    lengthname = 'length'+str(j)
                    if hasattr(constraint,atomname):
                        catoms.append(getattr(constraint,atomname))
                        clengths.append(getattr(constraint,lengthname))
                cline = '      {:d} '.format(i)
                for j in range(alen):
                    cline += ' {:d} '.format(int(catoms[j]))
                for j in range(alen,alen_max):
                    cline += ' <> '
                cline += constraint.type
                for j in range(clen):
                    cline += ' {:10.8f}'.format(float(clengths[j].value_in_unit(units.angstroms)))
                for j in range(clen,clen_max):
                    cline += ' <>'
            cline += '\n'

            dlines.append(cline)

        # now need to add the constraints specified through settles.  Only one settles per molecule
        for rigidwater in moleculetype.rigidwaters:
            i += 1
            # Assumes the water arrangement O, H, H, which might not always be the case.  Consider adding detection.
            cline = '      {:d} {:d} {:d} {:d} '.format(i, rigidwater.atom1, rigidwater.atom2, rigidwater.atom3)
            for j in range(alen_max-3):
                cline += '0 '
            cline += ' HOH '
            dOH = rigidwater.dOH.value_in_unit(units.angstroms)
            dHH = rigidwater.dHH.value_in_unit(units.angstroms)
            angle = 2.0*math.asin(0.5*dHH/dOH)*(180/math.pi)    # could automate conversion. . .
            cline += ' {:.8f} {:.8f} {:.8f}'.format(angle,dOH,dOH)
            cline += '\n'
            for j in range(alen,alen_max):
                cline += ' 0.0'
            dlines.append(cline)

        hlines.append("    ffio_constraints[{:d}] {{\n".format(i))
        if (i==0):
            hlines.append("      :::\n")
        else:
            letters = ['i','j','k','l','m','n','o','p','q']
            for j in range(alen_max):
                hlines.append('      i_ffio_a{:s}\n'.format(letters[j]))
            hlines.append('      s_ffio_funct\n')
            for j in range(clen_max):
                hlines.append('      r_ffio_c{:d}\n'.format(j+1))
            hlines.append("      :::\n")

        dlines.append("      :::\n")
        dlines.append("    }\n")
        hlines.extend(dlines)
        return hlines

    def write(self):

#        Write this topology to file
#        Write out this topology in Desmond format
#        Args:
#            filename: the name of the file to write out to

        lines = list()
        pos = 0
        name = ''

        logger.warning("MacroModel atom type is not defined in other files, is set to 1 for all cases as it must be validly defined for desmond files to run.  However, it does not affect the energies.")

        # for all CMS files
        lines.append('{\n')
        lines.append('  s_m_m2io_version\n')
        lines.append('  :::\n')
        lines.append('  2.0.0\n')
        lines.append('}\n')

        #FIRST F_M_CT BLOCK

        logger.debug("Writing first f_m_ct...")
        lines.append('f_m_ct {\n')
        lines.append('  s_m_title\n')
        for c in self.atom_box_vars:
            lines.append('  {:s}\n'.format(c))
        lines.append('  s_ffio_ct_type\n')
        lines.append('  :::\n')

        #box vector
        bv = self.system.box_vector
        lines.append('  "Desmond file converted by InterMol"\n')
        for bi in range(3):
            for bj in range(3):
                lines.append('{:22.11f}\n'.format(float(bv[bi][bj].value_in_unit(units.angstroms))))
        lines.append('  full_system\n')       

        #M_ATOM
        apos = len(lines) #pos of where m_atom will be; will need to overwite later based on the number of atoms
        lines.append('m_atom\n')
        lines.append('    # First column is atom index #\n')
        for vars in self.atom_col_vars:
                if '_pdb_atom' not in vars:
                    lines.append('    {:s}\n'.format(vars))
        lines.append('    :::\n')

        i = 0
        nmol = 0
        totalatoms = []
        totalatoms.append(0)
        for moleculetype in self.system._molecule_types.values():
            for molecule in moleculetype.molecules:
                for atom in molecule.atoms:
                    i += 1
                    line = '    {:d}        {:d}'.format(i,1) #HAVE TO PUT THE 1 HERE OR ELSE DESMOND DIES, EVEN THOUGH IT DOESN'T USE IT
                    for j in range(3):
                        line += " {:10.8f}".format(float(atom._position[j].value_in_unit(units.angstroms)))
                    line +=   "     {:2d} {:4s}    {:2d}  {:2s}".format(
                        atom.residue_index,
                        '"{:s}"'.format(atom.residue_name),
                        atom.atomic_number,
                        '"{:s}"'.format(atom.name))
                    if np.any(atom._velocity):
                        for j in range(3):
                            line += " {:10.8f}".format(float(atom._velocity[j].value_in_unit(units.angstroms / units.picoseconds)))
                    else:
                         for j in range(3):
                            line += " {:10.8f}".format(0)
                    lines.append(line + '\n')
            totalatoms.append(i)

        lines[apos] = "  m_atom[{:d}] {{\n".format(i)
        lines.append('    :::\n')
        lines.append('  }\n')

        bpos = len(lines)
        i = 0

        #M_BOND
        hlines = list()
        dlines = list()
        hlines.append('  m_bond_placeholder\n')
        hlines.append('    i_m_from\n')
        hlines.append('    i_m_to\n')
        hlines.append('    i_m_order\n')
        hlines.append('    i_m_from_rep\n')
        hlines.append('    i_m_to_rep\n')
        hlines.append('    :::\n')

        i = 0
        nonecnt = 0
        for moleculetype in self.system._molecule_types.values():
            # sort the bondlist because Desmond requires the first time a bond is listed to have
            # the atoms in ascending order
            repeatmol = len(moleculetype.molecules)
            #MRS: need to be fixed; gromacs loads in one set of bonds per molecue; desmond loads in all

            # OrderedSet isn't indexable so get the first molecule by iterating.
            for molecule in moleculetype.molecules:
                atoms_per_molecule = len(molecule.atoms)
                # all should have the same, once we have info from one, break.
                break

            bondlist = sorted(list(moleculetype.bond_forces), key=lambda x: (x.atom1,x.atom2))
            for n in range(repeatmol):
                for bond in bondlist:
                    if bond and bond.order:
                        i += 1
                        dlines.append('    {:d} {:d} {:d} {:d} {:d} {:d}\n'.format(i,
                                        bond.atom1 + n*atoms_per_molecule + totalatoms[nmol],
                                        bond.atom2 + n*atoms_per_molecule + totalatoms[nmol],
                                        int(bond.order),
                                        1,
                                        1))
                    elif not bond:
                        nonecnt+=1
                if nonecnt > 0:
                    logger.debug('FOUND {:d} BONDS THAT DO NOT EXIST'.format(nonecnt))
            nmol +=1

        hlines[0] = '  m_bond[{:d}] {{\n'.format(i)
        if (i > 0):
            lines.extend(hlines)
            lines.extend(dlines)
            lines.append('    :::\n')
            lines.append('  }\n')
        lines.append('}\n')

        #WRITE OUT ALL FFIO AND F_M_CT BLOCKS

        for molecule_name, moleculetype in self.system.molecule_types.items():
            logger.debug('Writing molecule block {:s}...'.format(molecule_name))
            #BEGINNING BLOCK

            logger.debug("  Writing f_m_ct...")
            lines.append('f_m_ct {\n')
            lines.append('  s_m_title\n')
            for c in self.atom_box_vars:
                lines.append('  {:s}\n'.format(c))
            lines.append('  s_ffio_ct_type\n')
            lines.append('  :::\n')

            lines.append('  "' + molecule_name + '"\n')
            for bi in range(3):
                for bj in range(3):
                    lines.append('{:22.11f}\n'.format(float(bv[bi][bj].value_in_unit(units.angstroms))))
            lines.append('  solute\n')
            #M_ATOMS

            logger.debug("  Writing m_atoms...")
            apos = len(lines) #pos of where m_atom will be; will need to overwite later based on the number of atoms
            lines.append('m_atom\n')
            lines.append('    # First column is atom index #\n')
            for vars in self.atom_col_vars:
                if '_pdb_atom' not in vars:  # kludge, have better filter
                    lines.append('    {:s}\n'.format(vars))
            lines.append('    :::\n')

            i = 0
            for molecule in moleculetype.molecules:
                for atom in molecule.atoms:
                    i += 1
                    #NOT SURE WHAT TO PUT FOR MMOD TYPE; 1 is currently used.
                    #This can't be determined currently from the information provided,
                    # unless it is stored previous, nor is it used by desmond
                    line = '    {:d}        {:d}'.format(i,1)
                    for j in range(3):
                        line += " {:10.8f}".format(float(atom._position[j].value_in_unit(units.angstroms)))
                    line +=   "     {:2d} {:4s}    {:2d}  {:2s}".format(
                                atom.residue_index,
                                '"{:s}"'.format(atom.residue_name),
                                atom.atomic_number,
                                '"{:s}"'.format(atom.name))
                    if np.any(atom._velocity):
                        for j in range(3):
                            line += " {:10.8f}".format(float(atom._velocity[j].value_in_unit(units.angstroms / units.picoseconds)))
                    else:
                         for j in range(3):
                            line += " {:10.8f}".format(0)
                    lines.append(line + '\n')

            lines[apos] = '  m_atom[{:d}] {{\n'.format(i)
            lines.append('    :::\n')
            lines.append('  }\n')

            logger.debug("  Writing m_bonds...")

            hlines = list()
            dlines =  list()
            hlines.append('m_bond_placeholder\n')
            hlines.append('    i_m_from\n')
            hlines.append('    i_m_to\n')
            hlines.append('    i_m_order\n')
            hlines.append('    i_m_from_rep\n')
            hlines.append('    i_m_to_rep\n')
            hlines.append('    :::\n')

            i = 0
            nonecnt = 0

            repeatmol = len(moleculetype.molecules)

            for molecule in moleculetype.molecules:
                atoms_per_molecule = len(molecule.atoms)
                break

            bondlist = sorted(list(moleculetype.bond_forces), key=lambda x: x.atom1)
            for n in range(repeatmol):
                for bond in bondlist:
                    if bond and bond.order:
                        i += 1
                        dlines.append('    {:d} {:d} {:d} {:d} {:d} {:d}\n'.format(i,
                                        bond.atom1 + n*atoms_per_molecule,
                                        bond.atom2 + n*atoms_per_molecule,
                                        int(bond.order),
                                        1,
                                        1))
                    else:
                        nonecnt+=1
                if nonecnt > 0:
                    logger.debug('FOUND {:d} BONDS THAT DO NOT EXIST'.format(nonecnt))

            header = '  m_bond[{:d}] {{\n'.format(i)

            if (i>0):
                hlines = end_header_section(False,header,hlines)
                lines.extend(hlines)
                lines.extend(dlines)
                lines.append('    :::\n')
                lines.append('  }\n')

            # only need the first molecule
            molecule = next(iter(moleculetype.molecules))
            logger.debug("  Writing ffio...")
            lines.append('  ffio_ff {\n')
            lines.append('    s_ffio_name\n')
            lines.append('    s_ffio_comb_rule\n')
            lines.append('    i_ffio_version\n')
            lines.append('    :::\n')

            #Adding Molecule Name
            if "Viparr" in molecule_name:
                lines.append('    Generated by Viparr\n')
            else:
                lines.append('    "{:s}"\n'.format(molecule_name))

            #Adding Combination Rule
            if self.system.combination_rule == 'Multiply-C6C12':
                lines.append('    C6C12\n')   # this may not exist in DESMOND, or if so, need to be corrected
            elif self.system.combination_rule == 'Lorentz-Berthelot':
                lines.append('    ARITHMETIC/GEOMETRIC\n')
            elif self.system.combination_rule == 'Multiply-Sigeps':
                lines.append('    GEOMETRIC\n')

            #Adding Version
            lines.append('    1.0.0\n') #All files had this, check if version is 1.0.0

            lines += self.write_vdwtypes_and_sites(molecule)
            lines += self.write_bonds(moleculetype)
            lines += self.write_angles(moleculetype)
            lines += self.write_dihedrals(moleculetype)
            lines += self.write_torsion_torsion(moleculetype)
            lines += self.write_exclusions(moleculetype)
            lines += self.write_pairs(moleculetype)
            lines += self.write_constraints(moleculetype)

            #STILL NEED TO ADD RESTRAINTS

            lines.append("  }\n")
            lines.append("}\n")

        with open(self.cms_file, 'w') as fout:
            for line in lines:
                fout.write(line)
