from collections import OrderedDict
import logging
import warnings
import os
import math

import simtk.unit as units
from intermol.atom import Atom

from intermol.forces import *
import intermol.forces.forcefunctions as ff
from intermol.molecule import Molecule
from intermol.moleculetype import MoleculeType
from intermol.system import System
import cmap_parameters

import pdb #FOR DEBUGGING PURPOSES


logger = logging.getLogger('InterMolLog')

# driver helper functions
def load_desmond(infile, defines=None):
    """Load a DESMOND input file into a 'System'

    Args:
        infile:
        include_dir:
        defines:
    Returns:
        system:
    """
    parser = DesmondParser(infile, defines=defines)
    return parser.read()

def write_gromacs(outfile, system):
    """Unpacks a 'System' into a DESMOND input file

    Args:
        outfile:
        system:
    Returns:
        outfile:
    """
    parser = DesmondParser(outfile, system)
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
    line = list(line)
    in_quotes = False
    for i, char in enumerate(line):
        if char == '"':
            in_quotes == not in_quotes
        if char == ' ' and in_quotes:
            line[i] = '_'

    space_split = "".join(line).split()
    for i, sub in enumerate(space_split):
        sub = sub.replace('"', '')
        space_split[i] = sub.replace('_', ' ')
    return space_split

def to_lookup(forward_dict, key_rule=None, value_rule=None):
    temp_list = []

    for k, v in forward_dict.items():
        if key_rule:
            k = key_rule(k)
        if value_rule:
            v = value_rule(v)
        temp_list.append(v, k)
    
    return dict(temp_list)    

def to_type(forward_dict):
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
    lookup_desmond_combination_rules = to_lookup(desmond_combination_rules)

    desmond_pairs = {'LJ12_6_SIG_EPSILON': LjSigepsPair,
                     'LJ': LjDefaultPair,
                     'COULOMB': LjDefaultPair
                     }
    lookup_desmond_pairs = to_lookup(desmond_pairs) 
    desmond_pair_types = to_type(desmond_pairs)
     
    desmond_bonds = {'HARM_CONSTRAINED': HarmonicBond,
                     'HARM': HarmonicBond
                     }
    lookup_desmond_bonds = to_lookup(desmond_bonds)
    desmond_bond_types = to_type(desmond_bonds)


    def canonical_bond(self, name, kwds, bond, direction = 'into'):

        if direction == 'into':
            canonical_force_scale = self.canonical_force_scale_into
            phase = 'Read'
        else:
            canonical_force_scale = self.canonical_force_scale_from
            phase = 'Write'

        if bond in [HarmonicBond, HarmonicPotentialBond]:

            kwds['k'] = canonical_force_scale * kwds['k']

            if direction == 'into':
                if name == 'harm_constrainted':
                    kdws['c'] = True
            else:
                # harmonic potentials in Gromacs should be constrained
                optkwds = forcefunctions.optparamlookup(bond)
                if optkwds['c'] == True and not isinstance(bond, HarmonicPotentialBond):
                    name = 'harm_constrained'
                else:
                    name = 'harm'
                return name, kwds

        else:
            raise Exception("%sError: bondtype %s is not supported by Desmond" % (phase, bond.__class__.__name__))


    desmond_angles = {'HARM_CONSTRAINED': HarmonicAngle,
                      'HARM': HarmonicAngle,
                      'UB': UreyBradleyAngle
                      }
    lookup_desmond_angles = to_lookup(desmond_angles)
    desmond_angle_types = to_type(desmond_angles)

    def canonical_angle(self, params, angle, direction = 'into'):
        """
        Args:
            params:
            angle:
            direction:
        Returns:
        """

        if direction == 'into':
            canonical_force_scale = self.canonical_force_scale_into
            phase = 'Read'
        else:
            canonical_force_scale = self.canonical_force_scale_from
            phase = 'Write'

        if bond in [HarmonicAngle, UreyBradleyAngle]:

            kwds['k'] = canonical_force_scale * kwds['k']

            if direction == 'into':
                if name == 'harm_constrainted':
                    kdws['c'] = True
            else:
                # harmonic potentials in Gromacs should be constrained
                optkwds = forcefunctions.optparamlookup(bond)
                if optkwds['c'] == True and not isinstance(bond, HarmonicPotentialBond):
                    name = 'harm_constrained'
                else:
                    name = 'harm'
                return name, kwds


            if isinstance(angle,UreyBradleyAngle):
                dlines.append('      %d %d %d %d %s %10.8f %10.8f\n' % (i, angle.atom1, angle.atom2, angle.atom3, 'UB', float(angle.r.in_units_of(units.angstroms)._value), 0.5*float(angle.kUB.in_units_of(units.kilocalorie_per_mole*units.angstroms**(-2))._value)))
                i+=1
            if angle.c:
                dlines.append('      %d %d %d %d %s %10.8f %10.8f\n' % (i, angle.atom1, angle.atom2, angle.atom3, 'Harm_constrained', float(angle.theta.in_units_of(units.degrees)._value), 0.5*float(angle.k.in_units_of(units.kilocalorie_per_mole/units.radians**2)._value)))
            else:
                if isinstance(angle, HarmonicAngle) or isinstance(angle, UreyBradleyAngle):
                    dlines.append('      %d %d %d %d %s %10.8f %10.8f\n' % (i, angle.atom1, angle.atom2, angle.atom3, 'Harm', float(angle.theta.in_units_of(units.degrees)._value), 0.5*float(angle.k.in_units_of(units.kilocalorie_per_mole/units.radians**2)._value)))
                else:
                    raise Exception("WriteError: angletype %s is not supported by Desmond" % (angle.__class__.__name__))


        if direction == 'into':
            return angle, params
        else:
            a_type = self.lookup_desmond_angles[angle.__class__]
            if a_type:
                return a_type, params
            else:
                loger.warn("WriteError: found unsupported angle type {0}".format(
                    angle.__class__.__name__))

    desmond_dihedrals = {'IMPROPER_HARM': ImproperHarmonicDihedral,
                         'PROPER_TRIG': TrigDihedral,
                         'IMPROPER_TRIG': TrigDihedral,
                         'OPLS_PROPER': TrigDihedral,
                         'OPLS_IMPROPER': TrigDihedral
                         }
    lookup_desmond_dihedrals = {TrigDihedral: 'Trig',
                                ImproperHarmonicDihedral: 'IMPROPER_HARM'
                                }

    def canonical_dihedral(self, params, dihedral, direction='into'):

        d_type = self.lookup_desmond_angles[dihedral.__class__]
        if d_type:
            return d_type, params
        else:
            loger.warn("WriteError: found unsupported dihedral type {0}".format(
                dihedral.__class__.__name__))        


    def __init__(self, infile, defines=None):
        """
        Initializes a DesmondParse object which serves to read in a CMS file
        into the abstract representation.

        Args:
        defines: Sets of default defines to use while parsing.
        """
        self.includes = set()       # set storing includes
        self.defines = dict()        # list of defines
        self.comments = list()      # list of comments

        self.atomtypes = HashMap()
        self.bondtypes = HashMap()
        self.constrainttypes = HashMap()
        
        if defines:
            self.defines.union(defines)
        self.defines["FLEX_SPC"] = None
        self.defines["POSRE"] = None
        self.viparr = 1

        self.fblockpos = []
        self.a_blockpos = []
        self.b_blockpos = []
        self.ffio_blockpos = []

        self.stored_ffio_data = {}  # dictionary of stored ffio_entries

        self.paramlist = forcefunctions.build_paramlist('desmond')
        self.unitvars = forcefunctions.build_unitvars('desmond', self.paramlist)

        self.gromacs_bonds = {
            'harm': HarmonicBond,
            }

        # reverse the above dictionary
        #self.lookup_desmond_bonds = dict((v,k) for k,v in self.gromacs_bonds.items())
        self.lookup_desmond_bonds ={
            HarmonicBond : 'harm',
            HarmonicPotentialBond : 'harm'
            }

        self.canonical_force_scale_into = 2.0
        self.canonical_force_scale_from = 0.5

    def get_parameter_list_from_kwds(self, force, kwds):
        return forcefunctions.get_parameter_list_from_kwds(force, kwds, self.paramlist)

    def get_parameter_list_from_force(self, force):
        return forcefunctions.get_parameter_list_from_force(force, self.paramlist)

    def get_parameter_kwds_from_force(self, force):
        return forcefunctions.get_parameter_kwds_from_force(force, self.get_parameter_list_from_force, self.paramlist)

    def create_kwd_dict(self, forcetype_object, values, optvalues = None):
        kwd = forcefunctions.create_kwd_dict(self.unitvars, self.paramlist, forcetype_object, values, optvalues = optvalues)
        return kwd

    def create_forcetype(self, forcetype_object, paramlist, values, optvalues = None):
        return forcetype_object(*paramlist, **self.create_kwd_dict(forcetype_object, values, optvalues))

#LOAD FFIO BLOCKS IN FIRST (CONTAINS TOPOLOGY)

    def parse_ffio_block(self,lines,start,end):
            # read in a ffio_block that isn't ffio_ff and split it into the
            # commands and the values.
            # lots of room for additional error checking here, such as whether
            # each entry has the correct number of data values, whether they are the correct type, etc.

            # scroll to the next ffio entry
        while not 'ffio_' in lines[start]:
            # this is not an ffio block! or, we have reached the end of the file
            if ('ffio_' not in lines[start]): 
                start+=1
            if start >= end:
                return 'Done with ffio', 0, 0, 0, start

        lines[start].split()[1]
        components = re.split('\W', lines[start].split()[0]) # get rid of whitespace, split on nonword
        ff_type = components[0]
        ff_number = int(components[1])
        i = start+1
        entry_data = []
        while not ':::' in lines[i]:
            entry_data.append(lines[i].split()[0])
            i+=1
        i+=1 # skip the separator we just found
        entry_values = []
        while not ':::' in lines[i]:
            if lines[i].strip():  # skip the blank spaces.
                entry_values.append(lines[i])
            i+=1
        while '}' not in lines[i]:  # wait until we hit an end to the block
            i+=1
        i+=1 # step past the end of the block

        return ff_type,ff_number,entry_data,entry_values,i
    
    def store_ffio_data(self, ff_type, ff_number, entry_data, entry_values):

        self.stored_ffio_data[ff_type] = dict()
        self.stored_ffio_data[ff_type]['ff_type'] = ff_type
        self.stored_ffio_data[ff_type]['ff_number'] = ff_number
        self.stored_ffio_data[ff_type]['entry_data'] = entry_data
        self.stored_ffio_data[ff_type]['entry_values'] = entry_values

    def retrive_ffio_data(self, ff_type):
        return [self.stored_ffio_data[ff_type]['ff_type'],
                self.stored_ffio_data[ff_type]['ff_number'],
                self.stored_ffio_data[ff_type]['entry_data'],
                self.stored_ffio_data[ff_type]['entry_values']
                ]

    def parse_bonds(self, shared_args, sites_args):
        atomlist, entry_data, entry_values, ff_number, ff_type, vdwtypes, vdwtypeskeys = shared_args
        currentMolecule, i, lines, moleculeName, start, sysDirectiveAtm, verbose, currentMoleculeType = sites_args
        
        forces = []
        if len(self.b_blockpos) > 1:  #LOADING M_BONDS
            if self.b_blockpos[0] < start:
                npermol = len(currentMoleculeType.moleculeSet[0].atoms)
                forces = self.loadMBonds(atomlist, lines, self.b_blockpos[0], i, npermol, verbose)
                currentMoleculeType.bondForceSet = forces[0]
                self.b_blockpos.pop(0)
        logger.debug("Parsing [ bonds ]...")

        for j in range(ff_number):
            entries = entry_values[j].split()
            key = entries[3].upper()
            atoms = map(int, entries[1:3])
            bondingtypes = map(lambda atom: atomlist[atom-1].name, atoms)
            atoms.extend(bondingtypes)

            if key == "HARM_CONSTRAINED":
                constrained = 1
            elif key == "HARM":
                constrained = 0
            else:
                warnings.warn("ReadError: Found unsupported bond")
            optionalkwd = {'c': constrained}
            params = map(float, entries[4:6])

            #MRS: it would be best to create the keywords first, and then create the type
            #so that they keywords can be converted into canonical form. Otherwise we have
            #to operate on the structures
            
            newBondType = self.create_forcetype(HarmonicBondType, bondingtypes, params, optionalkwd)
            newBondForce = self.create_forcetype(HarmonicBond, atoms, params, optionalkwd)

            # KLUDGE:  This below should be done in a separate routine, and shouldn't have to be done twice
            newBondType.k *= self.canonical_force_scale_into
            newBondForce.k *= self.canonical_force_scale_into

            if newBondForce in currentMoleculeType.bondForceSet:
                oldBondForce = currentMoleculeType.bondForceSet.get(newBondForce)
                currentMoleculeType.bondForceSet.remove(newBondForce)
                newBondForce.order = oldBondForce.order
                currentMoleculeType.bondForceSet.add(newBondForce)

            if newBondType not in self.bondtypes:
	    	self.bondtypes.add(newBondType)

    def parse_vdwtypes(self, shared_args, sites_args):
        atomlist, entry_data, entry_values, ff_number, ff_type, vdwtypes, vdwtypeskeys = shared_args
        currentMolecule, i, lines, moleculeName, start, sysDirectiveAtm, verbose, currentMoleculeType = sites_args
        # molecule name is at sites, but vdwtypes come
        # before sites. So we store info in vdwtypes and
        # edit it later at sites. Eventually, we should
        # probably move to a model where we store sections
        # we can't use yet, and then process them in the
        # order we want.
        logger.debug("Parsing [ vdwtypes]...")
        
        for j in range(ff_number):
        	vdwtypes.append(entry_values[j].split()[3:]) #THIS IS ASSUMING ALL VDWTYPES ARE STORED AS LJ12_6_SIG_EPSILON
                vdwtypeskeys.append(entry_values[j].split()[1])

    def parse_sites(self, shared_args, sites_args):
        atomlist, entry_data, entry_values, ff_number, ff_type, vdwtypes, vdwtypeskeys = shared_args
        currentMolecule, i, lines, moleculeName, start, sysDirectiveAtm, verbose, currentMoleculeType = sites_args
        #correlate with atomtypes and atoms in GROMACS
        logger.debug("Parsing [ sites]...")
        
        #set indices to avoid continually calling list functions.    
        ivdwtype = entry_data.index('s_ffio_vdwtype')+1
        icharge = entry_data.index('r_ffio_charge')+1
        imass = entry_data.index('r_ffio_mass')+1
                
        if 'i_ffio_resnr' in entry_data:
            iresnum = entry_data.index('i_ffio_resnr')+1
            iresidue = entry_data.index('s_ffio_residue')+1
                  
        cgnr = 0
        for j in range(ff_number):
            split = entry_values[j].split()
            if split[1] == "atom":
                if ('i_ffio_resnr' in entry_data): 
                    atom = Atom(int(split[0]), split[ivdwtype],
                                int(split[iresnum]), 
                                split[iresidue])
                else:
                    # No residuenr, means we will have identical atoms sharing this.
                    atom = Atom(int(split[0]), split[ivdwtype])

                atom.atomtype = (0, split[ivdwtype])
                atom.charge = (0, float(split[icharge])*units.elementary_charge) 
                atom.mass = (0, float(split[imass]) * units.amu)
                stemp = float(vdwtypes[vdwtypeskeys.index(split[ivdwtype])][0]) * units.angstroms #was in angstroms
                etemp = float(vdwtypes[vdwtypeskeys.index(split[ivdwtype])][1]) * units.kilocalorie_per_mole #was in kilocal per mol
                atom.sigma = (0, stemp)
                atom.epsilon = (0, etemp)
                atom.cgnr = cgnr
                cgnr+=1

                newAtomType = None
                currentMolecule.add_atom(atom)
                if not System._sys._atomtypes.get(AbstractAtomType(atom.atomtype.get(0))): #if atomtype not in System, add it
                    if System._sys.combination_rule == 'Multiply-C6C12':
                        sigma = (etemp/stemp)**(1/6)
                        epsilon = (stemp)/(4*sigma**6)
                        newAtomType = AtomCType(split[ivdwtypes],             #atomtype/name
                                      split[ivdwtype],                             #bondtype
                                      -1,                               #atomic_number
                                      float(split[imass]) * units.amu,      #mass
                                      float(split[icharge]) * units.elementary_charge,  #charge--NEED TO CONVERT TO ACTUAL UNIT
                                      'A',                             #pcharge...saw this in top--NEED TO CONVERT TO ACTUAL UNITS
                                      sigma * units.kilocalorie_per_mole * angstroms**(6),
                                      epsilon * units.kilocalorie_per_mole * unit.angstro,s**(12))
                    elif (System._sys.combination_rule == 'Lorentz-Berthelot') or (System._sys.combination_rule == 'Multiply-Sigeps'):
                        newAtomType = AtomSigepsType(split[ivdwtype], #atomtype/name
                                      split[ivdwtype],                 #bondtype
                                      -1,                   #atomic_number
                                      float(split[imass]) * units.amu,  #mass--NEED TO CONVERT TO ACTUAL UNITS
                                      float(split[icharge]) * units.elementary_charge,  #charge--NEED TO CONVERT TO ACTUAL UNIT
                                      'A',                  #pcharge...saw this in top--NEED TO CONVERT TO ACTUAL UNITS
                                      stemp,
                                      etemp)
                    System._sys._atomtypes.add(newAtomType)

        if len(self.a_blockpos) > 1:  #LOADING M_ATOMS
            if self.a_blockpos[0] < start:
                # generate the new molecules for this block; the number of molecules depends on
                # The number of molecules depends on the number of entries in ffio_sites (ff_number)
                NewMolecules = self.loadMAtoms(lines, self.a_blockpos[0], i, currentMolecule, ff_number, sysDirectiveAtm, verbose)
                self.a_blockpos.pop(0)

        # now construct an atomlist with all the atoms
        index = 0
        for molecule in NewMolecules:
            System._sys.add_molecule(molecule)
            for atom in molecule.atoms:
                # does this need to be a deep copy?
                tmpatom = copy.deepcopy(atom)
                tmpatom.index = index
                atomlist.add(tmpatom)
                index +=1

        currentMoleculeType = System._sys._molecules[moleculeName]
        currentMoleculeType.nrexcl = 0 #PLACEHOLDER FOR NREXCL...WE NEED TO FIND OUT WHERE IT IS
                                       #MRS: basically, we have to figure out the furthest number of bonds out 
                                       # to exclude OR explicitly set gromacs exclusions. Either should work.
                                       # for now, we'll go with the latter
        return currentMoleculeType

    def parse_pairs(self, shared_args, sites_args):
        atomlist, entry_data, entry_values, ff_number, ff_type, vdwtypes, vdwtypeskeys = shared_args
        currentMolecule, i, lines, moleculeName, start, sysDirectiveAtm, verbose, currentMoleculeType = sites_args
        logger.debug("Parsing [ pairs]...")
        ljcorr = False
        coulcorr = False

        for j in range(ff_number):
            split = entry_values[j].split()
            newPairForce = None
            atoms = map(int,split[1:3])
            bondingtypes = map(lambda atom: atomlist[atom-1].name, atoms)
            params = atoms + bondingtypes
            key = split[3].upper()
            if key == "LJ12_6_SIG_EPSILON":
                newPairForce = self.create_forcetype(LjSigepsPair, params, map(float, split[4:6]))
            elif key == "LJ" or key == "COULOMB":
                # I think we just need LjSigepsPair, not LjPair
                newPairForce =  self.create_forcetype(LjDefaultPair, params, [0, 0])
                pairmatch = currentMoleculeType.pairForceSet.get(newPairForce)
                if key == "LJ":
                    ljcorr = float(split[4])
                    if pairmatch:
                        pairmatch.scaleLJ = ljcorr
                    else:
                        newPairForce.scaleLJ = ljcorr
                if key == "COULOMB":
                    coulcorr = float(split[4])  # we need this for the system, since it's not stored anywhere
                    if pairmatch:
                        pairmatch.scaleQQ = coulcorr
                    else:
                        newPairForce.scaleQQ = coulcorr
            else:
                warnings.warn("ReadError: didn't recognize type %s in line %s", split[3], entry_values[j])

            if coulcorr:
                System._sys.coulomb_correction = coulcorr  # need this for gromacs to have the global declared 
                #If we have difference between global and local, catch in gromacs.

            if ljcorr:
                System._sys.lj_correction = ljcorr  # need this for gromacs to have the global declared 
                #If we have difference between global and local, catch in gromacs.

            if newPairForce:
                currentMoleculeType.pairForceSet.add(newPairForce)
                # IMPORTANT: we are going to assume that all pairs are both LJ and COUL.
                # if COUL is not included, then it is because the charges are zero, and they will give the
                # same energy.  This could eventually be improved by checking versus the sites.

    def parse_angles(self, shared_args, sites_args):
        atomlist, entry_data, entry_values, ff_number, ff_type, vdwtypes, vdwtypeskeys = shared_args
        currentMolecule, i, lines, moleculeName, start, sysDirectiveAtm, verbose, currentMoleculeType = sites_args
        logger.debug("Parsing [ angles]...")
        for j in range(ff_number):
            split = entry_values[j].split()
            newAngleForce = None
            angleType = None
            optkwd = None
            atoms = map(int, split[1:4])
            bondingtypes = map(lambda atom: atomlist[atom-1].name, atoms)
            params = atoms + bondingtypes
            key = split[4].upper()
            # todo: integrate constraints and angle constraint description together better.
            if key == "HARM_CONSTRAINED":  # this needs to go first because HARM is a substring
                kwd = map(float, split[5:7])
                optkwd = {'c': True}
                angleType = HarmonicAngle
            elif key == "HARM":
                kwd = map(float, split[5:7])
                angleType = HarmonicAngle
            elif key == "UB":
                # Urey-Bradley is implemented in DESMOND differently, with the 
                # terms implemented in a new angle term independent of the harmonic term.
                # Instead, we will add everything together afterwards into a single term
                kwd = [0, 0] + map(float, split[5:7]) 
                angleType = UreyBradleyAngle
            else:
                warnings.warn("ReadError: found unsupported angle in: %s" %str(lines[i]))

            if angleType:
                newAngleForce = self.create_forcetype(angleType, params, kwd)
                newAngleForce.k *= self.canonical_force_scale_into

            if newAngleForce:

                # check to see if this angle already has terms associated with it.
                if (currentMoleculeType.angleForceSet.get(newAngleForce)):
                    # if it's already in, then reset the forces
                    oldAngleForce = currentMoleculeType.angleForceSet.get(newAngleForce)
                    # OK, we have two values.  One of them is UB.
                    if (oldAngleForce.k._value == 0):
                        # the old one is the UB term.  Copy the new angle parameters into it.
                        oldAngleForce.k = newAngleForce.k
                        oldAngleForce.theta = newAngleForce.theta
                        # overwrite the old one with the new one
                        newAngleForce = oldAngleForce
                    elif (newAngleForce.k._value == 0):
                        # the new one is the UB term.  Copy the old angle parameters into it
                        newAngleForce.k = oldAngleForce.k
                        newAngleForce.theta = oldAngleForce.theta
                    else:
                        warnings.warn("Duplicate angle type! Shouldn't reach this point of code!")

                # add it on
                currentMoleculeType.angleForceSet.add(newAngleForce)

    def parse_dihedrals(self, shared_args, sites_args):
        atomlist, entry_data, entry_values, ff_number, ff_type, vdwtypes, vdwtypeskeys = shared_args
        currentMolecule, i, lines, moleculeName, start, sysDirectiveAtm, verbose, currentMoleculeType = sites_args
        logger.debug("Parsing [ dihedrals]...")

        for j in range(ff_number):
            split = entry_values[j].split()
            newDihedralForce = None
            dihedralType = None
            atoms = map(int, split[1:5])
            bondingtypes = map(lambda atom: atomlist[atom-1].name, atoms)
            params = atoms + bondingtypes
            key = split[5].upper()
            #Improper Diehdral 2 ---NOT SURE ABOUT MULTIPLICITY
            # These two should be the same function.  Check differences (polymer or protein defn, etc).
            if key == "IMPROPER_HARM":
                kwd = [float(split[6]), 2*float(split[7])]
                dihedralType = ImproperHarmonicDihedral
            elif key == "PROPER_TRIG" or key == "IMPROPER_TRIG":
                if key == "IMPROPER_TRIG":
                    improper = True
                else:
                    improper = False

                # currently, I can't see a difference in Proper_trig and improper_trig angles.
                kwd = map(float, split[6:14]) + [improper]
                dihedralType = TrigDihedral
            elif key == "OPLS_PROPER" or key == "OPLS_IMPROPER":
                # Internal desmond formatting for OPLS_PROPER
                #
                # r_ffio_c0 = f0 (always 0)  split(6)
                # r_ffio_c1 = f1             split(7)
                # r_ffio_c2 = f2             split(8)
                # r_ffio_c3 = f3             split(9)
                # r_ffio_c4 = f4             split(10)
                # r_ffio_c5 = f5 (always 0)  split(11)
                # r_ffio_c6 = f6 (always 0)  split(12)
                opls_kwd = {key: value for key, value in zip("c1 c2 c3 c4".split(), map(lambda s: units.kilocalorie_per_mole*float(s), split[7:11]))}
                kwd = ConvertDihedralFromFourierToTrig(opls_kwd)
                kwd.update({'phi': 0 * units.degrees})
                dihedralType = TrigDihedral
            else:
                warnings.warn("ReadError: found unsupported dihedral in: %s" % str(line[i]))

            if dihedralType:
                if key == "OPLS_PROPER" or key == "OPLS_IMPROPER":
                    newDihedralForce = dihedralType(*params, **kwd)
                else:
                    newDihedralForce = self.create_forcetype(dihedralType, params, kwd)

            if newDihedralForce:
                try:
                    # we can have multiple parameters with DESMOND, and append if we do
                    dihedralmatch = currentMoleculeType.dihedralForceSet.get(newDihedralForce)
                    # this will fail if it's the wrong type of dihedral
                    if dihedralmatch:
                        dihedralmatch.sum_parameters(newDihedralForce) 
                except Exception as e:
                    logger.exception(e)
                currentMoleculeType.dihedralForceSet.add(newDihedralForce)
        #9 proper dihedrals, funct = 1
        #3 improper dihedrals, funct = 2
        #Ryckaert-Bellemans type dihedrals, funct = 3 and pairs are removed

    def parse_torsion_torsion(self, shared_args, sites_args):
        atomlist, entry_data, entry_values, ff_number, ff_type, vdwtypes, vdwtypeskeys = shared_args
        currentMolecule, i, lines, moleculeName, start, sysDirectiveAtm, verbose, currentMoleculeType = sites_args
        logger.debug("Parsing [ torsion-torsion]...")
        for j in range(ff_number):
            split = entry_values[j].split()
            newTorsionTorsionForce = None
            key = split[9].upper()
            if key == "CMAP":
                # we shouldn't need to try/accept because there are no units.
                newTorsionTorsionForce = TorsionTorsionCMAP(int(split[1]),
                                                            int(split[2]),
                                                            int(split[3]),
                                                            int(split[4]),
                                                            int(split[5]),
                                                            int(split[6]),
                                                            int(split[7]),
                                                            int(split[8]),
                                                            'cmap',
                                                            int(split[10]))
            else:
                warnings.warn("ReadError: found unsupported torsion-torsion type in: %s" % str(line[i]))
            if newTorsionTorsionForce:
                currentMoleculeType.torsiontorsionForceSet.add(newTorsionTorsionForce)

    def parse_exclusions(self, shared_args, sites_args):
        atomlist, entry_data, entry_values, ff_number, ff_type, vdwtypes, vdwtypeskeys = shared_args
        currentMolecule, i, lines, moleculeName, start, sysDirectiveAtm, verbose, currentMoleculeType = sites_args
        logger.debug("Parsing [ exclusions]...")
        for j in range(ff_number):
            temp = entry_values[j].split()
            temp.remove(temp[0])
            newExclusion = Exclusions(map(int,temp))  # convert to integers
            currentMoleculeType.exclusions.add(newExclusion)

    def parse_restraints(self, shared_args, sites_args):
        warnings.warn("Parsing [ restraints] not yet implemented")

    def parse_constraints(self, shared_args, sites_args):
        atomlist, entry_data, entry_values, ff_number, ff_type, vdwtypes, vdwtypeskeys = shared_args
        logger.debug("Parsing [ constraints]...")
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
            if 'HOH' in entry_values[j] or 'AH' in entry_values[j]:
                split = entry_values[j].split()
                tempatom = []
                templength = []
                for a in atompos:
                    if not '<>' in split[a]:
                        tempatom.append(int(split[a]))
                    else:
                        tempatom.append(None)
                for l in lenpos:
                    if not '<>' in split[l]:
                        if 'HOH' in entry_values[j] and l == lenpos[0]:
                            templength.append(float(split[l])*units.degrees) # Check units?
                        else:
                            templength.append(float(split[l])*units.angstroms) # Check units?
                    else:
                        templength.append(None*units.angstroms)
                if 'AH' in split[funct_pos]:
                    templen = int(list(split[funct_pos])[-1])
                elif 'HOH' in split[funct_pos]:
                    templen = 2    # Different desmond files have different options here.
                if templen == 1:
                    newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos])
                elif templen == 2:
                    newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos],tempatom[2],templength[1],None,templength[2])
                elif templen == 3:
                    newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos],tempatom[2],templength[1],tempatom[3],templength[2])
                elif templen == 4:
                    newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos],tempatom[2],templength[1],tempatom[3],templength[2],tempatom[4],templength[3])
                elif templen == 5:
                    newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos],tempatom[2],templength[1],tempatom[3],templength[2],tempatom[4],templength[3],tempatom[5],templength[4])
                elif templen == 6:
                    newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos],tempatom[2],templength[1],tempatom[3],templength[2],tempatom[4],templength[3],tempatom[5],templength[4],tempatom[6],templength[5])
                elif templen == 7:
                    newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos],tempatom[2],templength[1],tempatom[3],templength[2],tempatom[4],templength[3],tempatom[5],templength[4],tempatom[6],templength[5],tempatom[7],templength[6])
                elif templen == 8:
                    newConstraint = Constraint(tempatom[0],tempatom[1],templength[0],split[funct_pos],tempatom[2],templength[1],tempatom[3],templength[2],tempatom[4],templength[3],tempatom[5],templength[4],tempatom[6],templength[5],tempatom[7],templength[6],tempatom[8],templength[7])
            else:
                warnings.warn("ReadError: found unsupported constraint")
            if newConstraint:
                currentMoleculeType.constraints.add(newConstraint)

    def load_ffio_block(self, lines, moleculeName, start, end, sysDirective, sysDirectiveAtm, verbose = False):

#        Loading in ffio blocks from Desmond format
#        Args:
#            lines: list of all data in CMS format
#            moleculeName: name of current molecule
#            start: beginning of where ffio_ff starts for each molecule
#            end: ending of where ffio_ff ends for each molecule
#           sysDirective: help locate positions of specific data in ffio blocks
#           sysDirectiveAtm: help locate positions of specific data in m_atoms

        i = start
        j = start
        vdwtypeskeys = []
        vdwtypes = []
        sitetypes = []

        split = []
        constraints = []
        temp = []

        currentMolecule = None
        currentMoleculeType = None
        newAtomType = None
        newBondForce = None
        newBondType = None
        newPairForce = None
        newAngleForce = None
        newDihedralForce = None
        newTorsionTorsionForce = None
        stemp = None
        etemp = None
        sigma = None
        epsilon = None

        #There are several sections which require sites information to
        #process.  We keep a flag for this so that we are aware when
        #we have seen the sites
        bPreambleRead = False
        stored_ffio_types = []  # a list of stored ffio_type to keep track 
                              # of the ordering later
        atomlist = OrderedSet()
        namecol = 0
        combrcol = 0
        vdwtypercol = 0

        #DEFAULT VALUES WHEN CONVERTING TO GROMACS
        System._sys.nonbonded_function = 1
        System._sys.genpairs = 'yes'

        logger.debug('Parsing [ molecule %s]'%(moleculeName))
        logger.debug('Parsing [ ffio]')

        while i < end:
            if not bPreambleRead:
                # read the first section for the forces field info
                while not (':::' in lines[i]):
                    if 's_ffio_name' in lines[i]:
                        namecol = i-start-1
                    elif 's_ffio_comb_rule' in lines[i]:
                        combrcol = i-start-1
                    elif 's_ffio_vdw_func' in lines[i]:
                        vdwtypercol = i-start-1
                    i+=1
                i+=1 # skip the ':::'    
                # figure out combination rule
                combrule = lines[i+combrcol]
                if re.search("GEOMETRIC", combrule, re.IGNORECASE):
                    if re.search("ARITHMETIC", combrule, re.IGNORECASE):
                        System._sys.combination_rule = 'Lorentz-Berthelot'
                    else:
                        System._sys.combination_rule = 'Multiply-Sigeps'
                elif re.search("LJ12_6_C6C12", combrule, re.IGNORECASE):  # is this even valid?
                    System._sys.combination_rule = 'Multiply-C6C12'
                if (vdwtypercol > 0):
                    vdwrule = lines[i+vdwtypercol]
                # MISSING: need to identify vdw rule here -- currently assuming LJ12_6_sig_epsilon!

                # skip to the next ffio entry
                while not re.search('ffio',lines[i]):
                    i+=1
                bPreambleRead = True

            currentMolecule = Molecule(moleculeName)
            ff_type, ff_number, entry_data, entry_values, i = self.parse_ffio_block(lines,i,end)
            stored_ffio_types.append(ff_type)
            self.store_ffio_data(ff_type,ff_number, entry_data, entry_values)

        # Reorder so 'vdwtypes' is first, then 'sites'.  Could eventually get some simplification
        # by putting sites first, but too much rewriting for now.

        stored_ffio_types.insert(0, stored_ffio_types.pop(stored_ffio_types.index('ffio_sites')))
        stored_ffio_types.insert(0, stored_ffio_types.pop(stored_ffio_types.index('ffio_vdwtypes')))
    
        # now process all the data
        for type in stored_ffio_types:
            if type in sysDirective:
                ff_type, ff_number, entry_data, entry_values = self.retrive_ffio_data(type)
                shared_args = [atomlist, entry_data, entry_values, ff_number, ff_type, vdwtypes, vdwtypeskeys]
                sites_args = [currentMolecule, i, lines, moleculeName, start, sysDirectiveAtm, verbose, currentMoleculeType]
                if type == 'ffio_sites':
                    currentMoleculeType = sysDirective[type](shared_args, sites_args)
                else:
                    sysDirective[type](shared_args, sites_args)
            elif type == 'Done with ffio':
                continue    
            else:
                while '}' not in lines[i]:
                    i+=1   # not the most robust if there is nesting in a particular pattern

    def loadMBonds(self, atomlist, lines, start, end, npermol, verbose = False): #adds new bonds for each molecule in System

#        Loading in m_bonds in Desmond format
#        Args:
#            lines: list of all data in CMS format
#           start: beginning of where m_bonds starts for each molecule
#           end: ending of where m_bondsends for each molecule

        logger.debug("Parsing [ m_bonds]...")
        bg = False
        newBondForce = None
        split = []
        i = start
        bondForceSet = HashMap()
        forces = OrderedSet()
        while i < end:
            if ':::' in lines[i]:
                if bg:
                    break
                else:
                    bg = True
                    i+=1
            if bg:
                split = lines[i].split()
                atomi = int(split[1])
                atomj = int(split[2])
                bondingtypei = atomlist[atomi-1].name
                bondingtypej = atomlist[atomj-1].name
                params = [atomi, atomj, bondingtypei, bondingtypej]
                if atomi > npermol:  # we've collected the number of atoms per molecule.  Exit.
                    break
                order = int(split[3])
                kwd = [0, 0]
                optkwd = {'order': order, 'c': False}
                newBondForce = self.create_forcetype(HarmonicBond, params, kwd, optkwd)
                bondForceSet.add(newBondForce)
                forces.add(newBondForce)
            i+=1

        return [bondForceSet, forces]

    def loadMAtoms(self, lines, start, end, currentMolecule, slength, sysDirective, verbose = False): #adds positions and such to atoms in each molecule in System

#        Loading in m_atoms from Desmond format
#        Args:
#            lines: list of all data in CMS format
#           start: beginning of where m_atoms starts for each molecule
#           end: ending of where m_atoms ends for each molecule
#           currentMolecule
#           slength: number of unique atoms in m_atoms, used to calculate repetitions
#           sysDirective: help locate positions of specific data in m_atoms

        logger.debug("Parsing [ m_atom ]...")
        i = start
        xcol = None
        ycol = None
        zcol = None
        rincol = None
        rncol = None
        azcol = None
        pdbancol = None
        ancol = None
        vxcol = None
        vycol = None
        vzcol = None
        bg = False
        pdbaline = ""
        aline = ""

        mult = int(re.split('\W',lines[start].split()[0])[1])/slength

        while i < end:
            if ':::' in lines[i]:
                i+=1
                break
            else:
                match = sysDirective.match(lines[i])
                if match:
                    if match.group('First'):
                        start+=1
                    if match.group('xcoord'):
                        logger.debug("   Parsing [ xcoord]...")
                        xcol = i - start
                    elif match.group('ycoord'):
                        logger.debug("   Parsing [ ycoord]...")
                        ycol = i - start
                    elif match.group('zcoord'):
                        logger.debug("   Parsing [ zcoord]...")
                        zcol = i - start
                    elif match.group('rindex'):
                        logger.debug("   Parsing [ rindex]...")
                        rincol = i - start
                    elif match.group('rname'):
                        logger.debug("   Parsing [ rname]...")
                        rncol = i - start
                    elif match.group('aZ'):
                        logger.debug("   Parsing [ atomic number ]...")
                        azcol = i - start
                    elif match.group('pdbaname'):
                        logger.debug("   Parsing [ pdb atom name]...")
                        pdbancol = i - start
                    elif match.group('aname'):
                        logger.debug("   Parsing [ atom name]...")
                        ancol = i - start
                    elif match.group('xvelocity'):
                        logger.debug("   Parsing [ xvelocity]...")
                        vxcol = i - start
                    elif match.group('yvelocity'):
                        logger.debug("   Parsing [ yvelocity]...")
                        vycol = i - start
                    elif match.group('zvelocity'):
                        logger.debug("   Parsing [ zvelocity]...")
                        vzcol = i - start
            i+=1

        atom = None

        newMoleculeAtoms = []
        j = 0
        logger.debug("   Parsing atoms...")

        molecules = []    
        while j < mult:
            newMolecule = copy.deepcopy(currentMolecule)
            for atom in newMolecule.atoms:
                if ':::' in lines[i]:
                    break
                else:
                    aline = split_with_quotes(lines[i])
                    atom.residue_index = int(aline[rincol])
                    atom.residue_name = aline[rncol].strip()
                    try:
                        atom.atomic_number = int(aline[azcol])
                    except Exception as e:
                        logger.exception(e) # EDZ: just pass statement before, now exception is recorded, but supressed
                    atom.position = [float(aline[xcol]) * units.angstroms,
                                     float(aline[ycol]) * units.angstroms,
                                     float(aline[zcol]) * units.angstroms]
                    if vxcol == vycol == vzcol == None:
                        atom.velocity = [0.0 * units.angstroms * units.picoseconds**(-1),
                                        0.0 * units.angstroms * units.picoseconds**(-1),
                                        0.0 * units.angstroms * units.picoseconds**(-1)]
                    else:
                        atom.velocity = [float(aline[vxcol]) * units.angstroms * units.picoseconds**(-1),
                                        float(aline[vycol]) * units.angstroms * units.picoseconds**(-1),
                                        float(aline[vzcol]) * units.angstroms * units.picoseconds**(-1)]
                    if (pdbancol):    
                        pdbaline = aline[pdbancol].strip()
                    if (ancol):
                        aline = aline[ancol].strip()
                    if re.match('$^',pdbaline) and not re.match('$^',aline):
                        atom.name = aline
                    elif re.match('$^',aline) and not re.match('$^',pdbaline):
                        atom.name = pdbaline
                    elif re.search("\d+",pdbaline) and not re.search("\d+",aline):
                        if re.search("\D+",pdbaline) and re.search("\w+",pdbaline):
                            atom.name = pdbaline
                        else:
                            atom.name = aline
                    elif re.search("\d+",aline) and not re.search("\d+",pdbaline):
                        if re.search("\D+",aline) and re.search("\w+",aline):
                            atom.name = aline
                        else:
                            atom.name = pdbaline
                    elif re.match('$^',pdbaline) and re.match('$^',aline):
                        atom.name = "None"
                    else:
                        atom.name = aline  #doesn't matter which we choose, so we'll go with atom name instead of pdb
                    i+=1

            molecules.append(newMolecule)        
            j+=1

        return molecules


    def load_box_vector(self, lines, start, end, verbose = False):

#       Loading Box Vector
#       Create a Box Vector to load into the System
#        Args:
#            lines: all the lines of the file stored in an array
#            start: starting position
#            end: ending position

        i = start
        v = np.zeros([3,3])*units.angstroms
        while (i<end):
            if 'r_chorus_box_ax' in lines[i]:
                startboxlabel = i-start
            if ':::' in lines[i]:
                endlabel = i
                break
            i+=1
        startbox = startboxlabel+endlabel
        nvec = 0
        for i in range(startbox,startbox+9):
            j = (nvec)/3
            k = (nvec)%3
            v[j,k] = float(re.sub(r'\s', '', lines[i])) * units.angstrom
            nvec += 1

        System._sys.box_vector = v

    def read_file(self, filename, verbose=True):

#        Load in data from file

#       Read data in Desmond format

#        Args:
#            filename: the name of the file to write out to

        lines = list()

        molnames = []
        fl = open(filename, 'r')
        lines = list(fl)
        fl.close()
        i,j=0,0

        for line in lines:
            if re.search("f_m_ct",line,re.VERBOSE):
                if j > 0:
                    self.fblockpos.append(i)
                j+=1
            if re.search("m_atom",line,re.VERBOSE) and not (re.search("i_m",line)
            or re.search("s_m",line)):
                if j > 1:
                    self.a_blockpos.append(i)
                j+=1
            if re.search("m_bond",line,re.VERBOSE):
                if j > 2:
                    self.b_blockpos.append(i)
                j+=1
            if re.search("ffio_ff",line,re.VERBOSE):
                if j > 2:
                    self.ffio_blockpos.append(i)
                j+=1
            i+=1
        i-=1
        self.fblockpos.append(i)
        self.a_blockpos.append(i)
        self.b_blockpos.append(i)
        self.ffio_blockpos.append(i)

        sysDirectiveTopD = {'ffio_vdwtypes': self.parse_vdwtypes,
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

        sysDirectiveTop = re.compile(r"""
          ((?P<vdwtypes>\s*ffio_vdwtypes)
          |
          (?P<sites>\s*ffio_sites)
          |
          (?P<bonds>\s*ffio_bonds)
          |
          (?P<pairs>\s*ffio_pairs)
          |
          (?P<angles>\s*ffio_angles)
          |
          (?P<dihedrals>\s*ffio_dihedrals)
          |
          (?P<torsiontorsion>\s*ffio_torsion_torsion)
          |
          (?P<constraints>\s*ffio_constraints)
          |
          (?P<exclusions>\s*ffio_exclusions)
          |
          (?P<restraints>\s*ffio_restraints))
        """, re.VERBOSE)

        sysDirectiveStr = re.compile(r"""
          ((?P<xcoord>\s*r_m_x_coord)
          |
          (?P<ycoord>\s*r_m_y_coord)
          |
          (?P<zcoord>\s*r_m_z_coord)
          |
          (?P<rindex>\s*i_m_residue_number)
          |
          (?P<rname>\s*s_m_pdb_residue_name)
          |
          (?P<aZ>\s*i_m_atomic_number)
          |
          (?P<pdbaname>\s*s_m_pdb_atom_name)
          |
          (?P<aname>\s*s_m_atom_name)
          |
          (?P<xvelocity>\s*r_ffio_x_vel)
          |
          (?P<yvelocity>\s*r_ffio_y_vel)
          |
          (?P<zvelocity>\s*r_ffio_z_vel)
          |
          (?P<First>\s*[#][\s+]First[\s+]column[\s+]is[\s+]atom[\s+]index[\s+][#]))
        """, re.VERBOSE)

        #LOADING Ffio blocks
        logger.debug("Reading ffio block...")
        #MRS: warning -- currently no check to avoid duplicated molecule names. Investigate.
        i = 0
        j = 0
        while i < (len(self.ffio_blockpos)-1):
            j = self.fblockpos[i]
            while not re.match(r'\s*[:::]',lines[j]):
                j+=1
            # make sure we have reasonable molecular names. 
            molname = lines[j+1].strip()
            molname = molname.replace("\"","")  # get rid of quotation marks so we can find unique names
            if molname == "":
                molname = "Molecule_"+str(len(molnames)+1)
            molnames.append(molname)    
            self.load_ffio_block(lines, molname, self.ffio_blockpos[i], self.fblockpos[i+1]-1, sysDirectiveTopD, sysDirectiveStr,  verbose)
            i+=1
        i = 0

        #LOAD RAW BOX VECTOR-Same throughout cms

        logger.debug("Reading Box Vector...")
        self.load_box_vector(lines, self.fblockpos[0], self.a_blockpos[0], verbose)

    def write_vdwtypes_and_sites(self, molecule):

        #-ADDING VDWTYPES AND SITES
        i = 0
        vdwtypes = []
        sites = []
        sig = None
        ep = None
        stemp = None
        etemp = None
        combRule = System._sys.combination_rule
        for atom in molecule.atoms:
            i+=1
            if atom.residue_index:
                sites.append(' %3d %5s %9.8f %9.8f %2s %1d %4s\n' % (
                        i, 'atom',
                        atom._charge[0].in_units_of(units.elementary_charge)._value,
                        atom._mass[0].in_units_of(units.atomic_mass_unit)._value,
                        atom.atomtype[0], atom.residue_index, atom.residue_name))
            else:
                sites.append(' %3d %5s %9.8f %9.8f %2s\n' % (
                        i, 'atom',
                        atom._charge[0].in_units_of(units.elementary_charge)._value,
                        atom._mass[0].in_units_of(units.atomic_mass_unit)._value,
                        atom.atomtype[0]))

            sig = float(atom.sigma[0].in_units_of(units.angstroms)._value)
            ep = float(atom.epsilon[0].in_units_of(units.kilocalorie_per_mole)._value)
            if combRule == 'Multiply-C6C12':   #MRS: seems like this should be automated more?
                stemp = ep * (4 * (sig**6))
                etemp = stemp * (sig**6)
            elif combRule in ['Lorentz-Berthelot','Multiply-Sigeps']:
                stemp = sig
                etemp = ep
            if ' %2s %18s %8.8f %8.8f\n' % (atom.atomtype[0], "LJ12_6_sig_epsilon", float(stemp), float(etemp)) not in vdwtypes:
                vdwtypes.append(' %2s %18s %8.8f %8.8f\n' % (atom.atomtype[0], "LJ12_6_sig_epsilon", float(stemp), float(etemp)))

        lines = []
        logger.debug("   -Writing vdwtypes...")
        lines.append("    ffio_vdwtypes[%d] {\n"%(len(vdwtypes)))
        lines.append("      s_ffio_name\n")
        lines.append("      s_ffio_funct\n")
        lines.append("      r_ffio_c1\n")
        lines.append("      r_ffio_c2\n")
        lines.append("      :::\n")
        i = 0
        for v in vdwtypes:
            i+=1
            lines.append('      %d%2s'%(i,v))
        lines.append("      :::\n")
        lines.append("    }\n")

        logger.debug("   -Writing sites...")
        lines.append("    ffio_sites[%d] {\n"%(len(sites)))
        lines.append("      s_ffio_type\n")
        lines.append("      r_ffio_charge\n")
        lines.append("      r_ffio_mass\n")
        lines.append("      s_ffio_vdwtype\n")
        if len(sites[0].split()) > 5:   # fix this to explicitly ask if resnr is in here rather than length
            lines.append("      i_ffio_resnr\n")
            lines.append("      s_ffio_residue\n")
        lines.append("      :::\n")
        for s in sites:
            lines.append('   %s'%(s))

        lines.append("      :::\n")
        lines.append("    }\n")

        return lines

    def write_bonds(self, moleculetype):

        #-ADDING BONDS
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
        bondlist = sorted(moleculetype.bondForceSet.itervalues(), key=lambda x: (x.atom1,x.atom2))
        for bond in bondlist:
            i += 1
            kwds = self.get_parameter_kwds_from_force(bond)
            try:
                name = self.lookup_desmond_bonds[bond.__class__]
            except:
                raise Exception("WriteError: bondtype %s is not supported by Desmond" % (bond.__class__.__name__))    
            name, kwds = self.canonical_bond(name, kwds, bond.__class__, direction = 'from')
            line = '      %d %d %d %s' %(i, bond.atom1, bond.atom2, name)
            bond_params = self.get_parameter_list_from_kwds(bond, kwds)
            u = self.unitvars[bond.__class__.__name__]
            for j, p in enumerate(bond_params):
                line += "%15.8f" % (p.value_in_unit(u[j]))
            line += '\n'    
            dlines.append(line)
        header = "    ffio_bonds[%d] {\n" % (i)
        hlines = endheadersection(i==0,header,hlines)

        dlines.append("      :::\n")
        dlines.append("    }\n")
        hlines.extend(dlines)
        return hlines

    def write_angles(self, moleculetype):
        #-ADDING ANGLES
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
        anglelist = sorted(moleculetype.angleForceSet.itervalues(), key=lambda x: (x.atom1,x.atom2,x.atom3))
        for angle in anglelist:
            i+=1
            kwds = self.get_parameter_kwds_from_force(angle)
            try:
                name = self.lookup_desmond_bonds[angle.__class__]
            except:
                raise Exception("WriteError: angletype %s is not supported by Desmond" % (angle.__class__.__name__))
            name, kwds = self.canonical_angle(name, kwds, angle.__class__, direction = 'from')
            line = '      %d %d %d %d %s' %(i, angle.atom1, angle.atom3, angle.atom3, name)
            angle_params = self.get_parameter_list_from_kwds(angle, kwds)
            u = self.unitvars[angle.__class__.__name__]
            for j, p in enumerate(angle_params):
                line += "%15.8f" % (p.value_in_unit(u[j]))
            line += '\n'    
            dlines.append(line)
            
        header = "    ffio_angles[%d] {\n" % (i)
        hlines = endheadersection(i==0,header,hlines)

        dlines.append("      :::\n")
        dlines.append("    }\n")
        hlines.extend(dlines)
        return hlines

    def write_dihedrals(self, moleculetype):
        #-ADDING DIHEDRALS
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
            hlines.append("      r_ffio_c%d\n" %(ih))
        hlines.append("      :::\n")

        i = 0
        #sorting by first index
        dihedrallist = sorted(moleculetype.dihedralForceSet.itervalues(), key=lambda x: x.atom1)
        # first, identify the number of terms we will print
        for dihedral in dihedrallist:
            i+=1
            dlines.append('      %d %d %d %d %d ' % (
                    i, dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4))
            if isinstance(dihedral, ImproperHarmonicDihedral):
                dlines.append('%s %10.8f %10.8f %1d %1d %1d %1d %1d %1d\n' % (
                        'Improper_Harm', float(dihedral.xi.in_units_of(units.degrees)._value),
                        0.5*float(dihedral.k.in_units_of(units.kilocalorie_per_mole/units.radians**2)._value),
                        0,0,0,0,0,0))
            elif isinstance(dihedral, TrigDihedral):
                if dihedral.improper:
                    dtype = 'Improper_Trig'
                else:
                    dtype = 'Proper_Trig'
                dlines.append('%s %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f\n' % (
                        dtype,
                        float(dihedral.phi.in_units_of(units.degrees)._value),
                        float(dihedral.fc0.in_units_of(units.kilocalories_per_mole)._value),
                        float(dihedral.fc1.in_units_of(units.kilocalories_per_mole)._value),
                        float(dihedral.fc2.in_units_of(units.kilocalories_per_mole)._value),
                        float(dihedral.fc3.in_units_of(units.kilocalories_per_mole)._value),
                        float(dihedral.fc4.in_units_of(units.kilocalories_per_mole)._value),
                        float(dihedral.fc5.in_units_of(units.kilocalories_per_mole)._value),
                        float(dihedral.fc6.in_units_of(units.kilocalories_per_mole)._value)))
            else:
                logger.error("WriteError: found unsupported dihedral. \n    %s\n    Printing zero-energy placeholder."
                             % str(dihedral))
                dlines.append('%s 0.0 %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f\n' %
                              (i, dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4, 'Proper_Trig',
                               0, 0, 0, 0, 0, 0, 0))
        header = "    ffio_dihedrals[%d] {\n" % (i)
        hlines = endheadersection(i==0,header,hlines)

        dlines.append("      :::\n")
        dlines.append("    }\n")
        hlines.extend(dlines)
        return hlines

    def write_torsion_torsion(self, moleculetype):

        # adding TORSION-TORSION terms
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

        for torsiontorsion in moleculetype.torsiontorsionForceSet.itervalues():
            i+=1
            # only type of torsion/torsion is CMAP currently
            dlines.append('      %d %d %d %d %d %d %d %d %d %s %d\n' % (
                    i,
                    int(torsiontorsion.atom1), int(torsiontorsion.atom2),
                    int(torsiontorsion.atom3), int(torsiontorsion.atom4),
                    int(torsiontorsion.atom5), int(torsiontorsion.atom6),
                    int(torsiontorsion.atom7), int(torsiontorsion.atom8),
                    'cmap', torsiontorsion.chart))
        header = "    ffio_torsion_torsion[%d] {\n"%(i)

        hlines = endheadersection(i==0,header,hlines)
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
    #ADDING EXCLUSIONS
        i = 0
        logger.debug("   -Writing exclusions...")
        hlines = list()
        dlines = list()
        hlines.append("    ffio_exclusions_placeholder\n")
        hlines.append("      i_ffio_ai\n")
        hlines.append("      i_ffio_aj\n")
        hlines.append("      :::\n")

        if moleculetype.nrexcl == 0:

            # Should probably be determined entirely by the bonds,
            # since settles now adds bonds.  For now, leave this in
            # for Desmond to Desmond conversion, where nrexcl is not
            # determined.  Probably should switch eventually.

            for exclusion in moleculetype.exclusions.itervalues():
                i+=1
                dlines.append('      %d %d %d\n'%(i, int(exclusion.exclusions[0]), int(exclusion.exclusions[1])))

        else:
            if moleculetype.nrexcl > 4:
                warnings.warn("Can't handle more than excluding 1-4 interactions right now!")

            fullbondlist = []
            fullbondlist = sorted(moleculetype.bondForceSet.itervalues(), key=lambda x: (x.atom1,x.atom2))
            # exclude HarmonicPotential types, which do not have exclusions.
            bondlist = [bond for bond in fullbondlist if (not isinstance(bond, HarmonicPotentialBond))]

            # first, figure out the first appearance of each atom in the bondlist
            currentatom = 0
            atompos = []
            bondindex = 0
            nsize = len((moleculetype.moleculeSet)[0].atoms)+1
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
                        dlines.append('      %d %d %d\n' % (i, atom, a))


        header = "    ffio_exclusions[%d] {\n"%(i)
        hlines = endheadersection(i==0,header,hlines)

        dlines.append("      :::\n")
        dlines.append("    }\n")
        hlines.extend(dlines)
        return hlines

    def write_pairs(self, moleculetype):

        #-ADDING PAIRS
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
        for pair in moleculetype.pairForceSet.itervalues():
            i += 2
            atom = '      %d %d ' % (pair.atom1, pair.atom2)
            if isinstance(pair,LjDefaultPair) or isinstance(pair,LjqDefaultPair):
                dlines += '       %d %d %d LJ %10.8f <>\n' % (i-1, pair.atom1, pair.atom2, System._sys.lj_correction)
                dlines += '       %d %d %d Coulomb %10.8f <>\n' % (i, pair.atom1, pair.atom2, System._sys.coulomb_correction)
            elif isinstance(pair, LjSigepsPair) or isinstance(pair, LjCPair) or isinstance(pair, LjqCPair) or isinstance(pair, LjqSigepsPair):
                # Check logic here -- not clear that we can correctly determine which type it is.
                # Basically, I think it's whether scaleLJ is defined or not.
                if pair.scaleLJ:
                    dlines += '       %d %d %d LJ %10.8f <>\n' % (i-1, pair.atom1, pair.atom2, pair.scaleLJ)
                    i -= 1
                else:
                    if isinstance(pair, LjCPair) or isinstance(pair, LjqCPair):
                        #if isinstance(pair, LjqCPair):
                        epsilon = 0.25 * (pair.C6**2) / pair.C12    # (16*eps^2*sig^12 / 4 eps*sig^12) = 4 eps
                        sigma = (0.25 * pair.C6 / epsilon)**(1.0/6.0)  # (0.25 * 4 eps sig^6 / eps)^(1/6)
                        #else:
                        #    epsilon = pair.epsilon
                        #    sigma = pair.sigma

                        dlines += '       %d %d %d LJ12_6_sig_epsilon %10.8f %10.8f\n' % (i-1,
                                                                          pair.atom1,
                                                                          pair.atom2,
                                                                          sigma.value_in_unit(units.angstroms),
                                                                          epsilon.value_in_unit(units.kilocalorie_per_mole))
                        i -= 1

                if isinstance(pair, LjSigepsPair) or isinstance(pair, LjqCPair):
                    if pair.scaleQQ:
                        dlines += '       %d %d %d Coulomb %10.8f <>\n' % (i-1, pair.atom1, pair.atom2,
                                                                            pair.scaleQQ.value_in_unit(units.dimensionless))
                        i -= 1
                    else:
                        dlines += '       %d %d %d Coulomb %10.8f <>\n' % (i-1, pair.atom1, pair.atom2,
                                                                                System._sys.coulomb_correction)
                        i -= 1

            else:
                warnings.warn("Unknown pair type!")


        header = "    ffio_pairs[%d] {\n"%(i)
        hlines = endheadersection(i==0,header,hlines)

        dlines.append("      :::\n")
        dlines.append("    }\n")
        hlines.extend(dlines)
        return hlines

    def write_constraints(self, moleculetype):
        #ADDING CONSTRAINTS
        logger.debug("   -Writing constraints...")
        isHOH = False

        if (moleculetype.settles):
            alen = 2
            clen = 3
        else:
            alen = 0
            clen = 0

        alen_max = alen
        clen_max = clen

        for constraint in moleculetype.constraints.itervalues():
            if re.search('AH',constraint.type):
                alen = int(list(constraint.type)[-1])
                clen = alen
            elif re.match('HOH',constraint.type):
                alen = 2
                clen = 3
            if alen_max < alen:
                alen_max = alen
                clen_max = clen
        # we now know the maximum length of all constraint types

        # not sure we need to sort these, but makes it easier to debug
        i = 0
        constraintlist = sorted(moleculetype.constraints.itervalues(),key=lambda x: x.atom1)
        dlines = list()
        hlines = list()

        for constraint in constraintlist: #calculate the max number of atoms in constraint
            i+=1
            if re.search('HOH',constraint.type):
                cline = '      %d %d %d %d ' % (i,int(constraint.atom1),int(constraint.atom2),int(constraint.atom3))
                for j in range(alen_max-3):
                    cline += '0 '
                cline += constraint.type
                cline += ' %10.8f' % (float(constraint.length1.in_units_of(units.degrees)._value))
                cline += ' %10.8f' % (float(constraint.length2.in_units_of(units.angstroms)._value))
                cline += ' %10.8f' % (float(constraint.length2.in_units_of(units.angstroms)._value))
                for j in range(clen_max-3):
                    cline += ' <>'
            elif re.match('AH',constraint.type):
                alen = int(list(constraint.type)[-1])
                cline = '      %d ' % i
                cline += ' %d ' % int(constraint.atom1)
                cline += ' %d ' % int(constraint.atom2)
                if alen > 1:
                    cline += ' %d ' % int(constraint.atom3)
                if alen > 2:
                    cline += ' %d ' % int(constraint.atom4)
                if alen > 3:
                    cline += ' %d ' % int(constraint.atom5)
                if alen > 4:
                    cline += ' %d ' % int(constraint.atom6)
                if alen > 5:
                    cline += ' %d ' % int(constraint.atom7)
                if alen > 6:
                    cline += ' %d ' % int(constraint.atom8)
                for j in range(alen,alen_max):
                    cline += ' 0 '
                cline += constraint.type
                cline += ' %10.8f' % (float(constraint.length1.in_units_of(units.angstroms)._value))
                if alen > 1:
                    cline += ' %10.8f' % (float(constraint.length2.in_units_of(units.angstroms)._value))
                if alen > 2:
                    cline += ' %10.8f' % (float(constraint.length3.in_units_of(units.angstroms)._value))
                if alen > 3:
                    cline += ' %10.8f' % (float(constraint.length4.in_units_of(units.angstroms)._value))
                if alen > 4:
                    cline += ' %10.8f' % (float(constraint.length5.in_units_of(units.angstroms)._value))
                if alen > 5:
                    cline += ' %10.8f' % (float(constraint.length6.in_units_of(units.angstroms)._value))
                if alen > 6:
                    cline += ' %10.8f' % (float(constraint.length7.in_units_of(units.angstroms)._value))
                for j in range(alen,alen_max):
                    cline += ' 0.0'
            cline += '\n'
            dlines.append(cline)

        # now need to add the constraints specified through settles.  Only one settles per molecule
        if (moleculetype.settles):
            i += 1

            settles = moleculetype.settles
            # Assumes the water arrangement O, H, H, which might not always be the case.  Consider adding detection.
            cline = '      %d %d %d %d ' % (i,1,3,2)
            for j in range(alen_max-3):
                cline += '0 '
            cline += ' HOH '
            dOH = settles.dOH._value
            dHH = settles.dHH._value
            angle = 2.0*math.asin(0.5*dHH/dOH)*(180/math.pi)    # could automate conversion. . .
            cline += " %.3f %.5f %.5f " % (angle,dOH,dOH)
            cline += '\n'
            for j in range(alen,alen_max):
                cline += ' 0.0'
            dlines.append(cline)

        hlines.append("    ffio_constraints[%d] {\n"%(i))
        if (i==0):
            hlines.append("      :::\n")
        else:
            letters = ['i','j','k','l','m','n','o','p','q']
            for j in range(alen_max+1):
                hlines.append('      i_ffio_a%s\n'%letters[j])
            hlines.append('      s_ffio_funct\n')
            for j in range(clen_max):
                hlines.append('      r_ffio_c%d\n' %(j+1))
            hlines.append("      :::\n")

        dlines.append("      :::\n")
        dlines.append("    }\n")
        hlines.extend(dlines)
        return hlines

    def write_file(self, filename, verbose=True):

#        Write this topology to file
#        Write out this topology in Desmond format
#        Args:
#            filename: the name of the file to write out to

        lines = list()
        vdwtypes = []
        sites = []
        pos = 0
        name = ''
        sig = None
        ep = None
        stemp = None
        etemp = None

        logger.warn("MacroModel atom type is not defined in other files, is set to 1 for all cases as it must be validly defined for desmond files to run.  However, it does not affect the energies.")

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
        lines.append('  r_chorus_box_ax\n')
        lines.append('  r_chorus_box_ay\n')
        lines.append('  r_chorus_box_az\n')
        lines.append('  r_chorus_box_bx\n')
        lines.append('  r_chorus_box_by\n')
        lines.append('  r_chorus_box_bz\n')
        lines.append('  r_chorus_box_cx\n')
        lines.append('  r_chorus_box_cy\n')
        lines.append('  r_chorus_box_cz\n')
        lines.append('  s_ffio_ct_type\n')
        lines.append('  :::\n')

        #box vector
        bv = System._sys.box_vector
        lines.append('  "full system"\n')
        for bi in range(3):
            for bj in range(3):
                lines.append('%22s\n'%float(bv[bi][bj].in_units_of(units.angstroms)._value))
        lines.append('  full_system\n')

        #M_ATOM
        apos = len(lines) #pos of where m_atom will be; will need to overwite later based on the number of atoms
        lines.append('m_atom\n')
        lines.append('    # First column is atom index #\n')
        lines.append('    i_m_mmod_type\n')
        lines.append('    r_m_x_coord\n')
        lines.append('    r_m_y_coord\n')
        lines.append('    r_m_z_coord\n')
        lines.append('    i_m_residue_number\n')
        lines.append('    s_m_pdb_residue_name\n')
        lines.append('    i_m_atomic_number\n')
        lines.append('    s_m_atom_name\n')
        lines.append('    r_ffio_x_vel\n')
        lines.append('    r_ffio_y_vel\n')
        lines.append('    r_ffio_z_vel\n')
        lines.append('    :::\n')

        i = 0
        nmol = 0
        totalatoms = []
        totalatoms.append(0)
        for moleculetype in System._sys._molecules.itervalues():
            for molecule in moleculetype.moleculeSet:
                for atom in molecule.atoms:
                    i += 1
                    line = '    %d        %d' % (i,1) #HAVE TO PUT THE 1 HERE OR ELSE DESMOND DIES, EVEN THOUGH IT DOESN'T USE IT
                    for j in range(3):
                        line += " %10.8f" % (float(atom._position[j].value_in_unit(units.angstroms)))
                    line +=   "     %2d %4s    %2d  %2s" % (
                        atom.residue_index,
                        '"%s"'%atom.residue_name,
                        atom.atomic_number,
                        '"%s"'%atom.name)
                    if atom._velocity:
                        for j in range(3):
                            line += " %10.8f" % (float(atom._velocity[j].value_in_unit(units.angstroms / units.picoseconds)))
                    else:
                         for j in range(3):
                            line += " %10.8f" % (0)
                    lines.append(line + '\n')
            totalatoms.append(i)

        lines[apos] = '  m_atom[%d] {\n'%(i)
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
        for moleculetype in System._sys._molecules.itervalues():
            # sort the bondlist because Desmond requires the first time a bond is listed to have
            # the atoms in ascending order
            repeatmol = len(moleculetype.moleculeSet)
            #MRS: need to be fixed; gromacs loads in one set of bonds per molecue; desmond loads in all
            atompermol = len(moleculetype.moleculeSet[0].atoms)
            bondlist = sorted(moleculetype.bondForceSet.itervalues(), key=lambda x: (x.atom1,x.atom2))
            for n in range(repeatmol):
                for bond in bondlist:
                    if bond and bond.order:
                        i += 1
                        dlines.append('    %d %d %d %d %d %d\n'
                                      %(i,
                                        bond.atom1 + n*atompermol + totalatoms[nmol],
                                        bond.atom2 + n*atompermol + totalatoms[nmol],
                                        int(bond.order),
                                        1,
                                        1))
                    elif not bond:
                        nonecnt+=1
                if nonecnt > 0:
                    logger.debug('FOUND %d BONDS THAT DO NOT EXIST' % nonecnt)
            nmol +=1

        hlines[0] = '  m_bond[%d] {\n' % i    
        if (i > 0):
            lines.extend(hlines)
            lines.extend(dlines)
            lines.append('    :::\n')
            lines.append('  }\n')
            lines.append('}\n')

        solute = True
        endline = ''
        resName = ''

        #WRITE OUT ALL FFIO AND F_M_CT BLOCKS

        for moleculetype in System._sys._molecules.itervalues():
            logger.debug('Writing molecule block %s...'% moleculetype.name)
            #BEGINNING BLOCK
            
            logger.debug("  Writing f_m_ct...")
            lines.append('f_m_ct {\n')
            lines.append('  s_m_title\n')
            bpos = len(lines) #bpos temporarily used for position of s_m_entry_name (for TIP3)
            lines.append('  s_m_entry_name\n')
            lines.append('  i_ffio_num_component\n')
            lines.append('  r_chorus_box_ax\n')
            lines.append('  r_chorus_box_ay\n')
            lines.append('  r_chorus_box_az\n')
            lines.append('  r_chorus_box_bx\n')
            lines.append('  r_chorus_box_by\n')
            lines.append('  r_chorus_box_bz\n')
            lines.append('  r_chorus_box_cx\n')
            lines.append('  r_chorus_box_cy\n')
            lines.append('  r_chorus_box_cz\n')
            lines.append('  s_ffio_ct_type\n')
            lines.append('  :::\n')

            if solute:
                lines.append('  solute\n')
                endline = '  solute\n'
                solute = False
                del lines[bpos]
                del lines[bpos]
            else:
                for atom in molecule.atoms:
                    resName = atom.residue_name
                    break
                if re.match("T3P", resName) or re.search("WAT", resName):
                    #lines[bpos] = ('  s_m_entry_name\n')
                    lines.append('  "TIP3P water box"\n')
                    lines.append('  "TIP3P water box"\n')
                    lines.append('  1\n')
                    endline = '  solvent\n'
                else:
                    lines.append('  %s\n'%(moleculetype.name))
                    endline = '  ion\n'
                    del lines[bpos]
                    del lines[bpos] #deletes line for num component (only in TIP3)

            for bi in range(3):
                for bj in range(3):
                    lines.append('%22s\n'%float(bv[bi][bj].in_units_of(units.angstroms)._value))
            lines.append(endline)

            #M_ATOMS

            logger.debug("  Writing m_atoms...")
            apos = len(lines) #pos of where m_atom will be; will need to overwite later based on the number of atoms
            lines.append('m_atom\n')
            lines.append('    # First column is atom index #\n')
            lines.append('    i_m_mmod_type\n')
            lines.append('    r_m_x_coord\n')
            lines.append('    r_m_y_coord\n')
            lines.append('    r_m_z_coord\n')
            lines.append('    i_m_residue_number\n')
            lines.append('    s_m_pdb_residue_name\n')
            lines.append('    i_m_atomic_number\n')
            lines.append('    s_m_atom_name\n')
            lines.append('    r_ffio_x_vel\n')
            lines.append('    r_ffio_y_vel\n')
            lines.append('    r_ffio_z_vel\n')
            lines.append('    :::\n')

            i = 0
            for molecule in moleculetype.moleculeSet:
                for atom in molecule.atoms:
                    i += 1
                    #NOT SURE WHAT TO PUT FOR MMOD TYPE; 1 is currently used.
                    #This can't be determined currently from the information provided,
                    # unless it is stored previous, nor is it used by desmond
                    line = '    %d        %d' % (i,1)
                    for j in range(3):
                        line += " %10.8f" % (float(atom._position[j].value_in_unit(units.angstroms)))
                    line +=   "     %2d %4s    %2d  %2s" % (
                                atom.residue_index,
                                '"%s"'%atom.residue_name,
                                atom.atomic_number,
                                '"%s"'%atom.name)
                    if atom._velocity:
                        for j in range(3):
                            line += " %10.8f" % (float(atom._velocity[j].value_in_unit(units.angstroms / units.picoseconds)))
                    else:
                         for j in range(3):
                            line += " %10.8f" % (0)
                    lines.append(line + '\n')

            lines[apos] = '  m_atom[%d] {\n'%(i)
            lines.append('    :::\n')
            lines.append('  }\n')

            #M_BONDS
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

            repeatmol = len(moleculetype.moleculeSet)
            atompermol = len(moleculetype.moleculeSet[0].atoms)
            bondlist = sorted(moleculetype.bondForceSet.itervalues(), key=lambda x: x.atom1)
            for n in range(repeatmol): 
                for bond in bondlist:
                    if bond and bond.order:
                        i += 1
                        dlines.append('    %d %d %d %d %d %d\n'
                                      %(i,
                                        bond.atom1 + n*atompermol,
                                        bond.atom2 + n*atompermol,
                                        int(bond.order),
                                        1,
                                        1))
                    else:
                        nonecnt+=1
                if nonecnt > 0:
                    logger.debug('FOUND %d BONDS THAT DO NOT EXIST' % nonecnt)

            header = '  m_bond[%d] {\n'%i

            if (i>0):
                hlines = endheadersection(False,header,hlines)
                lines.extend(hlines)
                lines.extend(dlines)
                lines.append('    :::\n')
                lines.append('  }\n')

            #FFIO
            molecule =  moleculetype.moleculeSet[0]
            logger.debug("  Writing ffio...")
            lines.append('  ffio_ff {\n')
            lines.append('    s_ffio_name\n')
            lines.append('    s_ffio_comb_rule\n')
            lines.append('    i_ffio_version\n')
            lines.append('    :::\n')

            #Adding Molecule Name
            if re.search("Viparr", moleculetype.name):
                lines.append('    Generated by Viparr\n')
            else:
                lines.append('    %s\n' % moleculetype.name)

            #Adding Combination Rule
            if System._sys.combination_rule == 'Multiply-C6C12':
                lines.append('    C6C12\n')   # this may not exist in DESMOND, or if so, need to be corrected
            elif System._sys.combination_rule == 'Lorentz-Berthelot':
                lines.append('    ARITHMETIC/GEOMETRIC\n')
            elif System._sys.combination_rule == 'Multiply-Sigeps':
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

        fout = open(filename, 'w')
        for line in lines:
            fout.write(line)
        fout.close()
