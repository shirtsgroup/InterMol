import argparse
from collections import OrderedDict
import logging
import os
import sys
import warnings

import numpy as np
import parmed as pmd
from parmed.utils.six import string_types
import parmed.unit as units

import intermol.gromacs as gmx
import intermol.lammps as lmp
import intermol.desmond as des
import intermol.amber as amb
import intermol.charmm as crm
import intermol.tests


# Make a global logging object.
logger = logging.getLogger('InterMolLog')
logger.setLevel(logging.DEBUG)
logging.captureWarnings(True)
warning_logger = logging.getLogger('py.warnings')

# we convert them into these names 
canonical_energy_names = [
    'bond',

    'angle', 'urey-bradley',

    'dihedral', 'improper', 'proper', 'cmap',

    'h-bond',
    'vdw total', 'vdw-14', 'vdw (SR)', 'vdw (LR)',
    'coulomb total', 'coulomb-14', 'coulomb (SR)', 'coulomb (LR)',

    'nonbonded', 'bonded',

    'potential']

# these names we actually expect to be similar.
comparable_energy_terms = ['bond','angle','dihedral','bonded','nonbonded','potential']

def canonicalize_energy_names(energy_dict, canonical_keys):
    """Adjust the keys in energy_dict to the canonical names.

    Parameters
    ----------
    energy_dict : OrderedDict
    engine : str

    Returns
    -------
    normalized : OrderedDict

    """
    # TODO: Look into creating an `EnergyDict` class.
    normalized = OrderedDict.fromkeys(canonical_energy_names,
                                      0 * units.kilojoules_per_mole)

    for key, energy in energy_dict.items():
        canonical_key = canonical_keys.get(key)
        if canonical_key is None:
            continue
        elif not isinstance(canonical_key, string_types):
            for k in canonical_key:
                normalized[k] += energy.in_units_of(units.kilojoules_per_mole)
        else:
            normalized[canonical_key] += energy.in_units_of(units.kilojoules_per_mole)

    if 'Non-bonded' in canonical_keys:
        normalized['nonbonded'] = energy_dict['Non-bonded']
    else:
        normalized['nonbonded'] = normalized['vdw total'] + normalized['coulomb total'] + normalized['h-bond']

    # could be a problem here since desmond calls UB angle as stretch = bond
    normalized['bonded'] = (normalized['bond'] + normalized['angle'] + normalized['dihedral'])

    return normalized


def record_exception(logger, e_out, e_outfile, e):
    logger.exception(e)
    e_out.append(-1)
    e_outfile.append(-1)


def parse_args(args):
    parser = argparse.ArgumentParser(description='Perform a file conversion')

    # Input arguments.
    group_in = parser.add_argument_group('Choose input conversion format')
    group_in.add_argument('--des_in', metavar='file',
            help='.cms file for conversion from DESMOND file format')
    group_in.add_argument('--gro_in', nargs=2, metavar='file',
            help='.gro and .top file for conversion from GROMACS file format')
    group_in.add_argument('--lmp_in', metavar='file',
            help='input file for conversion from LAMMPS file format (expects'
                 ' data file in same directory and a read_data call)')
    group_in.add_argument('--amb_in', nargs=2, metavar='file',
            help='input files (.prmtop) and (.crd or equivalent) for conversion from AMBER file format')
    group_in.add_argument('--crm_in', metavar='file',
            help='input file for conversion from CHARMM file format')

    # Output arguments.
    group_out = parser.add_argument_group('Choose output conversion format(s)')
    group_out.add_argument('--desmond', action='store_true',
            help='convert to DESMOND')
    group_out.add_argument('--gromacs', action='store_true',
            help='convert to GROMACS')
    group_out.add_argument('--lammps', action='store_true',
            help='convert to LAMMPS')
    group_out.add_argument('--amber', action='store_true',
            help='convert to AMBER')
    group_out.add_argument('--charmm', action='store_true',
            help='convert to CHARMM')

    # Other options.
    group_misc = parser.add_argument_group('Other optional arguments')
    group_misc.add_argument('--odir', metavar='directory', default=os.getcwd(),
            help='specification of output directory (default: ./)')
    group_misc.add_argument('--oname', metavar='prefix', default='',
            help='specification of prefix for output filenames '
                 '(default: inputprefix_converted)')

    group_misc.add_argument('-e', '--energy', dest='energy', action='store_true',
            help='evaluate energy of input and output files for comparison')
    group_misc.add_argument('--inefile', dest='inefile', default='',
            help='optional run settings file for input energy evaluation (e.g. .cfg, .mdp, .input)')

    # desmond settings 
    group_misc.add_argument('-dp', '--despath', dest='desmond_path',
            metavar='path', default='',
            help='path for DESMOND binary, needed for energy evaluation')
    group_misc.add_argument('-ds', '--desmondsettings', dest='desmond_set',
            metavar='settings', default=None,
            help='Desmond .cfg settings file used for energy evaluation')

    # gromacs settings
    group_misc.add_argument('-gp', '--gropath', dest='gromacs_path',
            metavar='path', default='',
            help='path for GROMACS binary, needed for energy evaluation')
    group_misc.add_argument('-gs', '--gromacssettings', dest='gromacs_set',
            metavar='settings', default=None,
            help='Gromacs .mdp settings file used for energy evaluation')

    # lammps settings
    group_misc.add_argument('-lp', '--lmppath', dest='lammps_path',
            metavar='path', default='',
            help='path for LAMMPS binary, needed for energy evaluation')
    group_misc.add_argument('-ls', '--lammpssettings', dest='lmp_settings',
            #metavar='settings', default="pair_style lj/cut/coul/long 15.0 15.0\npair_modify tail yes\nkspace_style pppm 1e-8\n\n",
            metavar='settings', default='pair_style lj/cut/coul/long 9.0 9.0\npair_modify tail yes\nkspace_style pppm 1e-8\n\n',
            help='pair_style string to use in the output file. Default is a periodic Ewald simulation')

    # amber settings
    group_misc.add_argument('-ap', '--amberpath', dest='amber_path',
            metavar='path', default='',
            help='path for AMBER binary, needed for energy evaluation')
    group_misc.add_argument('-as', '--ambersettings', dest='amber_set',
            metavar='settings', default=None,
            help='AMBER .in file used for energy evaluation')


    # charmm settings
    group_misc.add_argument('-cp', '--charmmpath', dest='charmm_path',
            metavar='path', default='',
            help='path for CHARMM binary, needed for energy evaluation')
    group_misc.add_argument('-cs', '--charmmsettings', dest='charmm_settings',
            metavar='settings', default='nbond inbfrq -1 imgfrq -1 -\nelec ewald pmew fftx 48 ffty 48 fftz 48 kappa 0.20822755 order 4 -\nvdw vips cutnb 15. cutim 15. ctofnb 15. ctonnb 15.',
            help='CHARMM .in settings used for energy evaluation. Because of the unified input file, settings lines rather than a file are part')


    group_misc.add_argument('-f', '--force', dest='force', action='store_true',
            help='ignore warnings (NOTE: may lead to partially correct topologies)')
    group_misc.add_argument('-v', '--verbose', dest='verbose', action='store_true',
            help='high verbosity, includes DEBUG level output')
    group_misc.add_argument('-n', '--noncanonical', dest='noncanonical', action='store_true',
            help='print the \'noncanonical\' enery terms, the ones that only occur in 1 or 2 programs')

    # Prints help if no arguments given.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(args)


def main(args=None):
    logger.info('Beginning InterMol conversion\n')

    logger.info('InterMol is research software. If you make use of InterMol '
                'in scientific publications please cite the following '
                'reference:\n'
                'Shirts, M.R., Klein, C., Swails, J.M. et al. J Comput Aided Mol Des (2016). doi:10.1007/s10822-016-9977-1\n')

    if not args:
        args = vars(parse_args(args))

    if args.get('gromacs_path'):
        gmx.GMX_PATH = args['gromacs_path']
    if args.get('lammps_path'):
        lmp.LMP_PATH = args['lammps_path']
    if args.get('desmond_path'):
        des.DES_PATH = args['desmond_path']
    if args.get('amber_path'):
        amb.AMB_PATH = args['amber_path']
    if args.get('charmm_path'):
        crm.CRM_PATH = args['charmm_path']

    if args.get('verbose'):
        h.setLevel(logging.DEBUG)

    # Print warnings.
    warnings.simplefilter("always")
    if not args.get('force'):
        # Warnings will be treated as exceptions unless force flag is used.
        warnings.simplefilter("error")

    # --------------- PROCESS INPUTS ----------------- #
    if args.get('gro_in'):
        system, prefix, gro_in, top_in = _load_gromacs(gromacs_files=args['gro_in'])
    elif args.get('des_in'):
        cms_in = args['des_in']
        system, prefix = _load_desmond(cms_in=cms_in)
    elif args.get('lmp_in'):
        lammps_in = args['lmp_in']
        system, prefix = _load_lammps(lammps_in=lammps_in)
    elif args.get('amb_in'):
        system, prefix, prmtop_in, crd_in, amb_structure = _load_amber(amber_files=args['amb_in'])
    elif args.get('crm_in'):
        charmm_in = args['crm_in']
        system, prefix = _load_charmm(charmm_in=charmm_in)
    else:
        logger.error('No input file')
        sys.exit(1)

    # --------------- WRITE OUTPUTS ----------------- #
    if not args.get('oname'):
        oname = '{0}_converted'.format(prefix)
    else:
        oname = args['oname']
    oname = os.path.abspath(os.path.join(args['odir'], oname))

    output_status = dict()
    if args.get('gromacs'):
        _save_gromacs(system, oname, output_status)
    if args.get('lammps'):
        _save_lammps(system, oname, output_status, args)
    if args.get('desmond'):
        _save_desmond(system, oname, output_status)
    if args.get('amber'):
        _save_amber(system, oname, output_status)
    if args.get('charmm'):
        # currently, this only works if amb_in is used. Reason is that
        # charmm does not support RB dihedrals.
        if args.get('amb_in'):
            prms, rtfs = _save_charmm(amb_structure, oname, output_status)
        else:
            errstr = "Can't convert to CHARMM unless inputs are in AMBER"
            output_status['charmm'] = errstr
            logger.error(errstr)

    # --------------- ENERGY EVALUATION ----------------- #
    if args.get('energy'):
        # Run control file paths.
        tests_path = os.path.abspath(os.path.dirname(intermol.tests.__file__))

        # Default locations of setting control files.
        mdp_in_default = os.path.abspath(os.path.join(tests_path, 'gromacs', 'grompp.mdp'))
        cfg_in_default = os.path.abspath(os.path.join(tests_path, 'desmond', 'onepoint.cfg'))
        in_in_default = os.path.abspath(os.path.join(tests_path, 'amber', 'min.in'))

        # this hardcoding of which input energy files to use should not be here; it should be in the testing files.
        if args.get('gro_in'):
            if args.get('inefile'):
                if os.path.splitext(args.get('inefile'))[-1] != '.mdp':
                    logger.warning("GROMACS energy settings file does not end with .mdp")
                mdp_in = args['inefile']
            else:
                mdp_in = mdp_in_default
            input_type = 'gromacs'
            e_in, e_infile = gmx.energies(top_in, gro_in, mdp_in, gmx.GMX_PATH)
            e_in = canonicalize_energy_names(e_in, gmx.to_canonical)

        elif args.get('lmp_in'):
            if args.get('inefile'):
                logger.warning("LAMMPS energy settings should not require a separate infile")
                input_type = 'lammps'
            e_in, e_infile = lmp.energies(lammps_in, lmp.LMP_PATH)
            e_in = canonicalize_energy_names(e_in, lmp.to_canonical)

        elif args.get('des_in'):
            if args.get('inefile'):
                if os.path.splitext(args.get('inefile'))[-1] != '.cfg':
                    logger.warning("DESMOND energy settings file does not end with .cfg")
                cfg_in = args['inefile']
            else:
                cfg_in = cfg_in_default
            input_type = 'desmond'
            e_in, e_infile = des.energies(cms_in, cfg_in, des.DES_PATH)
            e_in = canonicalize_energy_names(e_in, des.to_canonical)

        elif args.get('amb_in'):
            if args.get('inefile'):
                if os.path.splitext(args.get('inefile'))[-1] != '.in':
                    logger.warning("AMBER energy settings file does not end with .in")
                in_in = args['inefile']
            else:
                in_in = in_in_default
            input_type = 'amber'
            e_in, e_infile = amb.energies(prmtop_in, crd_in, in_in, amb.AMB_PATH)
            e_in = canonicalize_energy_names(e_in, amb.to_canonical)

        elif args.get('crm_in'):
            if args.get('inefile'):
                logger.warning("Original CHARMM input file is being used, not the supplied input file")
            input_type = 'charmm'
            # returns energy file
            e_in, e_infile = crm.energies(args.get('crm_in'), crm.CRM_PATH)
            e_in = canonicalize_energy_names(e_in, crm.to_canonical)
        else:
            logger.warning('No input files identified! Code should have never made it here!')

        # Evaluate output energies.
        output_type = []
        e_outfile = []
        e_out = []

        if args.get('gromacs') and output_status['gromacs'] == 'Converted':
            output_type.append('gromacs')
            if args.get('gromacs_set'):
                mdp = args.get('gromacs_set')
            else:
                mdp = mdp_in_default
            try:
                out, outfile = gmx.energies('{0}.top'.format(oname),
                                            '{0}.gro'.format(oname),
                                            mdp, gmx.GMX_PATH)
                out = canonicalize_energy_names(out, gmx.to_canonical)
            except Exception as e:
                record_exception(logger, e_out, e_outfile, e)
                output_status['gromacs'] = e
            else:
                output_status['gromacs'] = potential_energy_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)

        if args.get('lammps') and output_status['lammps'] == 'Converted':
            output_type.append('lammps')
            try:
                out, outfile = lmp.energies('{0}.input'.format(oname),
                                            lmp.LMP_PATH)
                out = canonicalize_energy_names(out, lmp.to_canonical)
            except Exception as e:
                record_exception(logger, e_out, e_outfile, e)
                output_status['lammps'] = e
            else:
                output_status['lammps'] = potential_energy_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)

        if args.get('desmond') and output_status['desmond'] == 'Converted':
            output_type.append('desmond')
            if args.get('desmond_set'):
                cfg = args.get('desmond_set')
            else:
                cfg = cfg_in_default
            try:
                out, outfile = des.energies('{0}.cms'.format(oname),
                                            cfg, des.DES_PATH)
                out = canonicalize_energy_names(out, des.to_canonical)
            except Exception as e:
                record_exception(logger, e_out, e_outfile, e)
                output_status['desmond'] = e
            else:
                output_status['desmond'] = potential_energy_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)

        if args.get('amber') and output_status['amber'] == 'Converted':
            output_type.append('amber')
            amber_in = args.get('amber_set') or in_in_default

            try:
                out, outfile = amb.energies('{}.prmtop'.format(oname),
                                            '{}.rst7'.format(oname),
                                            amber_in,
                                            amb.AMB_PATH)
                out = canonicalize_energy_names(out, amb.to_canonical)
            except Exception as e:
                record_exception(logger, e_out, e_outfile, e)
                output_status['amber'] = e
            else:
                output_status['amber'] = potential_energy_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)

        if args.get('charmm') and output_status['charmm'] == 'Converted':
            output_type.append('charmm')
            if args.get('inefile'):
                logger.warning("CHARMM energy input file not used, information recreated from command line options")
            inpfile = os.path.join(oname, '{0}.inp'.format(oname))
            crm.write_input_file(inpfile,
                                 oname + '.psf',
                                 rtfs,
                                 prms,
                                 [],
                                 crm.pick_crystal_type(amb_structure.box),
                                 amb_structure.box,
                                 oname + '.crd',
                                 args.get('charmm_settings'), ignore_warnings=True)
            try:
                out, outfile = crm.energies(inpfile, crm.CRM_PATH)
                out = canonicalize_energy_names(out, crm.to_canonical)
            except Exception as e:
                record_exception(logger, e_out, e_outfile, e)
                output_status['charmm'] = e
            else:
                output_status['amber'] = potential_energy_diff(e_in, out)
                e_out.append(out)
                e_outfile.append(outfile)


        # Display energy comparison results.
        out = ['InterMol Conversion Energy Comparison Results', '',
               '{0} input energy file: {1}'.format(input_type, e_infile)]
        for out_type, file in zip(output_type, e_outfile):
            out.append('{0} output energy file: {1}'.format(out_type, file))
        out += summarize_energy_results(e_in, e_out, input_type, output_type, args.get('noncanonical'))
        logger.info('\n'.join(out))

    logger.info('Finished!')
    return output_status


def potential_energy_diff(e_in, e_out):
    """Returns difference in potential energy.

    arguments:
        e_in  - dictionary of energy groups from input file
        e_out - dictionary of energy groups from output file

    returns:
        potential energy difference in units of the input
    """
    energy_type = 'potential'
    input_energy = e_in[energy_type]
    diff = e_out[energy_type].in_units_of(input_energy.unit) - input_energy
    return diff._value


def find_match(key, dictionary, unit):
    """Helper function for `summarize_energy_results`. """
    if key in dictionary:
        return dictionary[key].value_in_unit(unit)
    else:
        return np.nan


def summarize_energy_results(energy_input, energy_outputs, input_type, output_types, print_noncanonical):
    """Creates a table comparing input and output energy groups.

    Args:
        energy_input (dict): energy groups from input file
        energy_output(list): containing dictionary of energy groups or -1 for
            each output file
        input_type (str): input engine
        output_types (list): containing output formats
        print_noncanonical (bool): if true, print noncanonical energy terms

    Returns:
        out (list of strings): which forms a summary table using "\n".join(out)
    """
    out = []
    # Remove failed evaluations (-1 in energy_outputs)
    failed_i = [i for i, x in enumerate(energy_outputs) if x == -1]
    failed = [output_types[i] for i in failed_i]
    output_types = [x for i, x in enumerate(output_types) if i not in failed_i]
    energy_outputs = [x for x in energy_outputs if x != -1]

    # Find all distinct labels
    labels = set(energy_input.keys())
    for e_out in energy_outputs:
        for key, value in e_out.items():
            labels.add(key)

    # Set up energy comparison table
    labels = sorted(list(labels))
    unit = energy_input[list(energy_input.keys())[0]].unit
    energy_all = [energy_input] + energy_outputs
    data = np.empty((len(labels), len(energy_all)))
    for i in range(len(data)):
        for j in range(len(energy_all)):
            data[i, j] = find_match(labels[i], energy_all[j], unit)

    # TODO: sort table
    out.append('')
    out.append('Energy group summary')
    out.append('=======================================================================')
    header = '%20s %18s' % ('type', 'input(%s)' % input_type)
    for otype in output_types:
        header += '%18s' % ('output (%s)' % (otype)) 
        header += '%18s' % ('diff (%s)' % (otype))  
    out.append(header)
    types = ['notcomparable','comparable']
    for t in types:
        if t == 'notcomparable':
            out.append('-----------------------------------------------------------------------')
            out.append('Not comparable energies: are likely not to be the same')
            out.append('-----------------------------------------------------------------------')
            # warning on urey_bradley being put in two different categories here?
        elif t == 'comparable':
            out.append('-----------------------------------------------------------------------')
            out.append('Comparable energy terms: these should be very close')
            out.append('-----------------------------------------------------------------------')
        for i in range(len(data)):
            if t == 'notcomparable':
                if labels[i] in comparable_energy_terms:
                    continue
            elif t == 'comparable':
                if labels[i] not in comparable_energy_terms:
                    continue
            are_all_zero = True   # if all terms are zero, don't bother printing
            for j in range(0, len(data[i])):
                if data[i][j] != 0:
                    are_all_zero = False
            if labels[i] not in canonical_energy_names and not print_noncanonical or are_all_zero:
                continue
            line = '%20s ' % labels[i]
            if np.isnan(data[i][0]):
                line += '%18s' % 'n/a'
            else:
                line += '%18.8f' % (data[i][0])
            for j in range(1, len(data[i])):
                if np.isnan(data[i][j]):
                    line += '%18s' % 'n/a'
                else:
                    line += '%18.8f' % (data[i][j])
                if np.isnan(data[i][j]) or np.isnan(data[i][0]):
                    line += '%18s' % 'n/a'
                else:
                    line += '%18.8f' % (data[i][j]-data[i][0])
            out.append(line)
    out.append('')
    # get differences in potential energy
    i = labels.index('potential')
    diff = data[i, 1::] - data[i, 0]
    out.append('---------------- Total Potential Energy Comparison --------------------')
    out.append('Input %s potential energy: %22.8f' % (input_type,data[i,0])) 
    # print the canonical energies
    for d, otype in zip(diff, output_types):
        out.append('Difference in potential energy from %s=>%s conversion: %18.8f'
                   % (input_type, otype, d))
    for fail in failed:
        out.append('Energy comparison for {0} output failed'.format(fail))
    out.append('=======================================================================')
    return out


def _load_gromacs(gromacs_files):
    prefix = os.path.splitext(os.path.basename(gromacs_files[0]))[0]
    # Find the top file since order of inputs is not enforced.
    top_in = [x for x in gromacs_files if x.endswith('.top')]
    assert(len(top_in) == 1)
    top_in = os.path.abspath(top_in[0])

    # Find the gro file since order of inputs is not enforced.
    gro_in = [x for x in gromacs_files if x.endswith('.gro')]
    assert(len(gro_in) == 1)
    gro_in = os.path.abspath(gro_in[0])
    system = gmx.load(top_in, gro_in)
    return system, prefix, gro_in, top_in


def _save_gromacs(system, oname, output_status):
    try:
        gmx.save('{0}.top'.format(oname), '{0}.gro'.format(oname), system)
    except Exception as e:
        logger.exception(e)
        output_status['gromacs'] = e
    else:
        output_status['gromacs'] = 'Converted'


def _load_desmond(cms_in):
    prefix = os.path.splitext(os.path.basename(cms_in))[0]
    system = des.load(cms_file=cms_in)
    return system, prefix


def _save_desmond(system, oname, output_status):
    try:
        des.save('{0}.cms'.format(oname), system)
    except Exception as e:
        logger.exception(e)
        output_status['desmond'] = e
    else:
        output_status['desmond'] = 'Converted'


def _load_lammps(lammps_in):
    prefix = os.path.splitext(os.path.basename(lammps_in))[0]
    system = lmp.load(in_file=lammps_in)
    return system, prefix


def _save_lammps(system, oname, output_status, args):
    try:
        lmp.save('{0}.input'.format(oname), system, nonbonded_style=args.get('lmp_settings'))
    except Exception as e:
        logger.exception(e)
        output_status['lammps'] = e
    else:
        output_status['lammps'] = 'Converted'


def _load_amber(amber_files):
    prefix = os.path.splitext(os.path.basename(amber_files[0]))[0]

    # Find the prmtop file since order of inputs is not enforced.
    prmtop_in = [x for x in amber_files if (x.endswith('.prmtop') or
                                            x.endswith('.parm7'))]
    assert len(prmtop_in) == 1
    prmtop_in = os.path.abspath(prmtop_in[0])

    # Find the crd file since order of inputs is not enforced, not is suffix
    crd_in = [x for x in amber_files if (x.endswith('.rst7') or
                                         x.endswith('.crd') or
                                         x.endswith('.rst') or
                                         x.endswith('.inpcrd'))]
    assert len(crd_in) == 1
    crd_in = os.path.abspath(crd_in[0])

    print(prmtop_in, crd_in)

    # some coordinate files don't actually load in the box.  We need to assume something.
    # if there is no box, we assume that it's vacuum, and we just need to make sure it's big.
    # revisit this at some point.
    amb_structure = pmd.load_file(prmtop_in, xyz=crd_in)
    if amb_structure.box is None:
        amb_structure = pmd.load_file(prmtop_in, xyz=crd_in, box = [70,70,70,90,90,90])
        
    # write out the files.  Should write them out in the proper directory (the one reading in)
    pathprefix = os.path.dirname(prmtop_in)
    fromamber_top_in = os.path.join(pathprefix, prefix + '_from_amber.top')
    fromamber_gro_in = os.path.join(pathprefix, prefix + '_from_amber.gro')
    # pmd.gromacs.GromacsTopologyFile.write(parmed_system, fromamber_top_in)
    # pmd.gromacs.GromacsGroFile.write(parmed_system, fromamber_gro_in, precision = 8)
    amb_structure.save(fromamber_top_in, overwrite=True)
    amb_structure.save(fromamber_gro_in, precision=8, overwrite=True)

    # now, read in using gromacs
    system = gmx.load(fromamber_top_in, fromamber_gro_in)
    return system, prefix, prmtop_in, crd_in, amb_structure


def _save_amber(system, oname, output_status):
    # To get to AMBER, we write to GROMACS and then convert with ParmEd
    _save_gromacs(system, oname, output_status)
    try:
        top = pmd.load_file(oname + '.top', xyz=oname + '.gro')
        top.save(oname + '.prmtop', overwrite=True)
        top.save(oname + '.rst7', overwrite=True)
    except Exception as e:
        logger.exception(e)
        output_status['amber'] = e
    else:
        output_status['amber'] = 'Converted'


def _load_charmm(charmm_in):
    prefix = os.path.splitext(os.path.basename(charmm_in))[0]
    # we need to find the parameter and structure files by reading the input file.
    with open(charmm_in) as cinf:
        lines = cinf.readlines()
    box = []
    rtfs = []
    prms = []
    strms = []
    topsuffixes = ['.rtf', '.top']
    prmsuffixes = ['.prm', '.par']

    for line in lines:
        # Find the psf file
        if '.psf' in line:
            psffile = line.split()[-1] # append the file at the end of the line
        if '.crd' in line:
            crdfile = line.split()[-1] # append the file at the end of the line
        if '.str' in line:
            strms.append(line.split()[-1]) # append the file at the end of the line
        if any(x in line for x in topsuffixes): # if the file is any of the topology suffixes: need to be read in before prms.
            rtfs.append(line.split()[-1])  # append the file at the end of the line
        if any(x in line for x in prmsuffixes): # if the file is any of the parameter suffixes
            prms.append(line.split()[-1])  # append the file at the end of the line
        if 'set box ' in line:   # will need to handle general variables
            boxlength_vars = line.split()
            boxval = np.float(boxlength_vars[2])
        if 'crystal define' in line:
            boxangle_vars = line.split()
            boxvecs = boxangle_vars[3:6]
            for i, b in enumerate(boxvecs):
                if '@box' in b:
                    boxvecs[i] = boxval
                else:
                    boxvecs[i] = np.float(b)
            boxangles = np.array(boxangle_vars[6:9],float)
            box = np.append(np.array(boxvecs),boxangles)

    #load in the parameters
    psf = pmd.load_file(psffile)
    parameterset = pmd.charmm.CharmmParameterSet()
    for tfile in rtfs:
        parameterset.read_topology_file(tfile)
    for pfile in prms:
        parameterset.read_parameter_file(pfile)
    for sfile in strms:
        parameterset.read_stream_file(sfile)
    psf.load_parameters(parameterset)

    if len(box) == 6:
        psf.box = box
    # now load in the coordinates
    crd = pmd.load_file(crdfile)
    try:
        if len(crd.coordinates.shape) == 3:
            coords = crd.coordinates[0]
        else:
            coords = crd.coordinates
    except:
        logger.error('No coordinates in %s' % (crd))

    if coords.shape != (len(psf.atoms), 3):
        logger.error('Mismatch in number of coordinates (%d) and '
            '3*number of atoms (%d)' % (len(coords), 3*len(psf.atoms)))
        # Set the coordinates now, since creating the parm may re-order the
        # atoms in order to maintain contiguous molecules
    psf.coordinates = coords

    # copy the box over to the .psf
    if hasattr(crd, 'box') and crd.box is not None:
        if len(crd.box.shape) == 1:
            crdbox = crd.box
        else:
            # Trajectory
            crdbox = crd.box[0]

        if len(crdbox) == 3:
            psf.box = list(crdbox) + [90.0, 90.0, 90.0]
        elif len(crdbox) == 6:
            psf.box = list(crdbox)
        else:
            logger.error('Unexpected box array shape')

    #Make GROMACS topology
    parmed_system = pmd.gromacs.GromacsTopologyFile.from_structure(psf)

    # write out the files.  Should write them out in the proper directory (the one reading in)
    pathprefix = os.path.dirname(charmm_in)
    fromcharmm_top_in = os.path.join(pathprefix, prefix + '_from_charmm.top')
    fromcharmm_gro_in = os.path.join(pathprefix, prefix + '_from_charmm.gro')
    pmd.gromacs.GromacsTopologyFile.write(parmed_system, fromcharmm_top_in)
    pmd.gromacs.GromacsGroFile.write(parmed_system, fromcharmm_gro_in, precision=8)

    # now, read in using gromacs
    system = gmx.load(fromcharmm_top_in, fromcharmm_gro_in)
    return system, prefix


def _save_charmm(amb_structure, oname, output_status):
    try:
        parmed_system = pmd.charmm.CharmmPsfFile.from_structure(amb_structure)
        charmm_output_psf = '{0}.psf'.format(oname)
        charmm_output_crd = '{0}.crd'.format(oname)

        # We need these arrays for enery output.
        charmm_output_rtf = '{0}.rtf'.format(oname)
        charmm_output_prm = '{0}.prm'.format(oname)
        prms = [charmm_output_prm]
        rtfs = [charmm_output_rtf]
        pmd.charmm.CharmmParameterSet.write(
            pmd.charmm.CharmmParameterSet.from_structure(amb_structure),
            top=charmm_output_rtf,
            par=charmm_output_prm)
        amb_structure.save(charmm_output_psf, format='psf', overwrite=True)
        pmd.charmm.CharmmCrdFile.write(parmed_system, charmm_output_crd)
    except Exception as e:
        logger.exception(e)
        output_status['charmm'] = e
    else:
        output_status['charmm'] = 'Converted'

    return prms, rtfs

if __name__ == '__main__':
    # Specifies lowest severity log messages to handle.
    logger.setLevel(logging.DEBUG)
    h = logging.StreamHandler()
    h.setLevel(logging.INFO)  # Ignores DEBUG level for now.
    f = logging.Formatter("%(levelname)s %(asctime)s %(message)s",
                          "%Y-%m-%d %H:%M:%S")
    h.setFormatter(f)
    logger.addHandler(h)

    # Redirect warnings module messages to logging system.
    logging.captureWarnings(True)
    warning_logger = logging.getLogger('py.warnings')
    warning_logger.addHandler(h)

    main()
