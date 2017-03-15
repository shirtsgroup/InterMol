import parmed.unit as units

import intermol.forces.forcedata as forcedata

"""
functions for manipulating the data to extract keywords and unit-ed parameter lists from forces, 

name in code              example                           description 
force_type     -      HarmonicBondType      - class, contains atom types and parameters
force_class    -        HarmonicBond        - class, XType, also adds atoms
force          - (instance of HarmonicBond) - an instance of HarmonicBond
"""


def specify(program_units, unitset, dumself=None, shouldEval=True):
    """Takes the dict of units, and a set of dimensions and replaces the dimensions with the appropriate units.
    """
    specified_unitset = []
    for unit in unitset:
        specified_unit = []
        for chunk in unit.split():
            if chunk in program_units:
               chunk = program_units[chunk]
            specified_unit.append(chunk)
        rejoined_unit = ' '.join(specified_unit)
        if shouldEval:
            specified_unitset.append(eval(rejoined_unit))
        else:
            specified_unitset.append(rejoined_unit)
    return specified_unitset


def build_paramlist(program):
    """Create a paramlist specific for a given program. """
    change_list = eval('forcedata.' + program + '_paramlist')
    tmp_paramlist = forcedata.master_paramlist.copy()
    tmp_paramlist.update(change_list)

    paramlist = tmp_paramlist.copy()
    # add type and underscore names
    for name, paramset in tmp_paramlist.items():
        paramlist[capifyname(name)] = tmp_paramlist[name]
        paramlist[capifyname(name + '_type')] = tmp_paramlist[name]

    return paramlist


def capifyname(forcename):
    """
    Return name of the class in camelCase.
    """
    return forcename.replace('_',' ').title().replace(' ','')


def build_unitvars(program, paramlist, dumself=None):
    """
    Takes a string program name (one of the supported programs), and a 'self' object
    it looks like the keyword is not being used, but it is used in the line eval(unit). 
    The test name 'dumself' needs to match what is in the force data arrays. Currently only used for lammps.
    """
    unitvars = dict()
    unitdefs = forcedata.ProgramUnitSets[program]
    for name, uset in forcedata.master_unitlist.items():
        unitset = specify(unitdefs, uset, dumself)

        # reorder the units if necessary according to the order in the given paramlist
        original_params = forcedata.master_paramlist[name]
        program_params = paramlist[name]
        tmp_unitset = []
        if original_params != program_params:
            for i, op in enumerate(original_params):
                if op in program_params:
                    tmp_unitset.insert(program_params.index(op),unitset[i])
            unitset = tmp_unitset    

        if name in forcedata.ProgramUnitLists:
            # In case the units need to be defined differently.
            unitset = forcedata.ProgramUnitLists[name]
        unitvars[capifyname(name)] = unitset
        typename = name  + '_type'
        unitvars[typename] = unitset
        unitvars[capifyname(typename)] = unitset
    return unitvars


def get_parameter_list_from_force(force, paramlist):
    """Create a function that returns the paramters of a function type.

    First, we need make some additions to the parameter list dictionary,
    which we do once when the forcedata script is imported.  Useful to
    put the forces here as well.  We won't make this a function for now
    since it's needed in this module.
    """

    # We passed in an instance
    name = force.__class__.__name__
    pvars = []
    for param in paramlist[name]:
        paramstring = 'force.' + param
        pvars.append(eval(paramstring))
    return pvars


def get_parameter_list_from_kwds(force, kwds, paramlist):
    """ """
    # We passed in an instance, not a class
    name = force.__class__.__name__
    ordered = []
    for p in paramlist[name]:
        ordered.append(kwds[p])
    return ordered


def get_parameter_kwds_from_force(force, forceparams, paramlist):
    """ """
    kwds = dict()

    force_params = forceparams(force)
    for i, p in enumerate(paramlist[force.__class__.__name__]):
        kwds[p] = force_params[i]
    return kwds


def create_kwds_from_entries(unitvars, paramlist, entries, force_type, offset=0):
    """Create a keyword dictionary given an array of information from a file format

    requires the master set of units, the master set of parameter
    lists, an object (either a force_class or force_type), the
    list of information to be converted into a keyword, and an offset.

    Args:
        offset (int): how far over from the first entry we translate
    """
    kwds = dict()
    typename = force_type.__name__
    u = unitvars[typename]
    params = paramlist[typename]
    for i, p in enumerate(params):
        kwds[p] = float(entries[offset+i]) * u[i]
    return kwds


def optparamkeylookup(force_type):
    """Given a force_type object, determine the key associated with the
    optional parameters.

    """
    name = force_type.__name__.lower()
    for key, params in forcedata.AbstractOptParams.items():
        if key in name:
            return key


def optforceparams(force_type, forcetype_object=None):
    """Return the dictionary of optional paramters of an abstract force type.

    If no object is given, we fill with blanks.
    """
    pvars = dict()
    #MRS: should be able to get rid of the evals?  Apparently, no unit tests for this code yet, will get rid
    #when they are added.

    for i, param in enumerate(forcedata.AbstractOptParams[force_type]):
        if forcetype_object:
            pvars[param] = eval(forcetype_object.__class__.__name__ + '.' + param)
        else:
            pvars[param] = eval(forcedata.AbstractOptParamsDefaults[force_type][i])
    return pvars


def optparamlookup(force_type_object, object_default=False):
    """A wrapper for optforceparams that takes a force_type object and returns
    the optional parameter dictionary.
    """
    force_type = optparamkeylookup(force_type_object)
    if object_default:
       return optforceparams(force_type, force_type_object)
    else:
       return optforceparams(force_type)


def create_kwd_dict(unitvars, paramlist, force_type_object, values, optvalues=None):
    """ """
    name = force_type_object.__name__
    unitlist = unitvars[name]
    kwdlist =  paramlist[name]
    optkwddict = optparamlookup(force_type_object)
    
    arglist = [unit*value for unit, value in zip(unitlist, values)]
    kwd = {key: value for key, value in zip(kwdlist, arglist)}

    if optvalues:
       optkwddict.update(optvalues)
       kwd.update(optkwddict)
    return kwd
