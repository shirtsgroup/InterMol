from collections import OrderedDict
import numpy as np
import pdb

def combine_energy_results(e_in, e_out):
    """
    """
    data = OrderedDict()
    matches = []
    diff = []
    for e_type, value in e_in.iteritems():
        if e_type in e_out:
            matches.append(e_type)
            data[e_type] = (value,
                    e_out[e_type].in_units_of(value.unit),
                    value - e_out[e_type].in_units_of(value.unit))
            diff.append(data[e_type][2]._value)


    for e_type, value in e_in.iteritems():
        if e_type not in matches:
            data[e_type] = (value,
                    np.nan,
                    np.nan)

    e_in_unit = value.unit # assuming all e_in terms are in the same unit...
    for e_type, value in e_out.iteritems():
        if e_type not in matches:
            matches.append(e_type)
            data[e_type] = (np.nan,
                    value.in_units_of(e_in_unit),
                    np.nan)

    diff = np.asarray(diff)
    rms = np.sqrt(np.mean(diff ** 2))
    return data, rms

def print_energy_summary(results):
    """
    """
    data, rms = results
    print "======================================================================="
    print "Summary statistics"
    print "%20s %18s %18s %18s" % ("Type", "Input", "Output", "Diff")
    for (name, values) in data.iteritems():
        print "%20s %18.8f %18.8f %18.8f" % (name,
                values[0]._value if hasattr(values[0], '_value') else values[0], # in
                values[1]._value if hasattr(values[1], '_value') else values[1], # out
                values[2]._value if hasattr(values[2], '_value') else values[2]) # diff

    print " "
    print "RMS signed error: %18.8f" % (rms)
    return rms
