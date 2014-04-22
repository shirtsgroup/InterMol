import numpy as np
import pdb

def combine_energy_results(e_in, e_out):
    """
    """
    data = dict()
    for e_type, value in e_in.iteritems():
        if e_type in e_out:
            data[e_type] = (value,
                    e_out[e_type].in_units_of(value.unit),
                    value - e_out[e_type])

    diff = list()
    for value in data.values():
        diff.append(value[2]._value)
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
                values[0]._value,  # in
                values[1]._value,  # out
                values[2]._value)  # diff

    print " "
    print "RMS signed error: %18.8f" % (rms)
    return rms
