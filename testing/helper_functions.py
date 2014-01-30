import numpy as np

def combine_energy_results(e_in, e_out):
    """
    """
    data = dict()
    for e_type, value in e_in.iteritems():
        if e_type in e_out:
            data[e_type] = (value, 
                    e_out[e_type], 
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

    print "%20s %15s %15s %15s" % ("Type", "Input", "Output", "Diff")
    for (name, values) in data.iteritems():
        print "%20s %15.8f %15.8f %15.8f" % (name, values[0]._value, values[1]._value, values[2]._value)

    print " "
    print "RMS signed error: %10.5f" % (rms)
