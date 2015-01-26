from collections import OrderedDict
import numpy as np
import pdb

def get_diff(e_in, e_out):
    '''returns difference in potential energy
    
    arguments:
        e_in  - dictionary of energy groups from input file
        e_out - dictionary of energy groups from output file
    
    returns:
        potential energy difference in units of the input
    '''
    type = 'Potential' # getting difference in potential energy
    input = e_in[type]
    diff = e_out[type].in_units_of(input.unit) - input
    return diff._value

def find_match(key, dict, unit):
    '''helper function for multiple_energy_results()'''
    if key in dict:
        return dict[key].in_units_of(unit)._value
    else:
        return np.nan

def summarize_energy_results(energy_input, energy_outputs, input_type, output_types):
    ''' creates a table comparing input and output energy groups

    args:
        energy_input  - dictionary of energy groups from input file
        energy_output - list containing dictionary of energy groups or -1 for each output file
        input_type    - input format string
        output_types  - list containing output formats
    returns:
        out - list of strings which forms a summary table using "\n".join(out) 
    '''
    out = []
    # remove failed evaluations (-1 in energy_outputs)
    failed_i = [i for i,x in enumerate(energy_outputs) if x == -1]
    failed = [output_types[i] for i in failed_i]
    output_types = [x for i,x in enumerate(output_types) if i not in failed_i]
    energy_outputs = [x for x in energy_outputs if x != -1]
    # find all distinct labels
    labels = set(energy_input.keys())
    for e_out in energy_outputs:
        for key, value in e_out.iteritems():
            labels.add(key)
    # set up energy comparison table
    labels = list(labels)
    unit = energy_input[energy_input.keys()[0]].unit
    energy_all = [energy_input] + energy_outputs
    data = np.empty((len(labels),len(energy_all)))
    for i in range(len(data)):
        for j in range(len(energy_all)):
            data[i,j] = find_match(labels[i], energy_all[j], unit)
    # sort table <= To do
    # print table
    out.append('')
    out.append('Energy Group Summary')
    out.append('=======================================================================')
    header = '%20s %18s ' % ('Type', 'Input (%s)' % input_type)
    for otype in output_types:
        header += '%37s' %('Output (%s) Diff (%s)' % (otype, otype))
    out.append(header)
    for i in range(len(data)):
        line = '%20s ' % labels[i]
        line += '%18.8f ' % data[i][0]
        for j in range(1,len(data[i])):
            line += '%18.8f %18.8f' % (data[i][j],data[i][j]-data[i][0])
        out.append(line)
    out.append('')
    # get differences in potential energy
    i = labels.index('Potential')
    diff = data[i,1::] - data[i,0]
    for d, otype in zip(diff, output_types):
        out.append('Difference in potential energy from %s=>%s conversion: %18.8f' 
                    % (input_type, otype, d))
    for f in failed:
        out.append('Energy comparison for {0} output FAILED'.format(f))
    out.append('=======================================================================')
    return out

def combine_energy_results(e_in, e_out):
    '''deprecated
    '''
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
    '''deprecated
    '''
    data, rms = results
    print '======================================================================='
    print 'Summary statistics'
    print '%20s %18s %18s %18s' % ('Type', 'Input', 'Output', 'Diff')
    for (name, values) in data.iteritems():
        print '%20s %18.8f %18.8f %18.8f' % (name,
                values[0]._value if hasattr(values[0], '_value') else values[0], # in
                values[1]._value if hasattr(values[1], '_value') else values[1], # out
                values[2]._value if hasattr(values[2], '_value') else values[2]) # diff

    print ' '
    print 'RMS signed error: %18.8f' % (rms)
    return rms
