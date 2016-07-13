import numpy as np
import pylab
import os

##### DEPRECRATED

# NOTE THIS MODULE IS TOO GENERAL WHICH MAKES IT BOTH CUMBERSOME AND REQUIRING
# FAR TOO MUCH INFORMATION TO BE USED IN SPECIFIC CASES

# INSTEAD WE WILL COMPUTE ANY QUANTITIES OR PERFORM PLOTTING
# BY WRITING SPECIFIC SCRIPTS, E.G. PLOT ION DENSITY, COMPUTE TOTAL ENERGY,  ...



# cd = DECSKS/bin/postprocessing
# outputs are stored in DECSKS/etc/outputs/sim_number
rel_path = '../../etc/outputs/'

def HighPrecisionE(number):
    """Converts a number into a string object
    while retaining a chosen degree of precision. This
    is designed to evade the truncation that is involved
    with str() so that outputs can store numbers with high
    precision

    inputs:
    number -- (number)

    outputs:
    string object with chosen precision in scientific notation
    """

    return "%.22e" % number

def readfile(filename):

    filepath = rel_path + filename
    with open(filename, 'r') as infile:
        data = np.load(filename)

    return data

def readfiles(n_start, n_stop, quantity = None, ndim = 1, dims = None, sim_name = None):

    # should read in specs like dimensions from affiliated "about"" files
    # for each quantity instead of having so many arguments

    assert (type(n_stop) == int)
    assert (type(n_start) == int)
    assert (type(quantity) == str)
    assert (type(ndim) == int)
    assert (type(dims) == tuple)

    print os.getcwd()
    dir_path = rel_path + sim_name + '/'
    filename = sim_name + '_--_' + quantity + '_'
    fileformat = '.npy'

    # if n_start != 0, then we would like to return an array that still
    # has it so the first index corresponds to time step, so we pad with as many zero
    # entries as needed to get to n_start where data is stored in the output array
    # this amounts to increasing the size of the first dimension in the ndarray data
    # by exactly n_start


    total_n = n_stop - n_start + 1
    output_array_dims = [n_start + total_n] # first dimension enumerates time steps
    for dim in dims:
        output_array_dims.append(dim)

    data = np.zeros(output_array_dims)
    print "reading in files for quantity %s for times n = %d to n = %d" % (quantity, n_start, n_stop)
    print "\n..."
    for n in range(n_start, n_stop + 1):
        # generate specific filepath for
        time_label = '%06d' % n
        filepath = dir_path + filename + time_label + fileformat
        print data.shape
        data[n] = readfile(filepath) # stores in first dim, if ndim = 2, this is equivalent to data[n,:,:] for example
        print "read-in and stored the file %s" % (filename + time_label + fileformat)
    print "\nall files read in and returned"

    return data
