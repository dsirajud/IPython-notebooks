import numpy as np
import os
import errno

def make_sure_path_exists(path):
    """
    this code is authored by Heikki Toivonen
    as posted at:

    http://stackoverflow.com/questions/
        273192/how-to-check-if-a-directory-exists
        -and-create-it-if-necessary

    there is a race error that plagues naive methods
    that causes the OSError, this catches that

    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def tofiles(sim_name, n = None, **kwargs):

    if n == None:
        time_label = ''
    else:
        time_label = '_%06d' % n

    # directory where outputs will be stored: /etc/outputs/sim_name
    dir_path = './etc/outputs/' + sim_name + '/'
    make_sure_path_exists(dir_path)

    # **kwargs should be provided as key/value pairs where each key is name of a quantity (used in the filename)
    # and each value is the ndarray (data) to be stored to file
    for key in kwargs.keys():

        # key name is the quantity name
        quantity = key

        # its value is the ndarray of data
        data = kwargs[key]

        # create filename, will be stored as numpy binary data file, .npy ( http://docs.scipy.org/doc/numpy/neps/npy-format.html )
        filename = quantity + time_label + '.npy'

        # create filepath, this module will usually be called from the main.py directory, DECSKS/
        filepath = dir_path + filename

        with open(filepath, 'w') as outfile:
            np.save(outfile, data)

    return None
