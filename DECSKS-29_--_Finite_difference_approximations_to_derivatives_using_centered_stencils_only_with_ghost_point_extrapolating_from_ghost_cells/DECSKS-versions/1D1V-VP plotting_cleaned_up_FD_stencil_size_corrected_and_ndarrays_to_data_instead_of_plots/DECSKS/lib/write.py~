import numpy as np
import os

def tofiles(sim_name, n, **kwargs):

    time_label = '%06d' % n
    # **kwargs should be provided as key/value pairs where each key is name of a quantity (used in the filename)
    # and each value is the ndarray (data) to be stored to file

    for key in kwargs.keys():

        # key name is the quantity name
        quantity = key

        # its value is the ndarray of data
        data = kwargs[key]

        # create filename, will be stored as numpy binary data file, .npy ( http://docs.scipy.org/doc/numpy/neps/npy-format.html )
        filename = sim_name + '_--_' + quantity + '_' + time_label + '.npy'

        # create filepath, this module will usually be called from the main.py directory, DECSKS/
        rel_path = './etc/outputs/' + sim_name + '/'
        filepath = rel_path + filename

        with open(filepath, 'w') as outfile:
            np.save(outfile, data)
