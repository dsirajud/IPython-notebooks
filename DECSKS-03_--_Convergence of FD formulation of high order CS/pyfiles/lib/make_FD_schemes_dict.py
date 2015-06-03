import linecache # to specify pointer of line reader
import numpy as np
import numpy.linalg as LA
 

def store(dn, LTE):

    FD_schemes = {}
    infile_suffix = '_FD_coefficients.dat'

    infilename = 'f' + str(dn) + '_' + 'LTE_' + str(LTE) + infile_suffix
    infilepath = './etc/' + infilename

    # create empty subdictionary for given dn
    FD_schemes['dn' + str(dn)] = {}

    # store schemes in subdictionaries of dn
    # which are themselves keys: 'handedness', 'asymmetry'
    FD_schemes = read(
                      infilepath,
                      FD_schemes,
                      str(dn))

    return FD_schemes

def read(infilename,
                FD_schemes,
                dn):

    infile = open(infilename, 'r')
    lines = infile.readlines()

    dn = 'dn' + dn # dn key

    # create empty subdictionaries for each handedness
    FD_schemes[dn]['forward'] = {}
    FD_schemes[dn]['central'] = {}
    FD_schemes[dn]['backward'] = {}
    for i in range(len(lines)):

        if lines[i][0] == 'F':
            handedness = 'forward'
            asymmetry = lines[i][1:].strip()
            # create empty subdictionary for given asymmetry
            FD_schemes[dn][handedness][asymmetry] = {}
        elif lines[i][0] == 'C':
            handedness = 'central'
            asymmetry = lines[i][1] # always '0' if handedness = 'central'
            # create empty subdictionary for given asymmetry
            FD_schemes[dn][handedness][asymmetry] = {}
        elif lines[i][0] == 'B':
            handedness = 'backward'
            asymmetry = lines[i][1:].strip()
            # create empty subdictionary for given asymmetry
            FD_schemes[dn][handedness][asymmetry] = {}
        elif lines[i].split(' ')[0] == 'number':
            numlines = eval(lines[i][lines[i].find(':')+1:].strip())

            # set switch on for storage and intialize empties
            store_data_switch = 1
            w = []
            stencil = []

        elif store_data_switch == 1:
            for j in range(i+1, i + numlines+1):
                # note: first line in linecache is [1], not [0]
                line = linecache.getline(infilename,j)
                pairs = line.split(', ')

                w.append(eval(pairs[0]))
                stencil.append(eval(pairs[1]))

            # store as values in subdictionary keys 'w', 'stencil'
            FD_schemes[dn][handedness][asymmetry]['w'] = w
            FD_schemes[dn][handedness][asymmetry]['stencil'] = stencil

            # reset switch for next scheme dict key
            store_data_switch = 0

        else:
            pass

    return FD_schemes
