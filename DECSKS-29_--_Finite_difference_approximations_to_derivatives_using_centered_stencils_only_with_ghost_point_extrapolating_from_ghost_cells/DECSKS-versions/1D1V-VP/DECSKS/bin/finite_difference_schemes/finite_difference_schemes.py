import linecache

def store(dn_max = 4):
    """store all finite difference schemes from
    tables generated in
    generate_table_of_finite_difference_schemes.py
    in a dictionary FD_schemes

    inputs:
    dn_max -- (int) maximum derivative in .dat files
              should correspond to same dn_max as
              in tables generated

    outputs:
    FD_schemes -- (dict) all FD schemes equipped with
                  list of weights w and stencil
    """
    FD_schemes = {}
    infile_suffix = '_FD_coefficients.dat'

    for dn in range(1,dn_max+1):
        infilename = 'f' + str(dn) + infile_suffix

        # create empty subdictionary for given dn
        FD_schemes['dn' + str(dn)] = {}

        # store schemes in subdictionaries of dn
        # which are the keys: 'handedness', 'asymmetry'
        FD_schemes = read(
                                 infilename,
                                 FD_schemes,
                                 str(dn))

    return FD_schemes

def read(infilename,
         FD_schemes,
         dn):
    """reads infile, creates empty subdictionaries
    inside FD_schemes, and stores all data inside
    the corresponding dictionary objects

    inputs:
    infilename -- (str) file name for a single derivative
                  table, e.g. f1_FD_coefficients.dat,
                              f2_FD_coefficients.dat, ...

    FD_schemes -- (dict) empty dictionary
    dn -- (str) single derivative whose schemes are to be stored

    outputs:
    FD_schemes -- (dict) same dictionary loaded with
                  all schemes for the dn'th derivative
    """
    infile = open(infilename, 'r')
    lines = infile.readlines()

    dn = 'dn' + dn # dn key (str)

    # create empty subdictionaries for each handedness
    FD_schemes[dn]['forward'] = {}
    FD_schemes[dn]['central'] = {}
    FD_schemes[dn]['backward'] = {}

    for i in range(len(lines)):

        if lines[i][0] == 'F':
            handedness = 'forward'
            asymmetry = lines[i][1]

            # create empty subdictionary for given asymmetry
            FD_schemes[dn][handedness][asymmetry] = {}
        elif lines[i][0] == 'C':
            handedness = 'central'
            asymmetry = lines[i][1]

            # create empty subdictionary for given asymmetry
            FD_schemes[dn][handedness][asymmetry] = {}
        elif lines[i][0] == 'B':
            handedness = 'backward'
            asymmetry = lines[i][1]

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
            # let i increment until it reaches i + numlines and repeat
            pass

    return FD_schemes
