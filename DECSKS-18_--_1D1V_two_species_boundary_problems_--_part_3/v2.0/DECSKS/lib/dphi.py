import numpy as np
import linecache

#-----------------------------------------------------------#
# the following are for the sole purpose of reading a single
# dn at given LTE, the original purpose was to have dn1, LTE6
# so that dphi can be computed

def read_FD_scheme(dn, LTE):
    """store finite difference scheme for dn'th derivative
    from tables generated in dat files located in
    ./etc/finite_difference_schemes
    in a consoolidated dictionary called FD_schemes_dn

    inputs:
    dn -- (int) derivative number in .dat file containing
          difference coefficients

    LTE -- (int) local truncation error order

    outputs:
    FD_scheme -- (dict) FD scheme equipped with
                  list of weights w and stencil for
                  the specified order dn at specified LTE
    """
    FD_scheme = {}
    rel_path = './../etc/finite_difference_schemes/'
    infile_suffix = '_FD_coefficients.dat'

    infilename = 'f' + str(dn) + '_LTE' + str(LTE) + infile_suffix
    infilepath = rel_path + infilename

    # create empty subdictionary for given key 'dn#'
    FD_scheme['dn' + str(dn)] = {}

    # create empty subdictionary for given key 'LTE#'
    FD_scheme['dn' + str(dn)]['LTE' + str(LTE)] = {}


    # store schemes in subdictionaries of dn
    # which are the keys: 'handedness', 'asymmetry'
    FD_scheme = store_FD_scheme(infilepath,
                             FD_scheme,
                             str(dn),
                             str(LTE))

    return FD_scheme

def store_FD_scheme(infilename,
                    FD_scheme,
                    dn,
                    LTE):
    """reads infile, creates empty subdictionaries
    inside FD_schemes, and stores all data inside
    the corresponding dictionary objects

    inputs:
    infilename -- (str) file name for a single derivative
                  table, e.g. f1_LTE6_FD_coefficients.dat,
                              f2_LTE7_FD_coefficients.dat, ...

    FD_scheme -- (dict) empty dictionary
    dn -- (str) single derivative whose schemes are to be stored

    outputs:
    FD_scheme -- (dict) same dictionary loaded with
                  all schemes for the dn'th derivative
    """
    infile = open(infilename, 'r')
    lines = infile.readlines()

    dn = 'dn' + dn # dn key (str)
    LTE = 'LTE' + LTE # LTE key (str)

    # create empty subdictionaries for each handedness
    FD_scheme[dn][LTE]['forward'] = {}
    FD_scheme[dn][LTE]['central'] = {}
    FD_scheme[dn][LTE]['backward'] = {}

    for i in range(len(lines)):

        if lines[i][0] == 'f':
            handedness = 'forward'
            asymmetry = lines[i].split(' ')[1].strip()

            # create empty subdictionary for given asymmetry
            FD_scheme[dn][LTE][handedness][asymmetry] = {}
        elif lines[i][0] == 'c':
            handedness = 'central'
            asymmetry = lines[i].split(' ')[1].strip()

            # create empty subdictionary for given asymmetry
            FD_scheme[dn][LTE][handedness][asymmetry] = {}
        elif lines[i][0] == 'b':
            handedness = 'backward'
            asymmetry = lines[i].split(' ')[1].strip()

            # create empty subdictionary for given asymmetry
            FD_scheme[dn][LTE][handedness][asymmetry] = {}
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
            FD_scheme[dn][LTE][handedness][asymmetry]['w'] = w
            FD_scheme[dn][LTE][handedness][asymmetry]['stencil'] = stencil
 
            # reset switch for next scheme dict key
            store_data_switch = 0

        else:
            # let i increment until it reaches i + numlines and repeat
            pass

    return FD_scheme

def assemble_finite_difference_weight_matrix_single_derivative(Nx, dn = 1, LTE = 6):
    """Assembles a matrix corresponding to the weights of in
    the finite difference computation of derivatives, i.e.
    assembles the weight matrix W, in

        Wf = df

    where f and df are vectors of length z.N.

    inputs:
    sim_params -- (dict) simulation parameters
    z -- (instance) phase space variable

    outputs:
    Wz -- (ndarray, ndim=3) Wz[dn, z.N, z.N] where each 2D matrix Wz[dn,:,:]
          is the weight matrix W corresponding to the dn-th derivative
    """
    imax = Nx - 1
    W = np.zeros([Nx, Nx])

    FD_scheme_dn = read_FD_scheme(dn, LTE)
    FD_scheme = FD_scheme_dn['dn' + str(dn)]['LTE' + str(LTE)]
    stencil_size = LTE + dn
    stencil_center = stencil_size // 2

    for i in range(Nx):
        if i < stencil_center:
            handedness = 'forward'
            asymmetry = str(i)
        elif imax - i < stencil_center:
            handedness = 'backward'
            asymmetry = str(imax - i)
        else:
            if np.mod(stencil_size,2) == 1:
                handedness = 'central'
                asymmetry = str(0)
            else:
                handedness = 'forward'
                asymmetry = str(stencil_center - 1)

        w = FD_scheme[handedness][asymmetry]['w']
        stencil = FD_scheme[handedness][asymmetry]['stencil']

        W[i, i + np.array(stencil)] = w # load all weights at once into W_dn

    return W
