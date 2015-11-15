import numpy as np
import linecache
import scipy.misc

def inputfile(filename):
    """Reads the input file and returns a dictionary containing
    simulation parameters, sim_params. Function called from
    directory ./DECSKS, input files located in ./DECSKS/etc

    inputs:
    filename -- (str) string with filename containing sim params

    outputs:
    sim_params -- (dict) dictionary containing simulation parameters
                         as well as a dictionary of splitting coeffs
                         needed for chosen split scheme
    """
    rel_path = './etc/'
    infile = open(filename, 'r')
    lines = infile.readlines()

    HOC = lines[6][lines[6].find('=')+1:].strip()
    HOC = HOC.upper()
    derivative_method = '.'.join(('DECSKS.lib.derivatives', HOC.lower()))
    split_function_handle = '_'.join(('DECSKS.lib.split', HOC.lower()))
    split_function_handle = '.'.join((split_function_handle, 'scheme'))

    N = eval(lines[7][lines[7].find('=')+1:].strip())
    print "%s based high order corrections, LTE[CS] = %d" % (HOC, N+1)

    WindowedFilter = lines[8][lines[8].find('=')+1:].strip()
    WindowedFilter = WindowedFilter.upper()

    Nx = eval(lines[16][lines[16].find('=')+1:].strip())
    ax = eval(lines[17][lines[17].find('=')+1:].strip())
    bx = eval(lines[18][lines[18].find('=')+1:].strip())

    Ny = eval(lines[22][lines[22].find('=')+1:].strip())
    ay = eval(lines[23][lines[23].find('=')+1:].strip())
    by = eval(lines[24][lines[24].find('=')+1:].strip())

    Nz = eval(lines[28][lines[28].find('=')+1:].strip())
    az = eval(lines[29][lines[29].find('=')+1:].strip())
    bz = eval(lines[30][lines[30].find('=')+1:].strip())

    Nvx = eval(lines[34][lines[34].find('=')+1:].strip())
    avx = eval(lines[35][lines[35].find('=')+1:].strip())
    bvx = eval(lines[36][lines[36].find('=')+1:].strip())

    Nvy = eval(lines[40][lines[40].find('=')+1:].strip())
    avy = eval(lines[41][lines[41].find('=')+1:].strip())
    bvy = eval(lines[42][lines[42].find('=')+1:].strip())

    Nvz = eval(lines[46][lines[46].find('=')+1:].strip())
    avz = eval(lines[47][lines[47].find('=')+1:].strip())
    bvz = eval(lines[48][lines[48].find('=')+1:].strip())

    Nt = eval(lines[52][lines[52].find('=')+1:].strip())
    at = eval(lines[53][lines[53].find('=')+1:].strip())
    bt = eval(lines[54][lines[54].find('=')+1:].strip())
    T = eval(lines[55][lines[55].find('=')+1:].strip())

    # Boundary condition dictionary storage
    # BC['z']['LBC'] and BC['z']['RBC'] give the BC for phase space var z
    BC = read_boundary_conditions(lines) # pass already read-in lines from params.dat

    # the following list contains strings identifying all evolved phase
    # space variables, a subset of ['x', 'y', 'z', 'vx', 'vy', 'vz']
    phasespace_vars = lines[57][lines[57].find(':')+1:].strip().split(',')

    total_dims = []
    active_dims = []
    # strip all whitespace in each entry
    for var in range(len(phasespace_vars)):
        phasespace_vars[var] = phasespace_vars[var].strip()
        total_dims.append(eval('N' + phasespace_vars[var]))

        if ( (BC[phasespace_vars[var]]['lower'] == 'periodic') and (BC[phasespace_vars[var]]['upper'] == 'periodic') ):
            active_dims.append(eval('N' + phasespace_vars[var]) - 1)
        else:
            active_dims.append(eval('N' + phasespace_vars[var]))

    numdims = len(phasespace_vars)

    densities_list = lines[58][lines[58].find('=')+1:].strip().split(',')
    if len(densities_list) == 2: # if two species return dictionary of strings
        density = {}
        density['electron'] = densities_list[0]
        density['electron'] = density['electron'].lower()
        density['ion'] = densities_list[1]
        density['ion'] = density['ion'].lower()

    elif len(densities_list) == 1: # if one species return a string
        density = densities_list[0]

    split_scheme = lines[77][lines[77].find('=')+1:].strip()
    split_scheme = split_scheme.upper()
    print "using %s split scheme (note: only activated if more than 1D)" % split_scheme

    # splitting input filepath setup
    filename  = lines[78][lines[78].find(':')+1:].strip()
    filepath = rel_path + filename

    # get splitting coefficients for chosen scheme
    if split_scheme is not None:
        splitting = splitting_coefficients(filepath, split_scheme)
    else:
        splitting = None

    plot_params = plot_parameters(lines) # pass already read-in lines of params.dat

    # Table of Bernoulli numbers dat filepath setup
    filename = 'Table_of_Bernoulli_numbers.dat'
    filepath = rel_path + filename
    Bernoulli_numbers = Bernoulli(filepath)

    record_outputs = lines[97][lines[97].find(':')+1:].strip()
    record_outputs = record_outputs.lower()

    if record_outputs == 'yes':
        # output filepath setup
        filename = lines[99][lines[99].find(':')+1:].strip()
        filepath = rel_path + filename
        outfiles = output_files(filepath) # dictionary of opened files
    else:
        outfiles = None

    if HOC == 'FD':
        FD_schemes = read_FD_schemes(N)
        FD_scheme_dn1 = read_FD_scheme(1,6) # LTE = 6 currently, usage is
                                            # for 6th order FD Poisson
                                            # solver
        Xi = None
        xi = None

    else:
        FD_schemes = None
        FD_scheme_dn1 = None

        # Xi[q,i,j] = spectral differentiation matrix = (1j * xi[i,j]) ** q
        # xi = wave number vector = (x.N,)
        Xi, xi = assemble_spectral_derivative_operator(total_dims, active_dims,
                                                       ax, bx, avx, bvx,
                                                       N)

    # MISC STORAGE

    # for correctors c, need I_alternating = the negative of identity
    # matrix where each row is raised to the power of its row

    I_alternating = np.diag( (-np.ones(N))  ** np.arange(N) )

    # A matrices for Bernoulli number storage and matrix HOC application

    A_pos, A_neg = np.zeros([N,N]), np.zeros([N,N])
    for i in range(N):
        for j in range(i+1):
            A_pos[i,j] = Bernoulli_numbers[i-j] / scipy.misc.factorial(i-j)
            if (i - j) == 1:
                A_neg[i,j] = -A_pos[i,j]
            else:
                A_neg[i,j] = A_pos[i,j]

    A_matrix = {}
    # dictionary container
    # allow dictionary access to relevant matrix of Bernoulli numbers
    # by operating with str(int(np.sign(CFL.frac)))

    A_matrix['1'] = A_pos
    A_matrix['0'] = A_pos
    A_matrix['-1'] = A_neg

    # 6th order FD PBC Poisson solver matrices dictionary
    # (stored as keys 'D' [difference matrix] and 'B' [inhomogeneity])

    Poisson_6th_order_PBC_FD_solver_matrices = assemble_Poisson_6th_order_PBC_FD_solver_matrices(Nx, BC)

    sim_params = dict(
        N = N, HOC = HOC,
        derivative_method = derivative_method,
        split_function_handle = split_function_handle,
        WindowedFilter = WindowedFilter,
        Nx = Nx, ax = ax, bx = bx,
        Ny = Ny, ay = ay, by = by,
        Nz = Nz, az = az, bz = bz,
        Nvx = Nvx, avx = avx, bvx = bvx,
        Nvy = Nvy, avy = avy, bvy = bvy,
        Nvz = Nvz, avz = avz, bvz = bvz,
        Nt = Nt, at = at, bt = bt, T = T,
        phasespace_vars = phasespace_vars,
        numdims = numdims,
        active_dims = active_dims,
        total_dims = total_dims,
        density = density,
        split_scheme = split_scheme,
        splitting = splitting,
        plot_params = plot_params,
        Bernoulli_numbers = Bernoulli_numbers,
        record_outputs = record_outputs,
        outfiles = outfiles,
        FD_schemes = FD_schemes,
        FD_scheme_dn1 = FD_scheme_dn1,
        BC = BC,    # boundary conditions on all phase space variables
        I_alternating = I_alternating, # identity matrix with alternating signs according to row, used in computing correctors c
        A_matrix = A_matrix,     # Matrices of Bernoulli numbers for HOC
        Xi = Xi, # spectral differentiation operator matrix (1j*xi[i,j]) ** q
        xi = xi, # wave number vector
        Poisson_6th_order_PBC_FD_solver_matrices = Poisson_6th_order_PBC_FD_solver_matrices
        )

    infile.close()

    return sim_params

def splitting_coefficients(filepath, split_scheme):
    """Reads in the splitting coefficient for the specified
    scheme in input file (e.g. params.dat)
    inputs:

    split_scheme -- (str) designates which split scheme
    filepath -- (file) input file with splitting coeffs rel. path

    output:
    split_coeffs -- (dict) splitting coefficients for specified scheme

    usage:
    in split_schemes module, store and call as

    splitting = sim_params['splitting'] # grab from sim_params dict
    coeff = splitting['order']['coeffs'] = a, b, a, b, ...
    stage = splitting['order']['stages'] = 1, 1, 2, 2, ...
    access ith coefficient by
    splitting[ coeff[i]][int(stage[i])]
    """
    infile = open(filepath, 'r')
    lines = infile.readlines()
    infile.close()

    if split_scheme == 'LF2':
        coeffs = lines[6].strip().split(', ')
        stages = lines[9].strip().split(', ')
        a1 = eval(lines[14][lines[14].find('=')+1:].strip())
        a2 = eval(lines[15][lines[15].find('=')+1:].strip())
        b1 = eval(lines[16][lines[16].find('=')+1:].strip())
        b2 = eval(lines[17][lines[17].find('=')+1:].strip())

        order = dict(coeffs = coeffs, stages = stages)
        splitting = dict(order = order,
                            a = [None, a1, a2],
                            b = [None, b1, b2])

    elif split_scheme == 'Y4':
        coeffs = lines[23].strip().split(', ')
        stages = lines[26].strip().split(', ')
        a1 = eval(lines[31][lines[31].find('=')+1:].strip())
        a2 = eval(lines[32][lines[32].find('=')+1:].strip())
        a3 = eval(lines[33][lines[33].find('=')+1:].strip())
        a4 = eval(lines[34][lines[34].find('=')+1:].strip())
        b1 = eval(lines[36][lines[36].find('=')+1:].strip())
        b2 = eval(lines[37][lines[37].find('=')+1:].strip())
        b3 = eval(lines[38][lines[38].find('=')+1:].strip())
        b4 = eval(lines[39][lines[39].find('=')+1:].strip())

        order = dict(coeffs = coeffs, stages = stages)
        splitting = dict(order = order,
                            a = [None, a1, a2, a3, a4],
                            b = [None, b1, b2, b3, b4])

    elif split_scheme == 'O6-4':
        coeffs = lines[45].strip().split(', ')
        stages = lines[48].strip().split(', ')
        a1 = eval(lines[53][lines[53].find('=')+1:].strip())
        a2 = eval(lines[54][lines[54].find('=')+1:].strip())
        a3 = eval(lines[55][lines[55].find('=')+1:].strip())
        a4 = eval(lines[56][lines[56].find('=')+1:].strip())
        b1 = eval(lines[58][lines[58].find('=')+1:].strip())
        b2 = eval(lines[59][lines[59].find('=')+1:].strip())
        b3 = eval(lines[60][lines[60].find('=')+1:].strip())
        b4 = eval(lines[61][lines[61].find('=')+1:].strip())

        order = dict(coeffs = coeffs, stages = stages)
        splitting = dict(order = order,
                            a = [None, a1, a2, a3, a4],
                            b = [None, b1, b2, b3, b4])

    elif split_scheme == 'O11-6':
        coeffs = lines[67].strip().split(', ')
        coeffs += lines[68].strip().split(', ')

        stages = lines[71].strip().split(', ')
        stages += lines[72].strip().split(', ')

        a1 = eval(lines[78][lines[78].find('=')+1:].strip())
        a2 = eval(lines[79][lines[79].find('=')+1:].strip())
        a3 = eval(lines[80][lines[80].find('=')+1:].strip())
        a4 = eval(lines[81][lines[81].find('=')+1:].strip())
        a5 = eval(lines[82][lines[82].find('=')+1:].strip())
        a6 = eval(lines[83][lines[83].find('=')+1:].strip())
        b1 = eval(lines[85][lines[85].find('=')+1:].strip())
        b2 = eval(lines[86][lines[86].find('=')+1:].strip())
        b3 = eval(lines[87][lines[87].find('=')+1:].strip())
        b4 = eval(lines[88][lines[88].find('=')+1:].strip())
        b5 = eval(lines[89][lines[89].find('=')+1:].strip())
        b6 = eval(lines[90][lines[90].find('=')+1:].strip())

        order = dict(coeffs = coeffs, stages = stages)
        splitting = dict(order = order,
                            a = [None, a1, a2, a3, a4, a5, a6],
                            b = [None, b1, b2, b3, b4, b5, b6])

    elif split_scheme == 'O14-6':
        coeffs = lines[96].strip().split(', ')
        coeffs += lines[97].strip().split(', ')
        coeffs += lines[98].strip().split(', ')

        stages = lines[101].strip().split(', ')
        stages += lines[102].strip().split(', ')
        stages += lines[103].strip().split(', ')

        a1 = eval(lines[110][lines[110].find('=')+1:].strip())
        a2 = eval(lines[111][lines[111].find('=')+1:].strip())
        a3 = eval(lines[112][lines[112].find('=')+1:].strip())
        a4 = eval(lines[113][lines[113].find('=')+1:].strip())
        a5 = eval(lines[114][lines[114].find('=')+1:].strip())
        a6 = eval(lines[115][lines[115].find('=')+1:].strip())
        a7 = eval(lines[116][lines[116].find('=')+1:].strip())
        a8 = eval(lines[117][lines[117].find('=')+1:].strip())
        b1 = eval(lines[119][lines[119].find('=')+1:].strip())
        b2 = eval(lines[120][lines[120].find('=')+1:].strip())
        b3 = eval(lines[121][lines[121].find('=')+1:].strip())
        b4 = eval(lines[122][lines[122].find('=')+1:].strip())
        b5 = eval(lines[123][lines[123].find('=')+1:].strip())
        b6 = eval(lines[124][lines[124].find('=')+1:].strip())
        b7 = eval(lines[125][lines[125].find('=')+1:].strip())

        order = dict(coeffs = coeffs, stages = stages)
        splitting = dict(order = order,
                            a = [None, a1, a2, a3, a4, a5, a6, a7, a8],
                            b = [None, b1, b2, b3, b4, b5, b6, b7])

    return splitting

def plot_parameters(lines):
    """Reads in plot parameters from input file (e.g. params.dat)
    inputs:
    lines -- (list, str) lines from 'params.dat'.readlines()

    output:
    plot_params -- (dict) domain for plot [xmin, xmax, ymin, ymax]
    """

    # lines from filename = 'input_params.dat' in Input(*args) method
    xmin = eval(lines[91][lines[91].find('=')+1:].strip())
    xmax = eval(lines[92][lines[92].find('=')+1:].strip())
    ymin = eval(lines[94][lines[94].find('=')+1:].strip())
    ymax = eval(lines[95][lines[95].find('=')+1:].strip())

    plot_params = dict(xmin = xmin, xmax = xmax,
                       ymin = ymin, ymax = ymax)

    return plot_params

def Bernoulli(filepath):
    """Reads in Bernoulli numbers from data file
    inputs:
    filepath -- (str) relative path to file containing
                    'Table_of_Bernoulli_numbers.dat'

    output:
    B -- (ndarray, ndim=1), numpy 1D array of all Bernoulli
        numbers contained within 'Table_of_Bernoulli_numbers.dat'
    """

    infile = open(filepath, 'r')
    numbers = infile.readlines()
    infile.close()

    B = np.array([eval(number) for number in numbers])
    return B

def output_files(filepath):
    """Reads in output filenames from input file (e.g. params.dat), opens
    all files and returns a dictionary containing all files ready for writing

    inputs:
    filepath -- (str) relative path to file containing
                    'params_output.dat', which holds all output filenames
                    to be used and relative path information

    output:
    outfiles -- (dict) opened output files ready for writing

    """

    infile = open(filepath, 'r')
    lines = infile.readlines()

    rel_path = './'
    rel_path += lines[6][lines[6].find(':')+1:].strip()
    rel_path += lines[7][lines[7].find(':')+1:].strip()

    filename_I1 = lines[9][lines[9].find(':')+1:].strip()
    filename_I2 = lines[10][lines[10].find(':')+1:].strip()
    filename_IW = lines[12][lines[12].find(':')+1:].strip()
    filename_WE = lines[13][lines[13].find(':')+1:].strip()
    filename_CFLx = lines[15][lines[15].find(':')+1:].strip()
    filename_CFLv = lines[16][lines[16].find(':')+1:].strip()
    filename_S = lines[18][lines[18].find(':')+1:].strip()

    filepath_I1 = rel_path + filename_I1
    filepath_I2 = rel_path + filename_I2
    filepath_IW = rel_path + filename_IW
    filepath_WE = rel_path + filename_WE
    #    filepath_CFLx = rel_path + filename_CFLx
    #    filepath_CFLv = rel_path + filename_CFLv
    filepath_S = rel_path + filename_S

    outfile_I1 = open(filepath_I1, 'w')
    outfile_I2 = open(filepath_I2, 'w')
    outfile_IW = open(filepath_IW, 'w')
    outfile_WE = open(filepath_WE, 'w')
    #    outfile_CFLx = open(filepath_CFLx, 'w')
    #    outfile_CFLv = open(filepath_CFLv, 'w')
    outfile_S = open(filepath_S, 'w')

    outfiles = dict(I1 = outfile_I1,
                    I2 = outfile_I2,
                    IW = outfile_IW,
                    WE = outfile_WE,
    #                    CFLx = outfile_CFLx,
    #                    CFLv = outfile_CFLv,
                    S = outfile_S)

    return outfiles

def store_FD_schemes(infilename,
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

        if lines[i][0] == 'f':
            handedness = 'forward'
            asymmetry = lines[i].split(' ')[1].strip()

            # create empty subdictionary for given asymmetry
            FD_schemes[dn][handedness][asymmetry] = {}
        elif lines[i][0] == 'c':
            handedness = 'central'
            asymmetry = lines[i].split(' ')[1].strip()

            # create empty subdictionary for given asymmetry
            FD_schemes[dn][handedness][asymmetry] = {}
        elif lines[i][0] == 'b':
            handedness = 'backward'
            asymmetry = lines[i].split(' ')[1].strip()

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

def read_FD_schemes(N):
    """store all finite difference schemes from
    tables generated in dat files located in
    ./etc/finite_difference_schemes
    in a consoolidated dictionary called FD_schemes

    inputs:
    dn_max -- (int) maximum derivative in .dat files
              should correspond to same dn_max as
              in tables generated

    outputs:
    FD_schemes -- (dict) all FD schemes equipped with
                  list of weights w and stencil
    """
    FD_schemes = {}
    rel_path = './etc/finite_difference_schemes/'
    infile_suffix = '_FD_coefficients.dat'

    for dn in range(1,N):
        infilename = 'f' + str(dn) + infile_suffix
        infilepath = rel_path + infilename

        # create empty subdictionary for given dn
        FD_schemes['dn' + str(dn)] = {}

        # store schemes in subdictionaries of dn
        # which are the keys: 'handedness', 'asymmetry'
        FD_schemes = store_FD_schemes(infilepath,
                                 FD_schemes,
                                 str(dn))

    return FD_schemes

#-----------------------------------------------------------#
# the following are for the sole purpose of reading a single
# dn at given LTE, the original purpose was to have dn1, LTE6
# so that dphi can be computed

# TODO clean this up by combining above routines with this one
# TODO in a general function call and looping only if needed

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
    rel_path = './etc/finite_difference_schemes/'
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

def read_boundary_conditions(lines):
    """Assembles a dictionary permitting access to strings
    indicating the type of boundary condition on each
    phase space variable

        e.g. BC['z']['upper'] gives the string specified in
             etc/params.dat under the upper boundary condition
             for the variable z

    inputs:
    lines -- (list of strings) the read-in lines from etc/params.dat

    outputs:
    BC -- (dict) boundary conditions for every phase space variable
          at every boundary, see access example above

    """
    BC = {}

    BC['x'] = {}
    BC['y'] = {}
    BC['z'] = {}

    BC['vx'] = {}
    BC['vy'] = {}
    BC['vz'] = {}

    BC['x']['lower'] = lines[19][lines[19].find('=')+1:].strip()
    BC['x']['upper'] = lines[20][lines[20].find('=')+1:].strip()

    BC['y']['lower'] = lines[25][lines[25].find('=')+1:].strip()
    BC['y']['upper'] = lines[26][lines[26].find('=')+1:].strip()

    BC['z']['lower'] = lines[31][lines[31].find('=')+1:].strip()
    BC['z']['upper'] = lines[32][lines[32].find('=')+1:].strip()

    BC['vx']['lower'] = lines[37][lines[37].find('=')+1:].strip()
    BC['vx']['upper'] = lines[38][lines[38].find('=')+1:].strip()

    BC['vy']['lower'] = lines[43][lines[43].find('=')+1:].strip()
    BC['vy']['upper'] = lines[44][lines[44].find('=')+1:].strip()

    BC['vz']['lower'] = lines[49][lines[49].find('=')+1:].strip()
    BC['vz']['upper'] = lines[50][lines[50].find('=')+1:].strip()

    return BC

def assemble_spectral_derivative_operator(total_dims, active_dims,
                                          ax, bx, avx, bvx,
                                          N):
    """Returns a dictionary Xi with key/value pairs:

        Xi['x'] -- (ndarray, ndim=3, dtype=complex)
        Xi['vx'] -- (ndarray, ndim=3, dtype=complex)

    Each of these matrices correspond to a matrix with entries

      $$Xi = ((Xi)_{q,i,j}) = 1j * (Delta z xi_{i,j})^q$$

    USAGE NOTE: computing Xi * Ff, where Ff is numpy.fft.fft(f)
    and f.shape = (x.N, vx.N) produces the Fourier transform
    of the derivative coefficients $F[d] equiv D$, where
    D[q,i,j] corresponds to the qth order derivative coefficient
    at a phase space location [i,j]. The method
    lib.derivatives.trigonometric3D takes the row-wise inverse
    transform so that the tensor d[q,i,j] is generated.
    """

    Xi = {}
    xi = {}

    Lx = float(bx - ax)
    xwidth = float((bx - ax) / (total_dims[0] - 1))

    Lvx = float(bvx - avx)
    vxwidth = float((bvx - avx) / (total_dims[1] - 1))

    # build x-advection objects
    N_xi, N_cols = active_dims[0], active_dims[1]
    wave_index = np.arange(N_xi)
    xi_x = np.where(wave_index <= N_xi / 2,
              2*np.pi*wave_index / Lx,
              2*np.pi*(wave_index - N_xi) / Lx)

    xi['x'] = xi_x
    xi_2D = np.outer(xi_x, np.ones(N_cols)) # wavenumbers enumerated by row, copies in Nv columns

    dn = np.arange(1,N).reshape(N-1,1,1)
    Xi['x'] = (1j * xwidth * xi_2D) ** dn

    # build v-advection objects
    N_xi, N_rows = active_dims[1], active_dims[0]
    wave_index = np.arange(N_xi)
    xi_vx = np.where(wave_index <= N_xi / 2,
              2*np.pi*wave_index / Lvx,
              2*np.pi*(wave_index - N_xi) / Lvx)

    xi['vx'] = xi_vx
    xi_2D = np.outer(np.ones(N_rows), xi_vx) # wavenumbers enumerated by row, copies in Nv columns

    xi_2D = np.transpose (xi_2D) # dimensions are (vx.N, x.N) now to match v-advection cases
    dn = np.arange(1,N).reshape(N-1,1,1)
    Xi['vx'] = (1j * vxwidth * xi_2D) ** dn

    return Xi, xi

def assemble_Poisson_6th_order_PBC_FD_solver_matrices(Nx, BC):

    # Nx is the number of active nodes in configuration
    if ( (BC['x']['lower'] == 'periodic') and (BC['x']['upper'] == 'periodic') ):
        Nx -= 1

    Poisson_6th_order_PBC_FD_solver_matrices = {}
    D = np.zeros([Nx,Nx])
    for i in range(Nx):
        if i == 0:         # first row
            D[i,i] = -2
            D[i,i+1] = 1
            D[i,-1] = 1

        elif i == Nx - 1: # last row
            D[i,i] = -2
            D[i,i-1] = 1
            D[i,0] = 1
        else:              # interior rows
            D[i,i-1] = 1
            D[i,i] = -2
            D[i,i+1] = 1

    # Assemble FD matrix B
    B = np.zeros([Nx, Nx])
    for i in range(Nx):
        if i == 0:             # first row
            B[i,-2] = -1/240.
            B[i,-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        elif i == 1:           # second row
            B[i,-1] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        elif 1 < i < (Nx - 2): # 2 <= row < = third before last
            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        elif i == (Nx - 2): # second before last row
            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,0] = -1/240.

        elif i == (Nx - 1): # last row
            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,0] = 1/10.
            B[i,1] = -1/240.

    Poisson_6th_order_PBC_FD_solver_matrices['D'] = D
    Poisson_6th_order_PBC_FD_solver_matrices['B'] = B

    return Poisson_6th_order_PBC_FD_solver_matrices
