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
    infile = open(filename, 'r')
    lines = infile.readlines()

    # --------------------------------------------------------------------------
    # Domain specifications

    Nx = eval(lines[15][lines[15].find('=')+1:].strip())
    ax = eval(lines[16][lines[16].find('=')+1:].strip())
    bx = eval(lines[17][lines[17].find('=')+1:].strip())

    Ny = eval(lines[21][lines[21].find('=')+1:].strip())
    ay = eval(lines[22][lines[22].find('=')+1:].strip())
    by = eval(lines[23][lines[23].find('=')+1:].strip())

    Nz = eval(lines[27][lines[27].find('=')+1:].strip())
    az = eval(lines[28][lines[28].find('=')+1:].strip())
    bz = eval(lines[29][lines[29].find('=')+1:].strip())

    Nvx = eval(lines[33][lines[33].find('=')+1:].strip())
    avx = eval(lines[34][lines[34].find('=')+1:].strip())
    bvx = eval(lines[35][lines[35].find('=')+1:].strip())

    Nvy = eval(lines[39][lines[39].find('=')+1:].strip())
    avy = eval(lines[40][lines[40].find('=')+1:].strip())
    bvy = eval(lines[41][lines[41].find('=')+1:].strip())

    Nvz = eval(lines[45][lines[45].find('=')+1:].strip())
    avz = eval(lines[46][lines[46].find('=')+1:].strip())
    bvz = eval(lines[47][lines[47].find('=')+1:].strip())

    Nt = eval(lines[51][lines[51].find('=')+1:].strip())
    T = eval(lines[52][lines[52].find('=')+1:].strip())

    N = eval(lines[58][lines[58].find('=')+1:].strip())

    # --------------------------------------------------------------------------
    # Broadcast notification regarding start of simulation and order of solver

    print "\nStarting 1D1V Vlasov-Poisson simulation"
    print "\nadvection solver: LTE order %d" % (N+1)

    # --------------------------------------------------------------------------
    # Boundary conditions

    # stored as a dictionary of dictionaries, access as
    # BC['z']['upper'] and BC['z']['lower'] for z = {x, y, ...}

    BC = {}
    # main dictionary with key/values {'x' : {'lower' : value, 'upper : value},
    #                                 {'y' : {'lower' : value, 'upper : value},
    #                                 {'z' : {'lower' : value, 'upper : value},
    #                                 {'vx' : {'lower' : value, 'upper : value},
    #                                 {'vy' : {'lower' : value, 'upper : value},
    #                                 {'vz' : {'lower' : value, 'upper : value},


    # subdictionaries with key/values {'lower' : BC_value, and 'upper' : BC_value}
    BC['x'] = {}
    BC['x']['lower'] = lines[18][lines[18].find('=')+1:].strip()
    BC['x']['upper'] = lines[19][lines[19].find('=')+1:].strip()

    BC['y'] = {}
    BC['y']['lower'] = lines[24][lines[24].find('=')+1:].strip()
    BC['y']['upper'] = lines[25][lines[25].find('=')+1:].strip()

    BC['z'] = {}
    BC['z']['lower'] = lines[30][lines[30].find('=')+1:].strip()
    BC['z']['upper'] = lines[31][lines[31].find('=')+1:].strip()

    BC['vx'] = {}
    BC['vx']['lower'] = lines[36][lines[36].find('=')+1:].strip()
    BC['vx']['upper'] = lines[37][lines[37].find('=')+1:].strip()

    BC['vy'] = {}
    BC['vy']['lower'] = lines[42][lines[42].find('=')+1:].strip()
    BC['vy']['upper'] = lines[43][lines[43].find('=')+1:].strip()

    BC['vz'] = {}
    BC['vz']['lower'] = lines[48][lines[48].find('=')+1:].strip()
    BC['vz']['upper'] = lines[49][lines[49].find('=')+1:].strip()

    # --------------------------------------------------------------------------
    # Store number of active gridpoints for every phase space variable
    #
    # Note: for periodic BCs:  Nz_active = Nz - 1
    #       for all other BCs: Nz_active = Nz

    # TODO this is acknowledged as being redundant, but more specific than the lists
    # active_dims vs. total_dims
    if BC['x']['lower'] == 'periodic' and BC['x']['upper'] == 'periodic' and Nx is not None:
        Nx_active  = Nx - 1
    else:
        Nx_active = Nx

    if BC['y']['lower'] == 'periodic' and BC['y']['upper'] == 'periodic' and Ny is not None:
        Ny_active  = Ny - 1
    else:
        Ny_active = Ny

    if BC['z']['lower'] == 'periodic' and BC['z']['upper'] == 'periodic' and Nz is not None:
        Nz_active  = Nz - 1
    else:
        Nz_active = Nz

    if BC['vx']['lower'] == 'periodic' and BC['vx']['upper'] == 'periodic' and Nvx is not None:
        Nvx_active  = Nvx - 1
    else:
        Nvx_active = Nvx

    if BC['vy']['lower'] == 'periodic' and BC['vy']['upper'] == 'periodic' and Nvy is not None:
        Nvy_active  = Nvy - 1
    else:
        Nvy_active = Nvy

    if BC['vz']['lower'] == 'periodic' and BC['vz']['upper'] == 'periodic' and Nvz is not None:
        Nvz_active  = Nvz - 1
    else:
        Nvz_active = Nvz

    # --------------------------------------------------------------------------
    # High order correction (HOC) method applied to each phase space variable

    HOC = {}
    HOC['x'] = lines[68][lines[68].find(':')+1:].strip().upper()
    HOC['y'] = lines[69][lines[69].find(':')+1:].strip().upper()
    HOC['z'] = lines[70][lines[70].find(':')+1:].strip().upper()

    HOC['vx'] = lines[72][lines[72].find(':')+1:].strip().upper()
    HOC['vy'] = lines[73][lines[73].find(':')+1:].strip().upper()
    HOC['vz'] = lines[74][lines[74].find(':')+1:].strip().upper()


    # list of phase space variables used, in etc/params.dat must set unused
    # vars to have Nz as None, z = x, vx, y, ...
     # e.g. in 1D1V, phasespace_vars = ['x', 'vx']
    phasespace_vars = []
    if Nx is not None:
        phasespace_vars.append('x')
    if Ny is not None:
        phasespace_vars.append('y')
    if Nz is not None:
        phasespace_vars.append('z')
    if Nvx is not None:
        phasespace_vars.append('vx')
    if Nvy is not None:
        phasespace_vars.append('vy')
    if Nvz is not None:
        phasespace_vars.append('vz')

    print "will step through %d-dimensional solution in variables: %s" % (len(phasespace_vars), phasespace_vars)
    for var in phasespace_vars:
        print "high order correction method on %s: %s" % (var, HOC[var])

    # for periodic BCs, the number of active dims is not equal to the
    # total number of dims, we evolve "Nz-1" gridpoints, then assign
    # the Nth point by periodicity as equal to the 0th point. Hence,
    # a distinction is needed between active dims and total dims
    # where we note they are identical in all cases but periodic BCs.

    # TODO as mentioned above, this is now a redundant set of total grid points
    # as compared to active grid points. At some point, need to trace where
    # this is actually used in the code and replace or remove it

    # initialize lists
    total_dims = []
    active_dims = []

    # strip all whitespace in each entry
    for var in phasespace_vars:
        total_dims.append(eval('N' + var))

        if ( (BC[var]['lower'] == 'periodic') and (BC[var]['upper'] == 'periodic') ):
            active_dims.append(eval('N' + var) - 1)
        else:
            active_dims.append(eval('N' + var))

    # TODO this is a misleading name, should be numvars
    numdims = len(phasespace_vars)

    # --------------------------------------------------------------------------
    # Initial density specification
    #
    # the following establishes a difference between the number of densities
    # specified in etc/params.dat. Should there be two, the solver is a two
    # species Vlasov solver. If only one, then a cold background will be
    # automatically computed (TODO)


    densities_list = lines[79][lines[79].find(':')+1:].strip().split(', ')
    for i in range(len(densities_list)):
        densities_list[i] = densities_list[i].lower()

    if len(densities_list) == 2: # if two species return dictionary of strings
        density = {}
        density['electron'] = densities_list[0]
        density['electron'] = density['electron'].lower()
        density['ion'] = densities_list[1]
        density['ion'] = density['ion'].lower()

    elif len(densities_list) == 1: # if one species return a string
        density = densities_list[0]
        # TODO compute cold background, store both this and the above
        # in a common dictionary as above for two species.

    # --------------------------------------------------------------------------
    # Split scheme specification

    split_scheme = lines[98][lines[98].find('=')+1:].strip()
    split_scheme = split_scheme.upper()
    print "split scheme: %s" % split_scheme

    # filepath to splitting coefficient tables
    filename  = lines[99][lines[99].find(':')+1:].strip()
    filepath = './etc/' + filename

    # get splitting coefficients for chosen scheme
    if split_scheme is not None:
        splitting = splitting_coefficients(filepath, split_scheme)
    else:
        splitting = None

    # --------------------------------------------------------------------------
    # Plot window specification (used in lib.plots.Setup)

    xmin = eval(lines[113][lines[113].find('=')+1:].strip())
    xmax = eval(lines[114][lines[114].find('=')+1:].strip())
    ymin = eval(lines[116][lines[116].find('=')+1:].strip())
    ymax = eval(lines[117][lines[117].find('=')+1:].strip())

    plot_params = dict(xmin = xmin, xmax = xmax,
                       ymin = ymin, ymax = ymax)

    record_outputs = lines[120][lines[120].find(':')+1:].strip()
    record_outputs = record_outputs.lower()

    if record_outputs == 'yes':
        # output filepath setup
        filename = lines[121][lines[121].find(':')+1:].strip()
        filepath = './etc/' + filename
        outfiles = output_files(filepath) # dictionary of opened files
    else:
        outfiles = None

    # --------------------------------------------------------------------------
    # MISC STORAGE (e.g. stored matrices that are used routinely)
    #
    # dictionaries and matrices relevant for high order correction applications
    #
    # Constructing the finite different weight matricies, W.
    #-------------------------------------------------------
    #    requires: (dict) FD_schemes
    #
    #    Note: FD_schemes is only needed to construct W. W is what is used in
    #          the simulation. Hence, the building routine for FD_schemes
    #          is not optimized, since it happens before the simulation starts
    #          and hence is not a source of repeated computational cost.
    #
    # FD_schemes is a dictionary containing the families of every order derivative
    # needed for the indicated global error N in etc/params.dat, i.e. all schemes
    # of various degrees of asymmetry and handedness. For large N, this can be a
    # very large dictionary, see the function routine read_FD_schemes to see all
    # that gets stored inside. It is used to construct the difference coefficient
    # matrices W (for applying high order corrections). The other scheme
    # FD_scheme_dn1 is used to construct the matrix W_dn1 which is a difference
    # coefficient matrix for the first derivative (dn = 1) at LTE = 6, and used
    # in the finite difference 6th order Poisson solver (PBCs currently only).
    #---------------------------------------------------------------------------
    #
    # initialize all dictionaries whose keys correspond to phase space vars
    # and whose values contain the relevant ndarrays

    Xi = {}
    xi = {}
    W = {}

    # top level check: if any var has FD corrections, store FD_schemes and init W
    if 'FD' in HOC.values():
        # store finite difference schemes
        FD_schemes = read_FD_schemes(N)

    if HOC['x'] == 'FD':
        # first derivative with LTE = 6, used to find dphi = -E after phi is
        # found from a 6th order Poisson solve
        FD_scheme_dn1 = read_FD_scheme(1,6)
        W_dn1_LTE6 = assemble_finite_difference_weight_matrix_const_dn_const_LTE(Nx_active,
                                             FD_scheme_dn1,
                                             dn = 1,
                                             LTE = 6
                                             )

        # TODO if more than one or different spatial dimension
        # TODO than 'x' with FD corrections need to permit access to this
        # TODO dictionary W_dn1_LTE6 and have it be assembled.

    else:
        # else, Fourier Gauss solver is used, no need for this matrix
        W_dn1_LTE6 = None

    # variable-by-variable checks: assemble consistent objects needed
    # for the specified means of HOC from etc/params.dat

    # Note: the following is organized with the expectation that
    # higher dimensional implementations would be stepped through
    # as sets of 2D advection problems,  always paired as z and vz
    # i.e. not as mixed stepthroughs with x paired with vy for example

    for var in phasespace_vars:
        if HOC[var] == 'FD':
            W[var] = assemble_finite_difference_weight_matrix(
                eval('N' + var + '_active'),
                N,
                FD_schemes
                )
        elif HOC[var] == 'FOURIER':
            # ensure the correct number of grid points
            # is passed for the generalized velocity Nvz_active
            # for x,y,z, 'vz' = vx, vy, vz
            # for vx, vy, vz, 'vz' = ax, ay, az, which have
            # the same number of dims as x, y, z, respectively

            if var[0] == 'v':
                Nvz_active = eval('N' + var[1] + '_active')
            else:
                Nvz_active = eval('Nv' + var + '_active')

            Xi, xi = assemble_spectral_derivative_operator(Xi, xi,
                                                          var,
                                                          eval('a' + var),
                                                          eval('b' + var),
                                                          eval('N' + var),
                                                          eval('N' + var + '_active'),
                                                          Nvz_active,
                                                          N)

    # ---------------------------------------------------------------------
    # "Alternating" identity matrix


    # in lib.HOC.correctors, require an N x N  diagonal matrix with entries
    # (-1)^i, where i is the row number, for details see on github
    #
    #    dsirajud/IPython-notebooks/
    #       DECSKS-09 -- array-based implementation recast -- part 1.ipynb
    #
    # section "2D casting of correction coefficients c (vector) -> c (tensor)"

    I_alternating = np.diag( (-np.ones(N))  ** np.arange(N) )

    # obtain Bernoulli numbers (note: list only 23 numbers are listed)
    # for a correction up to global error order N, N-1 Bernoulli numbers
    # are needed. If higher than global error order 22 is desired, additional
    # Bernoulli numbes need to be entered in
    #
    #    etc/Table_of_Bernoulli_numbers.dat
    #

    # Store Bernoulli numbers from dat file etc/Table_of_Bernoulli_numbers.dat
    filename = 'Table_of_Bernoulli_numbers.dat'
    filepath = './etc/' + filename
    Bernoulli_numbers = Bernoulli(filepath)

    # "A" matrices for Bernoulli number storage and matrix HOC application
    # in lib.HOC.Beta_matrix, see notebook on github at
    # dsirajud/IPython-notebooks/
    # DECSKS-09 -- array-based implementation recast -- part 1.ipynb
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


    # ---------------------------------------------------------------------
    # 6th order finite difference Poisson solver for periodic BCs
    # (stored as keys 'D' [difference matrix] and 'B' [inhomogeneity])

    Poisson_6th_order_PBC_FD_solver_matrices = assemble_Poisson_6th_order_PBC_FD_solver_matrices(Nx, BC)

    # TODO specialize right now to just be x, vx. Figure out how to generalize later with higher dimensions
    compute_electric_field_function_handle_prefix = "DECSKS.lib.fieldsolvers.compute_electric_field_"
    compute_electric_field_function_handle = "".join((compute_electric_field_function_handle_prefix, HOC['x'].lower()))

    derivative_method = {}
    derivative_method_prefix = 'DECSKS.lib.derivatives'
    for var in phasespace_vars:
        derivative_method[var] = ".".join((derivative_method_prefix, HOC[var].lower()))

    sim_params = dict(
        N = N, HOC = HOC,
        derivative_method = derivative_method,
        Nx = Nx, ax = ax, bx = bx,
        Ny = Ny, ay = ay, by = by,
        Nz = Nz, az = az, bz = bz,
        Nvx = Nvx, avx = avx, bvx = bvx,
        Nvy = Nvy, avy = avy, bvy = bvy,
        Nvz = Nvz, avz = avz, bvz = bvz,
        Nt = Nt, T = T,
        phasespace_vars = phasespace_vars,
        numdims = numdims,
        active_dims = active_dims,
        total_dims = total_dims,
        density = density,
        split_scheme = split_scheme,
        splitting = splitting,
        plot_params = plot_params,
        record_outputs = record_outputs,
        outfiles = outfiles,
        BC = BC,    # boundary conditions on all phase space variables
        I_alternating = I_alternating, # identity matrix with alternating signs according to row, used in computing correctors c
        A_matrix = A_matrix,     # Matrices of Bernoulli numbers for HOC
        W = W,
        W_dn1_LTE6 = W_dn1_LTE6,
        Xi = Xi, # spectral differentiation operator matrix (1j*xi[i,j]) ** q
        xi = xi, # wave number vector
        Poisson_6th_order_PBC_FD_solver_matrices = Poisson_6th_order_PBC_FD_solver_matrices,
        compute_electric_field_function_handle = compute_electric_field_function_handle # determines if solver is FD or fourier based
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
# derivative of order dn at given LTE, the original purpose
# was to have dn = 1, LTE = 6 for a 6th order Poisson solver

def read_FD_scheme(dn, LTE):
    """store finite difference scheme for dn'th derivative
    from tables generated in dat files located in
    ./etc/finite_difference_schemes
    in a consolidated dictionary called FD_schemes_dn

    inputs:
    dn -- (int) derivative number in .dat file containing
          difference coefficients

    LTE -- (int) local truncation error order

    **Requires the generated dat file for dn = 1, LTE = 6

            etc/finite_difference_schemes/
            f1_LTE6_FD_coefficients.dat

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

def assemble_finite_difference_weight_matrix(zN,
                                             N,
                                             FD_schemes
                                             ):
    """Assembles a matrix corresponding to the weights of in
    the finite difference computation of derivatives, i.e.
    assembles the weight matrix W, giving the difference matrix d
    for the q-th derivative:

        1 / x.width ** q W[q,:,:].dot(f) =  d[q,:,:]

                                    i.e. W are the difference
                                    coefficients, which do not
                                    contain the width of the
                                    abscissa value, e.g. x.width

    where f and df are vectors of length z.N in the 1D case.

    inputs:
    zN -- (int) number of active grid points for the phase sapce variable z
    N -- (int) global error on the advection algorithm, specified in etc/params.dat
    FD_schemes -- (dict) dictionary containing all schemes for all dn, handedness,
        and asymmetry
    outputs:
    Wz -- (ndarray, ndim=3) Wz[dn, zN, zN] where each 2D matrix Wz[dn,:,:]
          is the weight matrix W corresponding to the dn-th derivative in
          the context of the above equation.
    """
    imax = zN - 1
    Wz = np.zeros([N, zN, zN]) # includes zeroeth order derivative container (not used)
    # i.e. Wz[q,:,:] is for the q-th derivative with this dummy zero index created

    for dn in range(1,N):
        W_dn = np.zeros([zN, zN])
        p = N - dn     # LTE of scheme on dn-th derivative decreases with dn
                       # given what is needed is LTE[z.width ** dn * dnf] = O(N)
                       # rather than LTE on just the derivative

        # local copy of all schemes pertaining to derivative order dn
        FD_schemes_dn = FD_schemes['dn' + str(dn)]
        stencil_size = p + dn
        stencil_center = stencil_size // 2

        for i in range(zN):

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

            FD_scheme = FD_schemes_dn[handedness][asymmetry]
            w = FD_scheme['w']
            stencil = FD_scheme['stencil']

            W_dn[i, i + np.array(stencil)] = w # load all weights at once into W_dn

        Wz[dn,:,:] = W_dn

    return Wz

def assemble_finite_difference_weight_matrix_const_dn_const_LTE(zN,
                                             FD_scheme_const_dn,
                                             dn = 1,
                                             LTE = 6
                                             ):
    """Assembles a matrix corresponding to the weights of in
    the finite difference computation of derivatives, i.e.
    assembles the weight matrix W, giving the difference matrix d
    for the q-th derivative:

        1 / x.width ** q W[q,:,:].dot(f) =  d[q,:,:]

                                    i.e. W are the difference
                                    coefficients, which do not
                                    contain the width of the
                                    abscissa value, e.g. x.width

    where f and df are vectors of length z.N in the 1D case.

    inputs:
    zN -- (int) number of active grid points for the phase sapce variable z
    N -- (int) global error on the advection algorithm, specified in etc/params.dat
    FD_schemes -- (dict) dictionary containing all schemes for all dn, handedness,
        and asymmetry
    outputs:
    Wz -- (ndarray, ndim=3) Wz[dn, zN, zN] where each 2D matrix Wz[dn,:,:]
          is the weight matrix W corresponding to the dn-th derivative in
          the context of the above equation.
    """
    imax = zN - 1
    W_dn_LTE = np.zeros([zN, zN])

    # local copy of all schemes pertaining to derivative order dn
    FD_scheme = FD_scheme_const_dn['dn' + str(dn)]['LTE' + str(LTE)]
    stencil_size = LTE + dn
    stencil_center = stencil_size // 2

    for i in range(zN):

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

        W_dn_LTE[i, i + np.array(stencil)] = w # load all weights at once into W_dn

    return W_dn_LTE

def assemble_spectral_derivative_operator(Xi, xi,
                                          z_str,
                                          az, bz,
                                          Nz, Nz_active,
                                          Nvz_active,
                                          N):
    """For 2D constructions, e.g. (x, vx). For higher dimensions,
    e.g. 4D (x,y, vx, v_perp) can reuse this with some minor
    changes. For 1D or 3D, a different (but, similar) routine
    needs to be coded. For 3D, the overall stepthrough will need
    to be deconstructed to be a split problem among 2D problems
    and a 1D problem.

    inputs:
    Xi -- (dict) to contain key/values:

        Xi['z'] -- (ndarray, ndim=3, dtype = complex), z = x, vx, ...

        this routine adds the key 'z' to the dictionary. Hence,
        the dictionary passed is at minimum an empty dictionary, but
        in general contains previous keys assigned by previuos calls
        to this same function

    xi -- (dict) contains key/values

        xi['z'] -- (ndarray, ndim=1, dtype = float64), z = x, vx, ...

        this routine adds the key 'z' to the dictionary. Hence,
        the dictionary passed is at minimum an empty dictionary, but
        in general contains previous keys assigned by previuos calls
        to this same function

    z_str -- (str) corresponding to phase space variable z affiliated
        with the objects Xi and xi.

    Nz -- (int) total number of gridpoints
    Nz_active -- (int) total number of active gridpoints for z
    Nvz_active -- (int) total number of active gridpoints for vz
    az -- (num) lower domain bound on z, used to compute width Lz
    bz -- (num) upper domain bound on z
    N -- (int) global error of scheme

    outputs: updates the dictionaries Xi, xi to have
    the key/value pair:

        Xi['z'] -- (ndarray, ndim=3, dtype=complex)
        xi['z'] -- (ndarray, ndim=1, dtype = float64)

    which corresponds to a matrix with entries

      $$Xi = ((Xi)_{q,i,j}) = 1j * (Delta z xi_{i,j})^q$$

    USAGE NOTE: computing Xi * Ff, where Ff is numpy.fft.fft(f)
    and f.shape = (x.N, vx.N) produces the Fourier transform
    of the derivative coefficients $F[d] equiv D$, where
    D[q,i,j] corresponds to the qth order derivative coefficient
    at a phase space location [i,j]. The method
    lib.derivatives.trigonometric3D takes the row-wise inverse
    transform so that the tensor d[q,i,j] is generated.
    """

    # catch any nonsense passes to this function, i.e. z
    # does not have a velocity, hence is not being advected
    if Nvz_active is None:
        return None

    # domain widths
    Lz = float(bz - az) # used in wave number vector, xi
    zwidth = Lz / (Nz - 1) # needed for 3D object, Xi = "(j zwidth xi)**q"

    # build wave vector xi for given z
    wave_index = np.arange(Nz_active)
    xi_z = np.where(wave_index <= Nz_active / 2,
              2*np.pi*wave_index / Lz,
              2*np.pi*(wave_index - Nz_active) / Lz)

    xi[z_str] = xi_z

    # Set up compound matrix Xi.
    # First, copy column vector xi along Nvz_active columns
    xi_2D = np.outer(xi_z, np.ones(Nvz_active))

    # set up vector extending in depth dimension so
    # broadcasting per ** operator produces the expected dims on Xi
    # i.e. Xi.shape = (N, z.N, vz.N)
    dn = np.arange(1,N).reshape(N-1,1,1)

    # with the previously formed objects with carefully chosen dims
    # we generate the required Xi object
    Xi[z_str] = (1j * zwidth * xi_2D) ** dn

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
