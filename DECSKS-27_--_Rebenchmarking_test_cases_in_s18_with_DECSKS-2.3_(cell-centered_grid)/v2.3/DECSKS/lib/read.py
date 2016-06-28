import numpy as np
import linecache
import scipy.misc

# lib.read for DECSKS-2.2 input decks

class InputError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def safe_eval(s):
    try:
        return eval(s)
    except NameError:
        return s.lower()

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

    Ny = eval(lines[19][lines[19].find('=')+1:].strip())
    ay = eval(lines[20][lines[20].find('=')+1:].strip())
    by = eval(lines[21][lines[21].find('=')+1:].strip())

    Nz = eval(lines[23][lines[23].find('=')+1:].strip())
    az = eval(lines[24][lines[24].find('=')+1:].strip())
    bz = eval(lines[25][lines[25].find('=')+1:].strip())

    Nvx = eval(lines[27][lines[27].find('=')+1:].strip())
    avx = eval(lines[28][lines[28].find('=')+1:].strip())
    bvx = eval(lines[29][lines[29].find('=')+1:].strip())

    Nvy = eval(lines[31][lines[31].find('=')+1:].strip())
    avy = eval(lines[32][lines[32].find('=')+1:].strip())
    bvy = eval(lines[33][lines[33].find('=')+1:].strip())

    Nvz = eval(lines[35][lines[35].find('=')+1:].strip())
    avz = eval(lines[36][lines[36].find('=')+1:].strip())
    bvz = eval(lines[37][lines[37].find('=')+1:].strip())

    Nt = eval(lines[39][lines[39].find('=')+1:].strip())
    T = eval(lines[40][lines[40].find('=')+1:].strip())

    N = eval(lines[46][lines[46].find('=')+1:].strip())

   # --------------------------------------------------------------------------
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

    # ==========================================================================
    # Boundary conditions dictionary -- contains dist. function BCs as well as phi

    BC = {}
    BC['f'] = {}
    BC['phi'] = {}

    # BC['f'] = BC dict on distribution function f

    #     BC['f']['x'] = {'lower' : lower_value, 'upper' : upper_value}
    #     BC['f']['y'] = {'lower' : lower_value, 'upper' : upper_value}
    #     BC['f']['z'] = {'lower' : lower_value, 'upper' : upper_value}
    #     BC['f']['vx'] = {'lower' : lower_value, 'upper' : upper_value}
    #     BC['f']['vy'] = {'lower' : lower_value, 'upper' : upper_value}
    #     BC['f']['vz'] = {'lower' : lower_value, 'upper' : upper_value}

    # BC['phi'] = BC dict on electric potential phi

    #     BC['phi']['x'] = {'lower' : lower_value, 'upper' : upper_value}
    #     BC['phi']['y'] = {'lower' : lower_value, 'upper' : upper_value}
    #     BC['phi']['z'] = {'lower' : lower_value, 'upper' : upper_value}
    #     BC['phi']['vx'] = {'lower' : lower_value, 'upper' : upper_value}
    #     BC['phi']['vy'] = {'lower' : lower_value, 'upper' : upper_value}
    #     BC['phi']['vz'] = {'lower' : lower_value, 'upper' : upper_value}
    #
    # subdict objects that give keyword descriptions that match method names in lib.boundaryconditions and lib.fieldsolvers
    # include, for var in phasespace_vars:
    #
    # BC['f'][var]['type'] and BC['phi'][var]['type']
    #
    # these are used to assemble function handle strings that select the corresponding routine needed for the specified BCs


    BC_infilename = './etc/' + lines[106][lines[106].find(':')+1:].strip()
    BC_infile = open(BC_infilename, 'r')
    BC_infile_lines = BC_infile.readlines()

    # DECSKS will throw an error if numbers are inputted as BCs in etc/params.dat

    # strings are stored as lowercase as they are used in an eval statement to access
    # the relevant method in lib.boundaryconditions. e.g. 'absorbing' is accessed as
    # either eval('lib.boundaryconditions.absorbing_lower_boundary') or
    # eval('lib.boundaryconditions.absorbing_upper_boundary') in lib.convect.remap_step

    BC['f']['x'] = {}
    BC['f']['x']['lower'] = safe_eval(BC_infile_lines[53][BC_infile_lines[53].find('=')+1:].strip())
    BC['f']['x']['upper'] = safe_eval(BC_infile_lines[54][BC_infile_lines[54].find('=')+1:].strip())

    BC['f']['y'] = {}
    BC['f']['y']['lower'] = safe_eval(BC_infile_lines[56][BC_infile_lines[56].find('=')+1:].strip())
    BC['f']['y']['upper'] = safe_eval(BC_infile_lines[57][BC_infile_lines[57].find('=')+1:].strip())

    BC['f']['z'] = {}
    BC['f']['z']['lower'] = safe_eval(BC_infile_lines[59][BC_infile_lines[59].find('=')+1:].strip())
    BC['f']['z']['upper'] = safe_eval(BC_infile_lines[60][BC_infile_lines[60].find('=')+1:].strip())

    BC['f']['vx'] = {}
    BC['f']['vx']['lower'] = safe_eval(BC_infile_lines[68][BC_infile_lines[68].find('=')+1:].strip())
    BC['f']['vx']['upper'] = safe_eval(BC_infile_lines[69][BC_infile_lines[69].find('=')+1:].strip())

    BC['f']['vy'] = {}
    BC['f']['vy']['lower'] = safe_eval(BC_infile_lines[71][BC_infile_lines[71].find('=')+1:].strip())
    BC['f']['vy']['upper'] = safe_eval(BC_infile_lines[72][BC_infile_lines[72].find('=')+1:].strip())

    BC['f']['vz'] = {}
    BC['f']['vz']['lower'] = safe_eval(BC_infile_lines[74][BC_infile_lines[74].find('=')+1:].strip())
    BC['f']['vz']['upper'] = safe_eval(BC_infile_lines[75][BC_infile_lines[75].find('=')+1:].strip())

    # make all BCs lowercase strings so they can be used to construct the function strings in lib.boundaryconditions module
    # whose names are all lowercase

    # if an accepted boundary condition synonym as been used, change value to the name it goes by in lib.boundaryconditions
    # check that all inputs for evolved phase space variables are recognized keywords and are compatible with the
    # boundary at which they are indicated
    for var in phasespace_vars:
        for boundary in ['lower', 'upper']:
            BC['f'][var][boundary] = BC['f'][var][boundary].lower()

            if BC['f'][var][boundary] == 'cutoff':
                if var in ['x', 'y', 'z']:
                    print "the following boundary condition was not accepted:\n"
                    print "distribution function %s boundary condition on %s: %s" % (boundary, var, BC['f'][var][boundary].upper())
                    print "\na cutoff condition on a configuration grid makes it unclear how to calculate or specify boundary conditions"
                    print "on the electric potential phi. The interpretation of a cutoff condition is that the numerical grid is a control"
                    print "volume inside an otherwise open system. Since we do not simulate what happens outside the control volume, "
                    print "the nature of a cutoff boundary condition on a configuration variable is inconsistent with the objective of"
                    print "specifying boundary conditions on phi and is not acceptable. The system must be closed (e.g. 'absorbing', 'collector')"
                    print "or periodic on a configuration variable. Please reconsider the context of the intended simulation (e.g. if boundary is desired,"
                    print " the corresponding boundary condition on the distribution function should be set to ABSORBING).\n"
                    raise InputError('a CUTOFF boundary condition on the distribution for a a configuration variable is inconsistent (can only be used on velocity variables).')
                else:
                    pass

            elif BC['f'][var][boundary] == 'collector':
                pass

            elif BC['f'][var][boundary] == 'absorbing':
                pass

            elif BC['f'][var][boundary] == 'symmetry':
                if boundary == 'upper':
                    print "the following boundary condition was not accepted:\n"
                    print "distribution function %s boundary condition on %s: %s" % (boundary, var, BC['f'][var][boundary].upper())
                    print "\nDECSKS only supports LOWER symmetry boundaries."
                    raise NotImplementedError('a symmetric UPPER boundary condition on the distribution function was specified in params_boundaryconditions.dat; however, DECSKS only has functionality to permit lower boundary symmetry.')

                elif boundary == 'lower':
                    print "\nCOURTESY NOTICE TO USER: the boundary condition %s was selected for the distribution function on %s at the %s boundary in params_boundaryconditions.dat; " % (BC['f'][var][boundary].upper(), var, boundary)
                    print "this is a recognized input synonym for a '%s' condition. Changing value stored to BC['f']['%s']['%s'] = '%s'\n" % ('SYMMETRIC', var, boundary, 'SYMMETRIC')
                    print "Please regard any warnings/error messages that cite the keyword '%s' with this change in mind\n" %  ('SYMMETRIC')
                    BC['f'][var][boundary] = 'symmetric'

            elif BC['f'][var][boundary] == 'symmetric':
                if boundary == 'lower':
                    pass
                elif boundary == 'upper':
                    print "the following boundary condition was not accepted:\n"
                    print "distribution function %s boundary condition on %s: %s" % (boundary, var, BC['f'][var][boundary].upper())
                    print "\nDECSKS only supports LOWER symmetry boundaries."
                    raise NotImplementedError('a symmetric UPPER boundary condition on the distribution function was specified in params_boundaryconditions.dat; however, DECSKS only has functionality to permit lower boundary symmetry.')

            elif BC['f'][var][boundary] == 'periodic':
                pass

            else: # inputs do not match any options
                print '\nThe invalid keyword %s was specified in params_boundaryconditions.dat on the variable %s at the %s boundary\n' % (BC['f'][var][boundary].upper(), var, boundary)
                raise InputError('inputs are restricted to those listed as options in params_boundaryconditions.dat')

    # above we have checked for valid input. Next, check for compatible inputs (if 'periodic' is selected, it must be selected for both
    # upper and lower bounds) and store a descriptor that toggles the correct orchestrator
    # function in lib.boundaryconditions module ('periodic' vs. 'nonperiodic')
    for var in phasespace_vars:
        if BC['f'][var]['lower'] == 'periodic' and BC['f'][var]['upper'] == 'periodic':
            BC['f'][var]['type'] = 'periodic'

        # check for invalid inputs (if 'periodic' was set at one boundary, it would need to be set at the opposite boundary as 'periodic' as well since
        # a periodic boundary condition effectively involves both boundaries)
        elif BC['f'][var]['lower'] == 'periodic' and BC['f'][var]['upper'] != 'periodic':
            print "\nThe following boundary conditions specified in params_boundaryconditions.dat:"
            print "\nlower boundary condition on f for the variable %s: %s" % (var, BC['f'][var]['lower'].upper())
            print "upper boundary condition on f for the variable %s: %s" % (var, BC['f'][var]['upper'].upper())

            print "\nare inconsistent. Cannot combine periodic and non-periodic boundary conditions on same variable for distribution function, check inputs in params_boundaryconditions.dat')"

            raise InputError('cannot combine periodic and non-periodic boundary conditions on same variable for distribution function, check inputs in params_boundaryconditions.dat')
        elif BC['f'][var]['lower'] != 'periodic' and BC['f'][var]['upper'] == 'periodic':
            print "\nThe following boundary conditions specified in params_boundaryconditions.dat:"
            print "\nlower boundary condition on f for the variable %s: %s" % (var, BC['f'][var]['lower'].upper())
            print "upper boundary condition on f for the variable %s: %s" % (var, BC['f'][var]['upper'].upper())

            print "\nare inconsistent. Cannot combine periodic and non-periodic boundary conditions on same variable for distribution function, check inputs in params_boundaryconditions.dat')"

            raise InputError('cannot combine periodic and non-periodic boundary conditions on same variable for distribution function, check inputs in params_boundaryconditions.dat')

        else:
        # by this point, we have already checked for consistency among the nonperiodic conditions
        # both bouundaries are non-periodic, and are a combination of: symmetric (lower), collector (lower or upper), absorbing (lower or upper), cutoff (lower or upper)
            BC['f'][var]['type'] = 'nonperiodic'

    distribution_function_boundarycondition_orchestrator_prefix = 'DECSKS.lib.boundaryconditions'

    # create a dictionary of function handles that call either
    # the 'periodic', 'nonperiodic', or 'symmetric' orchestrator in lib.boundaryconditions
    #
    # i.e. we form the string handle for each active variable var:
    #
    # distribution_function_boundarycondition_orchestrator_handle[var] =
    #
    #                 DECSKS.lib.boundaryconditions.periodic
    #                 DECSKS.lib.boundaryconditions.nonperiodic
    #                 DECSKS.lib.boundaryconditions.symmetric

    distribution_function_boundarycondition_orchestrator_handle = {}

    for var in phasespace_vars:
        distribution_function_boundarycondition_orchestrator_handle[var] = ".".join(
            (distribution_function_boundarycondition_orchestrator_prefix, BC['f'][var]['type']))


    # --------------------------------------------------------------------------
    # High order correction (HOC) method applied to each phase space variable

    # store as uppercase

    HOC = {}
    HOC['x'] = safe_eval(lines[56][lines[56].find(':')+1:].strip())
    HOC['y'] = safe_eval(lines[57][lines[57].find(':')+1:].strip())
    HOC['z'] = safe_eval(lines[58][lines[58].find(':')+1:].strip())

    HOC['vx'] = safe_eval(lines[60][lines[60].find(':')+1:].strip())
    HOC['vy'] = safe_eval(lines[61][lines[61].find(':')+1:].strip())
    HOC['vz'] = safe_eval(lines[62][lines[62].find(':')+1:].strip())

    # make all non-None inputs capitalized
    for key in HOC.keys():
        if HOC[key] is not None:
            HOC[key] = HOC[key].upper()
        else:
            pass

    # check for valid inputs
    for key in HOC.keys():
        if HOC[key] is not None:
            if type(HOC[key]) != str:
                raise InputError('A non-string entry was found as a high order correction specification. Only FD or FOURIER are accepted')
            elif HOC[key] != 'FD' and HOC[key] != 'FOURIER':
                print "\nThe following high order correction was specified in params.dat, but is not recognized:"
                print "\nHigh order correction on %s: %s\n" % (key, HOC[key].upper())
                print "only FD and FOURIER are accepted keywords\n"
                raise InputError('An unrecognized high order correction was specified. Only FD or FOURIER are accepted')

            elif HOC[key] == 'FOURIER' and BC['f'][key]['type'] != 'periodic': # Fourier corrections use trigonometric derivatives, which rely on periodicity of the underlying functions
                print "\nThe following boundary conditions specified in params_boundaryconditions.dat:"
                print "\nlower boundary condition on f for the variable %s: %s" % (key, BC['f'][key]['lower'].upper())
                print "upper boundary condition on f fore the variable %s: %s\n\n" % (key, BC['f'][key]['upper'].upper())

                print "are inconsistent with the high order correction specified in params.dat:"
                print "\nhigh order correction on %s: %s\n\n" % (key, HOC[var].upper())

                print "FOURIER high order corrections only make sense for periodic systems (if this is the intention, the BCs on f and phi must be set to PERIODIC in params_boundaryconditions.dat)\n"

                raise InputError('Fourier corrections on a variable only make sense for periodic systems. The boundary conditions on the distribution function were read-in as not periodic for this variable.')
            elif eval('N' + key) is None:
                raise InputError('a variable not involved in the simulation (its number of grid points was specified as None) must also have its high order correction method specified as None. While reading in the input deck, the aforementioned expectation was not met. Please revisit the entries (number of grid points) and high order correction specification.')

    # store lists containing number of total and active gridpoints
    # this is acknowledged as redundant given the above storing as Nx_active, Ny_active,
    # etc., but these objects are used in legacy methods inside DECSKS

    # initialize lists
    total_dims = [] # e.g. in 1D1V this could contain [Nx, Nvx]

    for var in phasespace_vars:
        total_dims.append(eval('N' + var))

    numdims = len(phasespace_vars)
    # --------------------------------------------------------------------------
    # Initial density specification (2 species)

    mu = safe_eval(lines[68][lines[68].find(':')+1:].strip())

    densities_list = lines[69][lines[69].find(':')+1:].strip().split(', ')
    for i in range(len(densities_list)):
        densities_list[i] = densities_list[i].lower()

    if len(densities_list) == 2: # if two species return dictionary of strings
        density = {}
        density['electrons'] = densities_list[0]
        density['electrons'] = density['electrons'].lower()
        density['ions'] = densities_list[1]
        density['ions'] = density['ions'].lower()
        print "\ntwo species simulation with initial densities:\n"
        print "electrons: %s" % density['electrons']
        print "ions:      %s\n" % density['ions']

    # --------------------------------------------------------------------------
    # split scheme specification

    split_scheme = lines[81][lines[81].find('=')+1:].strip()
    split_scheme = split_scheme.upper()
    print "split scheme: %s\n" % split_scheme

    # filepath to splitting coefficient tables
    filename  = lines[82][lines[82].find(':')+1:].strip()
    filepath = './etc/' + filename

    # get splitting coefficients for chosen scheme
    if split_scheme is not None:
        splitting = splitting_coefficients(filepath, split_scheme)
    else:
        splitting = None



    # --------------------------------------------------------------------------
    # check for validity on split scheme vs. boundary conditions
    #
    # i.e. check that if the problem is bounded, the user cannot use a split scheme that has negative time substeps
    #
    #    Schemes with only positive time substeps: LF2
    #    Schemes that contain negative time substeps: Y4, O6-4, O11-6, O14-6
    #

    for var in phasespace_vars:
        if BC['f'][var]['lower'] != 'periodic' and BC['f'][var]['upper'] != 'periodic':
            if split_scheme in ['LF2']:
                pass
            else: # a split scheme that involves negative time substeps has been selected
                print "\nThe following set of user specified information is not accepted by DECSKS:\n"
                print "\nin params.dat, the following was specified:"
                print "split scheme = %s:" % split_scheme
                print "\nand the boundary data was specified in params_boundaryconditions.dat:\n"
                print "distribution function lower boundary condition on %s: %s" % (BC['f'][var]['lower'],var)
                print "distribution function upper boundary condition on %s: %s" % (BC['f'][var]['upper'], var)
                print "\nThe split scheme involves negative time substeps, while the boundary conditions are non-periodic. The BOUNDED Vlasov-Poisson problem is irreversible. A split scheme with negative time substeps can only be used in periodic systems, which correspond to systems of infinite extent\n"
                raise InputError('The split scheme involves negative time substeps, while the boundary conditions are non-periodic. The BOUNDED Vlasov-Poisson problem is irreversible. A split scheme with negative time substeps can only be used in periodic systems, which correspond to systems of infinite extent. To rectify this, the user may wish to select periodic boundary conditions on the distribution function (hence phi).')

    # --------------------------------------------------------------------------
    # Plot window specification (used in lib.plots.Setup)

    xmin = eval(lines[96][lines[96].find('=')+1:].strip())
    xmax = eval(lines[97][lines[97].find('=')+1:].strip())
    ymin = eval(lines[99][lines[99].find('=')+1:].strip())
    ymax = eval(lines[100][lines[100].find('=')+1:].strip())

    plot_params = dict(xmin = xmin, xmax = xmax,
                       ymin = ymin, ymax = ymax)

    record_outputs = lines[103][lines[103].find(':')+1:].strip()
    record_outputs = record_outputs.lower()

    if record_outputs == 'yes':
        # output filepath setup
        print "\ncourtesy notice: selected outputs as indicated in etc/params_output.dat will be written to etc/outputs/\n"
        filename = lines[104][lines[104].find(':')+1:].strip()
        filepath = './etc/' + filename
        outfiles = output_files(filepath) # dictionary of opened files
    else:
        print "\ncourtesy notice: no output data will be written\n"
        outfiles = None

    # --------------------------------------------------------------------------
    # DICTIONARIES AND MATRICES RELEVANT FOR HIGH ORDER CORRECTION APPLICATIONS
    #

    # Constructing the finite different weight matrices, W.
    #-------------------------------------------------------
    #    requires: (dict) FD_schemes
    #
    #    Note: FD_schemes is only needed to construct W. W is what is used in
    #          the simulation. Hence, the building routine for FD_schemes
    #          is not optimized, since it happens before the simulation starts
    #          and is not a source of repeated computational cost.
    #
    # FD_schemes is a dictionary containing the families of every order derivative
    # needed for the indicated global error N in etc/params.dat, i.e. all schemes
    # of various degrees of asymmetry and handedness. For large N, this can be a
    # large dictionary, cf. the function routine read_FD_schemes to see all
    # that gets stored inside. It is used to construct the difference coefficient
    # matrices W (for applying high order corrections). The other scheme
    # FD_scheme_dn1 is used to construct the matrix W_dn1 which is a difference
    # coefficient matrix for the first derivative (dn = 1) at LTE = 6, and used
    # to compute the electric field E = "-dphi" = W_dn1.dot(phi),
    # where dphi is the first derivative# of the electric potential, as calculated by
    # the methods in lib.fieldsolvers package
    #---------------------------------------------------------------------------
    #
    # initialize all dictionaries whose keys correspond to phase space vars
    # and whose values contain the relevant ndarrays

    Xi = {}
    xi = {}
    W = {}

    # top level check: if any var has FD corrections, store FD_schemes and init FD weight matrix W
    # for 6th order first derivative
    if 'FD' in HOC.values():
        # store finite difference schemes
        FD_schemes = read_FD_schemes(N)

    # if FD on a configuration variable, need to differentiate phi to obtain the acceleration a ~ E = -dphi
    if HOC['x'] == 'FD' or HOC['y'] == 'FD' or HOC['z'] == 'FD':
        # first derivative with LTE = 6, used to find dphi = -E after phi is
        # found from a 6th order Poisson solve
        FD_scheme_dn1 = read_FD_scheme(1,6)
        W_dn1_LTE6 = assemble_finite_difference_weight_matrix_const_dn_const_LTE(Nx,
                                             FD_scheme_dn1,
                                             dn = 1,
                                             LTE = 6
                                             )

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
                eval('N' + var),
                N,
                FD_schemes
                )
        elif HOC[var] == 'FOURIER':
            # ensure the correct number of grid points
            # is passed for the generalized velocity Nvz_active
            # for x,y,z, 'vz' = vx, vy, vz
            # for vx, vy, vz, 'vz' = ax, ay, az, which have
            # the same number of dims as x, y, z, respectively
            # this is needed in the routine assemble_spectral_derivative_operator
            # so that the correctly dimensioned 2D arrays are returned

            if var[0] == 'v':
                # if a velocity variable, the velocity of this velocity is an acceleration
                # which has the same dimensions as the corresponding configuration variable
                # e.g. vx has velocity(vx) = ax which has the same dimensions as x
                Nvv = eval('N' + var[1])
            else:
                # if a configuration variable, the velocity is the physical velocity, which
                # must be a coresponding active variable
                # e.g. x has a velocity vx
                Nvv = eval('Nv' + var)

            # The 3D tensor Xi is used to compute trigonometric derivatives
            # by operating on a 2D array of Fourier wave components (transformed
            # row-wise for each column, where as usual the objects have been
            # transpoed if needed so that the variation (x or vx) is along
            # rows, not columns)
            #
            # Fourier transform (derivatives) = Xi * Fourier transform (f)
            #                derivatives = inverse transform (Xi * Fourier(f))
            #
            #
            # the object xi is used in legacy methods in DECSKS (pre-DECSKSv2.0)

            Xi, xi = assemble_spectral_derivative_operator(Xi, xi,
                                                          var,
                                                          eval('a' + var),
                                                          eval('b' + var),
                                                          eval('N' + var),
                                                          Nvv,
                                                          N)

    # ---------------------------------------------------------------------
    # "Alternating" identity matrix

    # in lib.HOC.correctors, require an diagonal matrix with shape = (Nz, Nz)
    # with entries as (-1)^i, where i is the row number, for details see on github
    #
    #    dsirajud/IPython-notebooks/
    #       DECSKS-09 -- array-based implementation recast -- part 1.ipynb
    #
    # section "2D casting of correction coefficients c (vector) -> c (tensor)"

    I_alternating = np.diag( (-np.ones(N))  ** np.arange(N) )

    # ---------------------------------------------------------------------
    # Bernoulli number storage, and forming the matrices A_pos, A_neg

    # obtain Bernoulli numbers (note: only 23 numbers are entered into the dat file ->
    # max global error is 23 - 1 = 22) for a correction up to global error order
    # N, N-1 Bernoulli numbers are needed. If higher than global error order 22 is
    # desired, additional Bernoulli numbes need to be entered in
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
    #
    # the A matrices are matrices containing scaled Bernoulli numbers (normalized by factorials)
    # that also factor in the sign (direction) information of the advecting density packets
    # (the different amounts to all odd coefficients having opposite sign)

    # The A matrices are used in the method lib.HOC.Beta_matrix (used to construct the array of the *magnitudes*
    # of the Nvz sets of N beta coefficients; note that the high order flux is further computed as a sum of
    # products that alternating with sign according to the parity of the derivative number, i.e. alternates signs
    # among odds and evens. These prefactors are applied at the end of the method lib.HOC.correctors by matrix
    # pre-multiplication of the matrix B with the alternating (in sight) identity matrix I formed above)

    # the method lib.HOC.Beta_matrix is called from inside lib.HOC.correctors (used to assemble the 2D array c of correctors)

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

    #--------------------------------------------------------------------------------------------#
    # ELECTRIC POTENTIAL PHI
    #--------------------------------------------------------------------------------------------#

    #--------------------------------------------------------------------------------------------#
    # Boundary conditions BC['phi'] dictionary and dictionary of boundary values, phi_BC
    #
    # BC['phi']['x', 'y', or 'z']['lower' or 'upper'] = string keyword that describes the BC
    # phi_BC['x', 'y', or 'z'] = boundary value vector phi_BC that appears in a Poisson solver
    #--------------------------------------------------------------------------------------------#

    phi_BC = {}
        # keys: 'x', 'y', 'z'
        # values: ndarrays of size eval('N' + var + '_active)

    BC['phi'] = {}
        # keys: 'x', 'y', 'z'
        # values / keys for subdict:  'lower', 'upper'
        #          values for subdict: string keyword that describes the BC at the key specification

    # --------------------------------------------------------------------------
    # PHI BOUNDARY CONDITIONS AND PHI BOUNDARY VALUES VECTORS FOR SOLVER Phi_BC['x', 'y', or 'z']

    # lines read in from boundaryconditions dat file were stored above in BC_infile_lines
    if HOC['x'] == 'FD':
        BC['phi']['x'] = {}
        BC['phi']['x']['lower'] = safe_eval(BC_infile_lines[209][BC_infile_lines[209].find('=')+1:].strip())
        BC['phi']['x']['upper'] = safe_eval(BC_infile_lines[210][BC_infile_lines[210].find('=')+1:].strip())
        phi_BC['x'] = np.zeros(Nx)
    elif HOC['x'] == 'FOURIER': # periodic fourier solver is used, a BC vector is not needed
        phi_BC['x'] = None

    if HOC['y'] == 'FD':
        BC['phi']['y'] = {}
        BC['phi']['y']['lower'] = safe_eval(BC_infile_lines[212][BC_infile_lines[212].find('=')+1:].strip())
        BC['phi']['y']['upper'] = safe_eval(BC_infile_lines[213][BC_infile_lines[213].find('=')+1:].strip())
        phi_BC['y'] = np.zeros(Ny)
    elif HOC['y'] == 'FOURIER': # periodic fourier solver is used, a BC vector is not needed
        phi_BC['y'] = None

    if HOC['z'] == 'FD':
        BC['phi']['z'] = {}
        BC['phi']['z']['lower'] = safe_eval(BC_infile_lines[215][BC_infile_lines[215].find('=')+1:].strip())
        BC['phi']['z']['upper'] = safe_eval(BC_infile_lines[216][BC_infile_lines[216].find('=')+1:].strip())
        phi_BC['z'] = np.zeros(Nz)
    elif HOC['z'] == 'FOURIER': # periodic fourier solver is used, a BC vector is not needed
        phi_BC['z'] = None

    # ensure all inputs stored above in BC['phi'] dict objects are uppercase and recognized
    for var in ['x', 'y', 'z']:
        if var in phasespace_vars:
            if HOC[var] == 'FOURIER':
                pass
            else: # HOC is FD which computes the Lorentz term through a potential phi (Fourier uses the electric field E)

                # LOWER BOUNDARY CHECKS
                if BC['phi'][var]['lower'] is None:
                    raise InputError('a NoneType was specified as a LOWER boundary condition on the electric potential phi for an active variable (a non-NoneType was specified for the number of grid points on this variable). If the variable is not meant to be evolved, set its number of grid points to None')

                elif type(BC['phi'][var]['lower']) != str:
                    raise InputError('a non-string type as a LOWER boundary condition on the electric potential phi for an active variable (a non-NoneType was specified for the number of grid points on this variable). If the variable is not intended to be active, set its number of grid points to None. Otherwise, a recognized string keyword must be specified on the boundary condition on phi for this variable.')

                else:
                    BC['phi'][var]['lower'] = BC['phi'][var]['lower'].upper()

                    if BC['phi'][var]['lower'] not in ['PERIODIC', 'SELF-CONSISTENT', 'SYMMETRIC', 'SYMMETRY', 'BIAS']:
                        print "\nThe following boundary conditions specified in params_boundaryconditions.dat is not a recognized keyword:\n\n"
                        print "lower boundary condition on phi for variable %s: %s" % (var, BC['phi'][var]['lower'].upper())

                        raise InputError('boundary condition indicated on phi is not an accepted keyword option')

                    elif (BC['phi'][var]['lower'] == 'SYMMETRIC' or BC['phi'][var]['lower'] == 'SYMMETRY') and BC['f'][var]['lower'] != 'symmetric':
                        print "\nThe following boundary conditions specified in params_boundaryconditions.dat is:\n\n"
                        print "lower boundary condition on phi for variable %s: %s\n" % (var, BC['phi'][var]['lower'].upper())
                        print "lower boundary condition on f for variable %s: %s" % (var, BC['f'][var]['lower'].upper())
                        print "upper boundary condition on f for variable %s: %s\n" % (var, BC['f'][var]['upper'].upper())

                        print "a SYMMETRIC boundary condition must be specified on both phi and f"
                        # by this point all synonyms have been normalized on BC['f'][var], 'symmetric' corresponds to the symmetry condition
                        raise InputError('a SYMMETRY boundary condition on phi was specified, but a symmetry boundary was not specified on the distribution function f at this same (lower) boundary. A symmetric domain requires a lower boundary condition to be SYMMETRIC on both phi and f.')

                    else:
                        pass

                # UPPER BOUNDARY CHECKS
                if BC['phi'][var]['upper'] is None:
                    raise InputError('a NoneType was specified as an upper boundary condition on the electric potential phi for an active variable (a non-NoneType was specified for the number of grid points on this variable). If the variable is not meant to be evolved, set its number of grid points to None')

                elif type(BC['phi'][var]['upper']) != str:
                    raise InputError('a non-string type as an upper boundary condition on the electric potential phi for an active variable (a non-NoneType was specified for the number of grid points on this variable). If the variable is not intended to be active, set its number of grid points to None. Otherwise, a recognized string keyword must be specified on the boundary condition on phi for this variable.')

                else:
                    BC['phi'][var]['upper'] = BC['phi'][var]['upper'].upper()

                    if BC['phi'][var]['upper'] not in ['PERIODIC', 'SELF-CONSISTENT', 'SYMMETRIC', 'SYMMETRY', 'BIAS']:
                        print "\nThe following boundary condition specified in params_boundaryconditions.dat is not a recognized boundary condition keyword:\n\n"
                        print "upper boundary condition on phi for variable %s: %s\n" % (var, BC['phi'][var]['upper'].upper())

                        raise InputError('boundary condition indicated on phi is not an accepted keyword option')

                    elif BC['phi'][var]['upper'] == 'SYMMETRIC' or BC['phi'][var]['upper'] == 'SYMMETRY':
                        print "\nThe following boundary condition specified in params_boundaryconditions.dat is not available:\n\n"
                        print "upper boundary condition on phi: %s\n" % BC['phi'][var]['upper'].upper()

                        raise NotImplementedError('a SYMMETRY boundary condition on phi as an UPPER boundary is specified in params_boundaryconditions.dat; only lower boundaries can support a symmetry boundary condition.')


                # CHECK FOR CONSISTENCY IN BOUNDARY CONDITIONS BETWEEN BOTH LOWER AND UPPER SPECIFICATIONS
                if BC['phi'][var]['lower'] == 'PERIODIC' and BC['phi'][var]['upper'] != 'PERIODIC':
                    print "\nThe following boundary conditions specified in params_boundaryconditions.dat are inconsistent together:\n\n"
                    print "lower boundary condition on phi for variable %s: %s" % (var, BC['phi'][var]['lower'].upper())
                    print "upper boundary condition on phi for variable %s: %s\n\n" % (var, BC['phi'][var]['upper'].upper())

                    raise InputError('PERIODIC boundary conditions on phi involve both lower and upper boundaries. The read-in of params_boundaryconditions.dat has the lower boundary condition as PERIODIC but the upper boundary condition is NOT. Both boundary conditions on phi must be set to PERIODIC if a periodic plasma is to be simulated.')

                elif BC['phi'][var]['lower'] != 'PERIODIC' and BC['phi'][var]['upper'] == 'PERIODIC':
                    print "\nThe following boundary conditions specified in params_boundaryconditions.dat are inconsistent together:\n\n"
                    print "lower boundary condition on phi for variable %s: %s" % (var, BC['phi'][var]['lower'].upper())
                    print "upper boundary condition on phi for variable %s: %s\n\n" % (var, BC['phi'][var]['upper'].upper())

                    raise InputError('PERIODIC boundary conditions on phi involve both lower and upper boundaries. The read-in of params_boundaryconditions.dat has the upper boundary condition as PERIODIC but the lower boundary condition is NOT. Both boundary conditions on phi must be set to PERIODIC if a periodic plasma is to be simulated.')

                elif BC['phi'][var]['lower'] == 'PERIODIC' and BC['phi'][var]['upper'] == 'PERIODIC':

                    if BC['f'][var]['type'] != 'periodic': # note that validity and consistency checks on inputs for the distribution function have already been done above
                        print "\nThe following boundary conditions specified in params_boundaryconditions.dat are inconsistent together:\n\n"
                        print "lower boundary condition on phi for variable %s: %s" % (var, BC['phi'][var]['lower'].upper())
                        print "upper boundary condition on phi for variable %s: %s\n" % (var, BC['phi'][var]['upper'].upper())
                        print "lower boundary condition on phi for variable %s: %s" % (var, BC['f'][var]['lower'].upper())
                        print "upper boundary condition on phi for variable %s: %s\n" % (var, BC['f'][var]['upper'].upper())
                        print "e.g. periodic boundaries on phi require periodic boundaries on f for the same variable\n"
                        raise InputError('PERIODIC boundary conditions on were specifed consistently for phi in params_boundaryconditions.dat; however, periodic boundary conditions must also be consistently specified on the distribution function. Revisit params_boundaryconditions.dat and ensure that both lower and upper boundaries on the distribution function f and the potential phi are set to PERIODIC if a periodic plasma is intended to be simulated.')
                    elif BC['f'][var]['type'] == 'periodic': # note that validity and consistency checks on inputs for the distribution function have already been done above
                        pass


                # CHECK FOR CONSISTENCY ON PHI BCS WITH HIGH ORDER CORRECTION METHOD SPECIFIED (note we have already checked this against the distribution function BCs)
                # here, we are only checking to see if that BCs on phi aren't periodic, to ensure that HOC is NOT set to fourier (relies on periodicity))
                # the following conditional check asks: "if (BCs on phi are not periodic) AND (HOC is FOURIER)"
                if ((BC['phi'][var]['lower'] == 'PERIODIC' and BC['phi'][var]['upper'] != 'PERIODIC') or (BC['phi'][var]['lower'] != 'PERIODIC' and BC['phi'][var]['upper'] == 'PERIODIC')) and HOC[var] == 'fourier':
                    print "\nThe following boundary conditions specified in params_boundaryconditions.dat are inconsistent with the specified high order correction method in params.dat: \n\n"
                    print "lower boundary condition on phi for variable %s: %s" % (var, BC['phi'][var]['lower'].upper())
                    print "upper boundary condition on phi for variable %s: %s\n\n" % (var, BC['phi'][var]['upper'].upper())
                    print "upper boundary condition on phi for variable %s: %s\n\n" % (var, HOC[var].upper())
                    print "\n\nFourier high order corrections require periodic boundary conditions on both phi and the distribution function f\n"

                    raise InputError('the high order correction is specified as FOURIER; however, the BCs on the electric potential phi are not periodic. FOURIER corrections require PERIODIC BCs on phi and the distribution function as the methods rely on periodicity')

    #--------------------------------------------------------------------------------------------#
    # BIAS values
    #--------------------------------------------------------------------------------------------#

    Bias = {} # this dictionary is created for reading in the bias values, it is not returned
              # in sim_params dict. If a bias condition is set on any boundary, this dictionary
              # assigns its value at that boundary in the vector phi_BC[var], phi_BC[var] is
              # returned (as usual, var = ['x', 'y', 'z'])

    Bias['x'] = {}
    Bias['y'] = {}
    Bias['z'] = {}

    Bias['x']['lower'] = safe_eval(BC_infile_lines[227][BC_infile_lines[227].find('=')+1:].strip())
    Bias['x']['upper'] = safe_eval(BC_infile_lines[228][BC_infile_lines[228].find('=')+1:].strip())
    Bias['y']['lower'] = safe_eval(BC_infile_lines[230][BC_infile_lines[230].find('=')+1:].strip())
    Bias['y']['upper'] = safe_eval(BC_infile_lines[231][BC_infile_lines[231].find('=')+1:].strip())
    Bias['z']['lower'] = safe_eval(BC_infile_lines[233][BC_infile_lines[233].find('=')+1:].strip())
    Bias['z']['upper'] = safe_eval(BC_infile_lines[234][BC_infile_lines[234].find('=')+1:].strip())

    # check for valid inputs on active variables for any boundary that is specified as BIAS
    for var in ['x', 'y', 'z']:
        if var in phasespace_vars:
            if HOC[var] == 'FOURIER':
                pass
            else:
                for boundary in ['lower', 'upper']:
                    if var in phasespace_vars:
                        if BC['phi'][var][boundary] == 'BIAS':
                            if Bias[var][boundary] is None: # if the BC is BIAS but the value input for the BIAS value is None
                                print "\nThe following specifications in params_boundaryconditions.dat are inconsistent:\n"
                                print "%s boundary condition on phi for variable %s: %s" % (boundary, var, BC['phi'][var][boundary].upper())
                                print "%s BIAS value on phi for variable %s: %s\n" % (boundary, var, Bias[var][boundary])
                                print "e.g. if a boundary condition on phi is set to BIAS for a variable, a number must be specifed under BIAS value\n"
                                raise InputError('A phi boundary condition on an active variable (number of grid points on this variable has been set as non-None) has been specified as BIAS; however, the corresponding BIAS value is NoneType. Must be a number.')
                            elif type(Bias[var][boundary]) == str:
                                print "\nThe following specifications in params_boundaryconditions.dat are inconsistent:\n"
                                print "%s boundary condition on phi for variable %s: %s" % (boundary, var, BC['phi'][var][boundary].upper())
                                print "%s BIAS value on phi for variable %s: %s\n" % (boundary, var, Bias[var][boundary])
                                print "e.g. if a boundary condition on phi is set to BIAS for a variable, a number must be specifed under BIAS value\n"

                                raise InputError('A phi boundary condition on an active variable (number of grid points on this variable has been set as non-None) has been specified as BIAS; however, the corresponding BIAS value is str type. Must be a number.')
                            else:
                                pass

    # E is calculated by the following call flow, first an ORCHESTRATOR is called:
    #
    #        E = lib.fieldsolvers.compute_electric_field_fourier  <--- solves with a Gauss' law solver directly
    #
    #                or
    #
    #        E = lib.fieldsolvers.compute_electric_field_fd       <--- solves a Poisson solver for phi, then differentiate to get E
    #
    # which can generally be called by eval operating on string handles that are themselves constructed
    # per 'lib.fieldsolvers.compute_electric_field_' + HOC[var].lower()
    #
    # If a finite difference routine is specified, a Poisson solve must be performed to obtain phi.
    # We call the relevant Poisson solver among the following options (L = lower boundary, U = upper boundary, DBC = Dirichlet BC, NBC = Neumann BC):
    #
    #      Poisson_6th_PBC
    #      Poisson_6th_LDBC_UDBC
    #      Poisson_6th_LDBC_UNBC
    #      Poisson_6th_LNBC_UDBC
    #      Poisson_6th_LDBC_LDBC
    #      Poisson_6th_UDBC_UNBC
    #

    # which are selected based on the boundary conditions the user has supplied in params_boundaryconditions.dat.
    #
    # finally, we compute and return:
    #
    #    E = - 1 / config_var.width * W_dn1_LTE6.dot(phi)
    #

    # --------------------------------------------------------------------------
    # fieldsolver orchestator handle string for electric field (periodic or non-periodic)
    #
    # currently only 1D1V, only one handle needed. When this will be generalized, can make a dict object with keys corresponding
    # to each active configuration variable

    compute_electric_field_orchestrator_handle = {}
    for var in ['x', 'y', 'z']:
        if var in phasespace_vars:
            # dictionary key labels the component of the electric field: 'x', 'y', 'z'
            compute_electric_field_orchestrator_handle[var] = "DECSKS.lib.fieldsolvers.compute_electric_field_" + HOC[var].lower()


    # ---------------------------------------------------------------------
    # initialize dictionaries for wall charge objects

    sigma = {}
    sigma_n = {}

    for var in ['x', 'y', 'z']:
        if var in phasespace_vars:
            sigma_n[var] = {}
            sigma[var] = {}

    # --------------------------------------------------------------------------
    # Dictionary for the specific electric potential phi function solver needed
    # according to the specified boundary conditions on phi

    for var in ['x', 'y', 'z']:
        if var in phasespace_vars:

            if HOC[var] == 'FOURIER':
                pass # uses electric field E, periodic boundary conditions only

            else: # is FD corrections, and electric potential phi in a Poisson solver, can be periodic or other BCs
                BC['phi'][var]['type'] = BC['phi'][var]['lower'] + '_' + BC['phi'][var]['upper']
                if BC['phi'][var]['type'] == 'PERIODIC_PERIODIC':
                    BC['phi'][var]['type'] = 'PBC'

                    if BC['f'][var]['lower'] != 'periodic' and BC['f'][var]['upper'] != 'periodic':
                        raise InputError('A boundary condition on phi was specified as BIAS; however, the corresponding boundary condition on f is not compatible (must be set to absorbing or equivalent synonym)')


                if BC['phi'][var]['type'] == 'BIAS_BIAS':
                    BC['phi'][var]['type'] = 'LDBC_UDBC'

                    # Dirichlet condition, phi = BIAS value, see notebook s24 for reason for factor of -2.0 in the derivation (from interpolation to the half-integer index)
                    # if grounded (phi = 0 at wall), the factor of -2.0 is of no consequence so this is general.
                    phi_BC[var][0] = -2.0 * float(Bias[var]['lower'])
                    # Dirichlet condition, phi = BIAS value
                    phi_BC[var][-1] = -2.0 * float(Bias[var]['upper'])

                    if BC['f'][var]['lower'] != 'absorbing' or BC['f'][var]['upper'] != 'absorbing': # all synonyms for 'absorbing' (except 'collector') have been seen by this point, and if encountered changed to 'absorbing'
                        raise InputError('A boundary condition on phi was specified as BIAS; however, the corresponding boundary condition on f is not compatible (must be set to absorbing or equivalent synonym)')

                elif BC['phi'][var]['type'] == 'BIAS_SELF-CONSISTENT':
                    BC['phi'][var]['type'] = 'LDBC_UNBC'

                    # Dirichlet condition, phi = BIAS value
                    phi_BC[var][0] = float(Bias[var]['lower'])
                    # Neumann condition, dphi = sigma_upper, translates to phi_BC[-1] = -6 var.width * sigma_upper (see https://github.com/dsirajud/IPython-notebooks/DECSKS-04...ipynb for details)
                    # phi_BC[-1] = - 6 * var.width * sim_params['sigma'][var]['upper'], changes with time step

                    if BC['f'][var]['lower'] != 'absorbing': # all synonyms for 'absorbing' (except 'collector') have been seen by this point, and if encountered changed to 'absorbing'
                        raise InputError('A lower boundary condition on phi was specified as BIAS; however, the corresponding boundary condition on f is not compatible (must be set to absorbing or equivalent synonym)')

                    if BC['f'][var]['upper'] == 'collector': # all synonyms for 'absorbing' (except 'collector') have been seen by this point, and if encountered changed to 'absorbing'
                        # initialize wall charge densities, sigma for the collector (f) /self-consistent (phi) conditions
                        sigma[var]['upper'] = 0    # initialize to zero charge at time zero
                        sigma_n[var]['upper'] = np.zeros(Nt + 1)  # this was put in at one point for plotting wall charge vs. time
                    else:
                        print "\nThe following boundary conditions specified in params_boundaryconditions.dat are inconsistent together:\n\n"
                        print "upper boundary condition on phi for variable %s: %s\n" % (var, BC['phi'][var]['upper'].upper())
                        print "upper boundary condition on f for variable %s: %s\n" % (var, BC['f'][var]['upper'].upper())
                        print "\ne.g. an upper boundary condition on phi as SELF-CONSISTENT must have the upper boundary condition on f as COLLECTOR"
                        print "\ne.g. an upper boundary condition on f as ASBORBING must have the upper boundary condition on phi as BIAS\n"

                        raise InputError('An upper boundary condition on phi was specified as SELF-CONSISTENT; however, the corresponding boundary condition on f is not compatible (must be set to collector)')

                elif BC['phi'][var]['type'] == 'SELF-CONSISTENT_BIAS':
                    BC['phi'][var]['type'] = 'LNBC_UDBC'

                    # Neumann condition, dphi = -sigma_lower, translates to phi_BC[0] = -6 var.width * sigma_lower (see https://github.com/dsirajud/IPython-notebooks/DECSKS-04...ipynb for details)
                    #phi_BC[var][0] = - 6 * var.width * sim_params['sigma'][var]['lower'], changes with time step
                    # Dirichlet condition, phi = BIAS value
                    phi_BC[var][-1] = float(Bias[var]['upper'])

                    # check upper boundary
                    if BC['f'][var]['upper'] == 'absorbing': # all synonyms for 'absorbing' (except 'collector') have been seen by this point, and if encountered changed to 'absorbing'
                        pass
                    else:
                        print "\nThe following boundary conditions specified in params_boundaryconditions.dat are inconsistent together:\n\n"
                        print "upper boundary condition on phi for variable %s: %s\n" % (var, BC['phi'][var]['upper'].upper())
                        print "upper boundary condition on f for variable %s: %s\n\n" % (var, BC['f'][var]['upper'].upper())
                        print "\ne.g. an upper boundary condition set on phi as BIAS must have the upper boundary condition on f as ABSORBING\n"

                        raise InputError('An upper boundary condition on phi was specified as BIAS; however, the corresponding boundary condition on f is not compatible (must be set to absorbing or equivalent synonym)')

                    # check lower boundary
                    if BC['f'][var]['lower'] == 'collector': # all synonyms for 'absorbing' (except 'collector') have been seen by this point, and if encountered changed to 'absorbing'
                        # initialize wall charge densities, sigma for the collector (f) /self-consistent (phi) conditions
                        sigma[var]['lower'] = 0    # initialize to zero charge at time zero
                        sigma_n[var]['lower'] = np.zeros(Nt + 1)  # this was put in at one point for plotting wall charge vs. time
                    else:
                        print "\nThe following boundary conditions specified in params_boundaryconditions.dat are inconsistent together:\n\n"
                        print "lower boundary condition on phi: %s" % BC['phi'][var]['lower'].upper()
                        print "lower boundary condition on f: %s\n" % BC['f'][var]['lower'].upper()
                        print "\ne.g. an lower boundary condition set on phi as SELF-CONSISTENT must have the lower boundary condition on f as COLLECTOR"
                        print "e.g. an lower boundary condition set on f as ABSORBING must have the lower boundary condition on phi as BIAS"
                        print "e.g. an lower boundary condition set on f as PERIODIC requires the upper boundary on f to be PERIODIC as well as both lower and upper boundary conditions on phi to be set to PERIODIC\n"
                        raise InputError('A lower boundary condition on phi was specified as SELF-CONSISTENT; however, the corresponding boundary condition on f is not compatible (must be set to collector if self-consistent boundary potentials are desired). Equivalently, phi is not compatible with f (e.g. if periodic boundaries on f were desired, the potential must also be periodic)')

                elif BC['phi'][var]['type'] == 'SYMMETRIC_BIAS' or BC['phi'][var]['type'] == 'SYMMETRY_BIAS':
                    BC['phi'][var]['type'] = 'LNBC_UDBC'

                    # Neumann condition, dphi = 0 for symmetry
                    phi_BC[var][0] = 0.
                    # Dirichlet condition, phi = BIAS value
                    phi_BC[var][-1] = float(Bias[var]['upper'])

                    if BC['f'][var]['upper'] != 'absorbing': # all synonyms for 'absorbing' (except 'collector') have been seen by this point, and if encountered changed to 'absorbing'
                        print "\nThe following boundary conditions specified in params_boundaryconditions.dat are inconsistent together:\n\n"
                        print "upper boundary condition on phi: %s" % BC['phi'][var]['upper'].upper()
                        print "upper boundary condition on f: %s\n\n" % BC['f'][var]['upper'].upper()
                        print "\ne.g. an upper boundary condition set on phi as BIAS must have the upper boundary condition on f as ABSORBING\n "
                        raise InputError('An upper boundary condition on phi was specified as BIAS; however, the corresponding boundary condition on f is not compatible (must be set to absorbing or equivalent synonym)')


                elif BC['phi'][var]['type'] == 'SYMMETRIC_SELF-CONSISTENT' or BC['phi'][var]['type'] == 'SYMMETRY_SELF-CONSISTENT':
                    BC['phi'][var]['type'] = 'LDBC_LNBC'

                    # We default to a LDBC_LNBC solver, both boundary conditions on left edge, entries 0 (Dirichlet) and 1 (Neumann)
                    # cf. DECSKS-04 notebook for more details:
                    #
                    #    https://github.com/dsirajud/IPython-notebooks/DECSKS-04...ipynb
                    #
                    # Dirichlet condition, set reference potential phi = 0
                    phi_BC[var][0] = 0. # reference potential set to zero
                    # Neumann condition, dphi = 0 for symmetry
                    phi_BC[var][1] = 0.


                    if BC['f'][var]['upper'] == 'collector': # all synonyms for 'absorbing' (except 'collector') have been seen by this point, and if encountered changed to 'absorbing'
                        # initialize wall charge densities, sigma for the collector (f) /self-consistent (phi) conditions
                        # By virtue of the setup, the above enforcements on the lower boundary ensures this unenforced upper Neumann BC is
                        # satisfied automatically given the relationship that Neumann BCs are fixed by due to the Poisson equation
                        #
                        # see github.com/dsirajud/IPython-Notebooks/DECSKS-04 for more information (final few sections of the notebook)
                        #
                        # Thus, we do not need to actually enforce the wall potential directly in terms of the charge accumulated for this boundary; however,
                        # we initialize and track the objects here so that the data can be accessed, analyzed or otherwise plotted, should the user wish
                        sigma[var]['upper'] = 0    # initialize to zero charge at time zero
                        sigma_n[var]['upper'] = np.zeros(Nt + 1)  # this was put in at one point for plotting wall charge vs. time
                    else:
                        print "\nThe following boundary conditions specified in params_boundaryconditions.dat are inconsistent together:\n\n"
                        print "upper boundary condition on phi: %s" % BC['phi'][var]['upper'].upper()
                        print "upper boundary condition on f: %s\n\n" % BC['f'][var]['upper'].upper()
                        print "\ne.g. an upper boundary condition set on phi as SELF-CONSISTENT must have the upper boundary condition on f as COLLECTOR\n "

                        raise InputError('An upper boundary condition on phi was specified as SELF-CONSISTENT; however, the corresponding boundary condition on f is not compatible (must be set to collector)')

                elif BC['phi'][var]['type'] == 'SELF-CONSISTENT_SELF-CONSISTENT':
                    BC['phi'][var]['type'] = 'LDBC_LNBC'

                    # We default to a LDBC_LNBC solver, both boundary conditions on left edge, entries 0 (Dirichlet) and 1 (Neumann)
                    # cf. DECSKS-04 notebook for more details:
                    #
                    #    https://github.com/dsirajud/IPython-notebooks/DECSKS-04...ipynb
                    #
                    # Dirichlet condition, set reference potential phi = 0
                    phi_BC[var][0] = 0. # reference potential set to zero
                    # Neumann condition, dphi = 0 for symmetry
                    #phi_BC[var][1] = - 6 * var.width * sim_params['sigma'][var]['lower'], changes with time step


                    if BC['f'][var]['lower'] == 'collector': # all synonyms for 'absorbing' (except 'collector') have been seen by this point, and if encountered changed to 'absorbing'
                        # initialize wall charge densities
                        sigma[var]['lower'] = 0    # initialize to zero charge at time zero
                        sigma_n[var]['lower'] = np.zeros(Nt + 1)  # this was put in at one point for plotting wall charge vs. time
                    else:
                        print "\nThe following boundary conditions specified in params_boundaryconditions.dat are inconsistent together:\n\n"
                        print "lower boundary condition on phi on variable %s: SELF-CONSISTENT" % var
                        print "lower boundary condition on f on variable %s: %s\n\n" % (var, BC['f'][var]['lower'].upper())
                        print "\ne.g. a lower boundary condition set on phi as SELF-CONSISTENT must have the lower boundary condition on f as COLLECTOR\n "

                        raise InputError('A lower boundary condition on phi was specified as SELF-CONSISTENT; however, the corresponding boundary condition on f is not compatible (must be set to collector)')

                    if BC['f'][var]['upper'] == 'collector': # all synonyms for 'absorbing' (except 'collector') have been seen by this point, and if encountered changed to 'absorbing'
                        # initialize wall charge densities, sigma for the collector (f) /self-consistent (phi) conditions
                        # By virtue of the setup, the above enforcements on the lower boundary ensures this unenforced upper Neumann BC is
                        # satisfied automatically given the relationship that Neumann BCs are fixed by due to the Poisson equation
                        #
                        # see github.com/dsirajud/IPython-Notebooks/DECSKS-04 for more information (final few sections of the notebook)
                        #
                        # Thus, we do not need to actually enforce the wall potential directly in terms of the charge accumulated for this boundary; however,
                        # we initialize and track the objects here so that the data can be accessed, analyzed or otherwise plotted, should the user wish
                        sigma[var]['upper'] = 0    # initialize to zero charge at time zero
                        sigma_n[var]['upper'] = np.zeros(Nt + 1)  # this was put in at one point for plotting wall charge vs. time
                    else:
                        print "\nThe following boundary conditions specified in params_boundaryconditions.dat are inconsistent together:\n\n"
                        print "upper boundary condition on phi: SELF-CONSISTENT"
                        print "upper boundary condition on f: %s\n\n" % BC['f'][var]['upper'].upper()
                        print "\ne.g an upper boundary condition set on phi as SELF-CONSISTENT must have the upper boundary condition on f as COLLECTOR\n "

                        raise InputError('An upper boundary condition on phi was specified as SELF-CONSISTENT; however, the corresponding boundary condition on f is not compatible (must be set to collector)')

              # else: boundary conditions have already been checked for valid inputs, no invalid input will be encountered

    # --------------------------------------------------------------------------
    # ELECTRIC POTENTIAL PHI FUNCTION HANDLE STRING and BOUNDARY CONDITION TYPE FUNCTION HANDLE STRING
    #
    # currently only 1D1V, only one handle needed. When this will be generalized, can make a dict objects with keys corresponding
    # to each active configuration variable
    #
    # The forms of each string call their associated method per the boundary conditions specified by the user in params_boundaryconditions.dat,
    # based on the boundary conditions specified by the user, one of the following will be created:
    #
    #    compute_electric_potential_phi_handle[var] =
    #
    #                     DECSKS.lib.fieldsolvers.Poisson_6th_PBC
    #                     DECSKS.lib.fieldsolvers.Poisson_6th_LDBC_UDBC
    #                     DECSKS.lib.fieldsolvers.Poisson_6th_LDBC_UNBC
    #                     DECSKS.lib.fieldsolvers.Poisson_6th_LNBC_UDBC
    #                     DECSKS.lib.fieldsolvers.Poisson_6th_LDBC_LNBC
    #                     DECSKS.lib.fieldsolvers.Poisson_6th_UDBC_UNBC (<-- available, but not used in any current combination of BCs)
    #
    #
    # and, one of the following
    #
    #    distribution_function_boundarycondition_handle[var]['lower'] =
    #
    #                     DECSKS.lib.boundaryconditions.absorbing_lower_boundary
    #                     DECSKS.lib.boundaryconditions.collector_lower_boundary
    #                     DECSKS.lib.boundaryconditions.symmetric_lower_boundary
    #
    #                     NOTE: if 'periodic' has been specified, everything is
    #                     handled in the orchestrator, distribution_function_boundarycondition_orchestrator
    #                     which would take on the string value = 'DECSKS.lib.boundaryconditions.periodic


    distribution_function_boundarycondition_prefix = 'DECSKS.lib.boundaryconditions'
    distribution_function_boundarycondition_handle = {}
    for var in phasespace_vars:
        if BC['f'][var]['type'] == 'periodic':
            pass
        else:
            distribution_function_boundarycondition_handle[var] = {}

            distribution_function_boundarycondition_handle[var]['lower'] = ".".join((distribution_function_boundarycondition_prefix,  BC['f'][var]['lower']))
            distribution_function_boundarycondition_handle[var]['lower'] = "_".join((distribution_function_boundarycondition_handle[var]['lower'],  'lower_boundary'))

            distribution_function_boundarycondition_handle[var]['upper'] = ".".join((distribution_function_boundarycondition_prefix,  BC['f'][var]['upper']))
            distribution_function_boundarycondition_handle[var]['upper'] = "_".join((distribution_function_boundarycondition_handle[var]['upper'],  'upper_boundary'))


    compute_electric_potential_phi_handle = {}
    compute_electric_potential_phi_prefix = "DECSKS.lib.fieldsolvers.Poisson_6th_"
    for var in ['x', 'y', 'z']:
        if var in phasespace_vars:
            if HOC[var] == 'FOURIER': # uses a Gauss law solver to find E directly, which is called by the orchestrator on the fieldsolver
                pass
            else: # computes the electric field E by differentiating phi in an orchestrator fieldsolver function (string handle constructed above)
                  # inside the orchestrator, a particular Poisson solver is called according with the boundary conditions indicated in params_boundaryconditions.dat
                compute_electric_potential_phi_handle[var] = compute_electric_potential_phi_prefix + BC['phi'][var]['type']
        else:
            pass

    # in the future, can generalize this to multiple dimensions by making this a dict with keys ['x', 'y', 'z']
    # currently just on 1D1V and expecting an 'x' variable to be evolved in configuration

    if 'x' not in phasespace_vars:
        raise NotImplementedError('Current 1D1V version of DECSKS is expecting x to be the active configuration variable. Please revise the intended simulation so that x is the symbol chosen in params.dat.')
    else:
        if HOC['x'] == 'FOURIER': # uses a Gauss solver to find E directly
            Poisson_6th_order_FD_solver_matrices = None

        else: # uses a Poisson solver to find phi, then differentiates to obtain E
            Poisson_6th_order_FD_solver_matrices = assemble_Poisson_6th_order_FD_solver_matrices(Nx, BC)

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
        total_dims = total_dims,
        density = density,
        mu = mu,
        split_scheme = split_scheme,
        splitting = splitting,
        plot_params = plot_params,
        record_outputs = record_outputs,
        outfiles = outfiles,
        BC = BC,    # boundary condition types on all phase space variables on distribution function f and phi
        phi_BC = phi_BC, # dictionary containing boundary value vector for electric potential used in Poisson solve, e.g. phi_BC['x']
        sigma = sigma,
        sigma_n = sigma_n, # this was put in for charge history plots
        distribution_function_boundarycondition_handle = distribution_function_boundarycondition_handle, # dictionary with keys (var in phasespace_vars), which are keys to a subdict with keys 'lower', 'upper'
        distribution_function_boundarycondition_orchestrator_handle = distribution_function_boundarycondition_orchestrator_handle, # dictionary with keys (var in phasespace_vars)
        compute_electric_potential_phi_handle = compute_electric_potential_phi_handle,
        compute_electric_field_orchestrator_handle = compute_electric_field_orchestrator_handle,
        I_alternating = I_alternating, # identity matrix with alternating signs according to row, used in computing correctors c
        A_matrix = A_matrix,     # Matrices of Bernoulli numbers for HOC
        W = W,
        W_dn1_LTE6 = W_dn1_LTE6,
        Xi = Xi, # spectral differentiation operator matrix (1j*xi[i,j]) ** q
        xi = xi, # wave number vector
        Poisson_6th_order_FD_solver_matrices = Poisson_6th_order_FD_solver_matrices
        )

    infile.close()

    # --------------------------------------------------------------------------
    # Before return, broadcast notification
    # regarding start of simulation and order of solver

    print "\nStarting 1D1V Vlasov-Poisson simulation"
    print "\nadvection solver: LTE order %d" % (N+1)
    print "\nwill step through %d-dimensional solution in variables: %s\n" % (len(phasespace_vars), phasespace_vars)
    for var in phasespace_vars:
        print "high order correction method on %s: %s" % (var, HOC[var])

    print "\n"
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
    splitting[ coeff[i] ][ int(stage[i]) ]
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

        number_of_stages = dict(a = 2, b = 2)
        order = dict(coeffs = coeffs, stages = stages)
        splitting = dict(order = order, number_of_stages = number_of_stages,
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

        number_of_stages = dict(a = 4, b = 4)
        order = dict(coeffs = coeffs, stages = stages)
        splitting = dict(order = order, number_of_stages = number_of_stages,
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

        number_of_stages = dict(a = 4, b = 4)
        order = dict(coeffs = coeffs, stages = stages)
        splitting = dict(order = order, number_of_stages = number_of_stages,
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

        number_of_stages = dict(a = 6, b = 6)
        order = dict(coeffs = coeffs, stages = stages)
        splitting = dict(order = order, number_of_stages = number_of_stages,
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

        number_of_stagess = dict(a = 8, b = 7)
        order = dict(coeffs = coeffs, stages = stages)
        splitting = dict(order = order, number_of_stages = number_of_stages,
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
    filename_S = lines[15][lines[15].find(':')+1:].strip()

    filepath_I1 = rel_path + filename_I1
    filepath_I2 = rel_path + filename_I2
    filepath_IW = rel_path + filename_IW
    filepath_WE = rel_path + filename_WE
    filepath_S = rel_path + filename_S

    outfile_I1 = open(filepath_I1, 'w')
    outfile_I2 = open(filepath_I2, 'w')
    outfile_IW = open(filepath_IW, 'w')
    outfile_WE = open(filepath_WE, 'w')
    outfile_S = open(filepath_S, 'w')

    outfiles = dict(I1 = outfile_I1,
                    I2 = outfile_I2,
                    IW = outfile_IW,
                    WE = outfile_WE,
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
                                          az,
                                          bz,
                                          Nz,
                                          Nvz,
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

    Nz -- (int) total number of gridpoints for z
    Nvz -- (int) total number of gridpoints for vz
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
    if Nvz is None:
        return None

    # domain widths
    Lz = float(bz - az) # used in wave number vector, xi
    zwidth = Lz / Nz # needed for 3D object, Xi = "(j zwidth xi)**q"

    # build wave vector xi for given z
    wave_index = np.arange(Nz)
    xi_z = np.where(wave_index <= Nz / 2,
              2*np.pi*wave_index / Lz,
              2*np.pi*(wave_index - Nz) / Lz)

    xi[z_str] = xi_z

    # Set up compound matrix Xi.
    # First, copy column vector xi along Nvz columns
    xi_2D = np.outer(xi_z, np.ones(Nvz))

    # set up vector extending in depth dimension so
    # broadcasting per ** operator produces the expected dims on Xi
    # i.e. Xi.shape = (N, z.N, vz.N)
    dn = np.arange(1,N).reshape(N-1,1,1)

    # with the previously formed objects with carefully chosen dims
    # we generate the required Xi object
    Xi[z_str] = (1j * zwidth * xi_2D) ** dn

    return Xi, xi

def assemble_Poisson_6th_order_FD_solver_matrices(Nx, BC):
    """
        forms the matrices D and B required for the 6th order finite
        difference based Poisson solver. The solvers come in several
        varieties (DBC = Dirichlet BC, NBC = Neumann BC, L = lower boundary
        U = upper boundary where L and U refer to lower and higher values
        of a configuration variable):

        Poisson_6th_PBC
        Poisson_6th_LDBC_UDBC
        Poisson_6th_LNBC_UDBC
        Poisson_6th_LDBC_UNBC

        For Neumann/Neumann conditions (is not a well-posed problem) we
        require recasting the NBC/NBC problem into an equivalent problem
        that is representative by the following Cauchy boundary condition
        setup:

        Poisson_6th_LDBC_LNBC
        Poisson_6th_UDBC_UNBC

        The matrices D and B for each method are slightly different. This routine
        determines the variety of solver needed as chosen from the list above,
        and assembles these matrices

        Inputs:
        Nx -- (int) this is Nx_active, the number of active grid points in a
              configuration variable x

        BC -- (dict) contains boundary condition information that has been
              determined based on user inputs in params_boundaryconditions.dat

        Outputs:

        D -- (ndarray, ndim = 2, shape = (Nx, Nx)) matrix of finite difference
             coefficients on phi

        B -- (ndarray, ndim = 2, shape = (Nx, Nx)) matrix of finite difference
        coefficients on the totaldensity n
    """

    Poisson_6th_order_FD_solver_matrices = {}

    # Nx is the number of active nodes in configuration
    if BC['phi']['x']['type'] == 'PBC':

        # assemble D, a matrix of difference coefficients on phi
        D = np.zeros([Nx,Nx])
        for i in range(Nx):
            if i == 0:
                D[i,-1] = 1
                D[i,i] = -2
                D[i,i+1] = 1
            elif i == Nx-1:
                D[i,0] = 1
                D[i,i] = -2
                D[i,i-1] = 1
            else:
                D[i,i-1] = 1
                D[i,i] = -2
                D[i,i+1] = 1

        # assemble B, a matrix of difference coefficients on the total density
        B = np.zeros([Nx,Nx])
        for i in range(Nx):
            if i == 0:
                B[i,-2] = -1/240.
                B[i,-1] = 1/10.
                B[i,i] = 97/120.
                B[i,i+1] = 1/10.
                B[i,i+2] = -1/240.

            if i == 1:
                B[i,-1] = -1/240.
                B[i,i-1] = 1/10.
                B[i,i] = 97/120.
                B[i,i+1] = 1/10.
                B[i,i+2] = -1/240.


            elif 1 < i < Nx-2:
                B[i,i-2] = -1/240.
                B[i,i-1] = 1/10.
                B[i,i] = 97/120.
                B[i,i+1] = 1/10.
                B[i,i+2] = -1/240.

            elif i == Nx-2:
                B[i,i-2] = -1/240.
                B[i,i-1] = 1/10.
                B[i,i] = 97/120.
                B[i,i+1] = 1/10.
                B[i,0] = -1/240.

            elif i == Nx-1:
                B[i,i-2] = -1/240.
                B[i,i-1] = 1/10.
                B[i,i] = 97/120.
                B[i,0] = 1/10.
                B[i,1] = -1/240.

    elif BC['phi']['x']['type'] == 'LDBC_UDBC':

        # assemble D, a matrix of difference coefficients on phi
        D = np.zeros([Nx,Nx])
        for i in range(Nx):
            if i == 0:
                D[i,i] = -3
                D[i,i+1] = 1
            elif i == Nx-1:
                D[i,i] = -3
                D[i,i-1] = 1
            else:
                D[i,i-1] = 1
                D[i,i] = -2
                D[i,i+1] = 1

        # assemble B, a matrix of difference coefficients on the total density
        B = np.zeros([Nx, Nx])
        for i in range(Nx):
            if i == 0:
                B[i,i] = 317/240.
                B[i,i+1] = -133/120.
                B[i,i+2] = 187/120.
                B[i,i+3] = -23/20.
                B[i,i+4] = 109/240.
                B[i,i+5] = -3/40.

            if i == 1:
                B[i,i-1] = 3/40.
                B[i,i] = 209/240.
                B[i,i+1] = 1/60.
                B[i,i+2] = 7/120.
                B[i,i+3] = -1/40.
                B[i,i+4] = 1/240.

            elif 1 < i < Nx-2:
                B[i,i-2] = -1/240.
                B[i,i-1] = 1/10.
                B[i,i] = 97/120.
                B[i,i+1] = 1/10.
                B[i,i+2] = -1/240.

            elif i == Nx-2:
                B[i,i-4] = 1/240.
                B[i,i-3] = -1/40.
                B[i,i-2] = 7/120.
                B[i,i-1] = 1/60.
                B[i,i] = 209/240.
                B[i,i+1] = 3/40.

            elif i == Nx-1:
                B[i,i] = 317/240.
                B[i,i-1] = -133/120.
                B[i,i-2] = 187/120.
                B[i,i-3] = -23/20.
                B[i,i-4] = 109/240.
                B[i,i-5] = -3/40.

    elif BC['phi']['x']['type'] == 'LNBC_UDBC':

        # assemble D, a matrix of difference coefficients on phi
        D = np.zeros((Nx,Nx))

        # LNBC row
        D[0,0] = -97/10.
        D[0,1] = 16.
        D[0,2] = -10
        D[0,3] = 5.
        D[0,4] = -3/2.
        D[0,5] = 1/5.

        # UDBC row
        D[-1,-1] = 1.

        # Poisson's equation rows
        for i in range(1,Nx-1):
            D[i,i-1] = 1
            D[i,i] = -2
            D[i,i+1] = 1

        # assemble B, a matrix of difference coefficients on the total density
        B = np.zeros((Nx,Nx))
        for i in range(B.shape[0]):
            if i == 0:
                B[i,i] = 317 / 240.
                B[i,i+1] = -133/120.
                B[i,i+2] = 187 / 120.
                B[i,i+3] = -23 / 20.
                B[i,i+4] = 109 / 240.
                B[i,i+5] = -3/40.

            elif i == 1:

                B[i, i-1] = 3 / 40.
                B[i, i] = 209 / 240.
                B[i,i+1] = 1 / 60.
                B[i,i+2] = 7 / 120.
                B[i,i+3] = -1 / 40.
                B[i,i+4] = 1 / 240.

            elif 2 <= i <= Nx-3:

                B[i,i-2] = -1/240.
                B[i,i-1] = 1/10.
                B[i,i] = 97/120.
                B[i,i+1] = 1/10.
                B[i,i+2] = -1/240.

            elif i == Nx-2:

                B[i,i+1] = 3 / 40.
                B[i,i] = 209 / 240.
                B[i,i-1] = 1 / 60.
                B[i,i-2] = 7 / 120.
                B[i,i-3] = -1 / 40.
                B[i,i-4] = 1 / 240.

            # else i == Nx-1: row of zeros

    elif BC['phi']['x']['type'] == 'LDBC_UNBC':

        # assemble D, a matrix of difference coefficients on phi
        D = np.zeros((Nx,Nx))

        # UDBC row
        D[0,0] = 1.

        # LNBC row
        D[-1,-1] = -97/10.
        D[-1,-2] = 16.
        D[-1,-3] = -10
        D[-1,-4] = 5.
        D[-1,-5] = -3/2.
        D[-1,-6] = 1/5.

        # Poisson's equation rows
        for i in range(1,Nx-1):
            D[i,i-1] = 1
            D[i,i] = -2
            D[i,i+1] = 1


        # assemble B, a matrix of difference coefficients on the total density
        B = np.zeros((Nx,Nx))
        for i in range(B.shape[0]):
            # i == 0 row contains all zeros

            if i == 1:

                B[i, i-1] = 3 / 40.
                B[i, i] = 209 / 240.
                B[i,i+1] = 1 / 60.
                B[i,i+2] = 7 / 120.
                B[i,i+3] = -1 / 40.
                B[i,i+4] = 1 / 240.

            elif 2 <= i <= Nx-3:

                B[i,i-2] = -1/240.
                B[i,i-1] = 1/10.
                B[i,i] = 97/120.
                B[i,i+1] = 1/10.
                B[i,i+2] = -1/240.

            elif i == Nx-2:

                B[i,i+1] = 3 / 40.
                B[i,i] = 209 / 240.
                B[i,i-1] = 1 / 60.
                B[i,i-2] = 7 / 120.
                B[i,i-3] = -1 / 40.
                B[i,i-4] = 1 / 240.

            if i == Nx-1:
                B[i,i-5] = -3/40.
                B[i,i-4] = 109 / 240.
                B[i,i-3] = -23 / 20.
                B[i,i-2] = 187 / 120.
                B[i,i-1] = -133/120.
                B[i,i] = 317 / 240.

    elif BC['phi']['x']['type'] == 'LDBC_LNBC':

        # assemble D, a matrix of difference coefficients on phi
        D = np.zeros((Nx,Nx))

        # LDBC row, (row 0)
        D[0,0] = 1.

        # LNBC row, (row 1)
        D[1,0] = -97/10.
        D[1,1] = 16.
        D[1,2] = -10
        D[1,3] = 5.
        D[1,4] = -3/2.
        D[1,5] = 1/5.

        # Poisson's equation rows
        for i in range(2,Nx):
            D[i,i-2] = 1
            D[i,i-1] = -2
            D[i,i] = 1

        # assemble B, a matrix of difference coefficients on the total density
        B = np.zeros((Nx,Nx))
        for i in range(1,B.shape[0]):
            # if i == 0: row of zeros, density is not involved (corresponds to DBC)

            if i == 1:
                B[i,i-1] = 317 / 240.
                B[i,i] = -133/120.
                B[i,i+1] = 187 / 120.
                B[i,i+2] = -23 / 20.
                B[i,i+3] = 109 / 240.
                B[i,i+4] = -3/40.

            if i == 2:
                B[i, i-2] = 3 / 40.
                B[i, i-1] = 209 / 240.
                B[i,i] = 1 / 60.
                B[i,i+1] = 7 / 120.
                B[i,i+2] = -1 / 40.
                B[i,i+3] = 1 / 240.

            elif 3 <= i <= Nx-2:
                B[i,i-3] = -1/240.
                B[i,i-2] = 1/10.
                B[i,i-1] = 97/120.
                B[i,i] = 1/10.
                B[i,i+1] = -1/240.

            elif i == Nx-1:
                B[i,i-5] = 1/240.
                B[i,i-4] = -1/40.
                B[i,i-3] = 7/120.
                B[i,i-2] = 1/60.
                B[i,i-1] = 209/240.
                B[i,i] = 3/40.

    elif BC['phi']['x']['type'] == 'UDBC_UNBC':

        # assemble D, a matrix of difference coefficients on phi
        D = np.zeros((Nx,Nx))

        # LDBC row, (row Nx-1)
        D[-1,-1] = 1.

        # LNBC row, (row Nx-2)
        D[-2,-1] = -97/10.
        D[-2,-2] = 16.
        D[-2,-3] = -10
        D[-2,-4] = 5.
        D[-2,-5] = -3/2.
        D[-2,-6] = 1/5.

        # Poisson's equation rows
        for i in range(Nx-2):
            D[i,i] = 1
            D[i,i+1] = -2
            D[i,i+2] = 1


        # assemble B, a matrix of difference coefficients on the total density
        B = np.zeros((Nx,Nx))
        for i in range(B.shape[0]):
            if i == 0:
                B[i,i] = 3/40.
                B[i,i+1] = 209/240.
                B[i,i+2] = 1/60.
                B[i,i+3] = 7/120.
                B[i,i+4] = -1/40.
                B[i,i+5] = 1/240.

            if 1 <= i < Nx-3:
                B[i,i-1] = -1/240.
                B[i,i] = 1/10.
                B[i,i+1] = 97/120.
                B[i,i+2] = 1/10.
                B[i,i+3] = -1/240.

            elif i == Nx-3:
                B[i,i-3] = 1/240.
                B[i,i-2] = -1/40.
                B[i,i-1] = 7/120.
                B[i,i] = 1/60.
                B[i,i+1] = 209/240.
                B[i,i+2] = 3/40.

            elif i == Nx-2:
                B[i,i+1] = 317 / 240.
                B[i,i] = -133/120.
                B[i,i-1] = 187 / 120.
                B[i,i-2] = -23 / 20.
                B[i,i-3] = 109 / 240.
                B[i,i-4] = -3/40.

            # else i == Nx - 1: row of zeros, density is not involved (corresponds to DBC)

    Poisson_6th_order_FD_solver_matrices['D'] = D
    Poisson_6th_order_FD_solver_matrices['B'] = B

    return Poisson_6th_order_FD_solver_matrices
