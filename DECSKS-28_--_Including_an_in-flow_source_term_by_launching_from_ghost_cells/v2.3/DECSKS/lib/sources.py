import numpy as np
import DECSKS

def setup(x, vx, t, sim_params):

    # get source distributions
    Se = get_distribution(vx, source = 'maxwellian', vD = -1 / np.sqrt(sim_params['mu'])) # left-moving
    Si = get_distribution(vx, source = 'maxwellian', mu = sim_params['mu'], tau = 1.0, vD = -1 / np.sqrt(sim_params['mu'])) #left-moving

    # print out source selections
    print "--------------------------------------------------------------------------------"
    print "\nAn upper source term has been specified in params.dat. The following sources have been selected:"
    print "\nelectrons: Maxwellian with drift velocity vD = %s" % 'Bohm velocity'
    print "ions     : Maxwellian with drift velocity vD = %s, mu = %s, tau = %g\n" % ('Bohm velocity', 'Hydrogen', 1.0)
    print "--------------------------------------------------------------------------------"

    # create 1D data for relevant velocity values
    #    Se_data = np.ones(len(vx.gridvalues[vx.gridvalues < 0]))
    #    Si_data = np.ones(len(vx.gridvalues[vx.gridvalues < 0]))

    Se_data = Se(vx.gridvalues[vx.gridvalues < 0])
    Si_data = Si(vx.gridvalues[vx.gridvalues < 0])

    # init CFL numbers for ghost cells (GCFL) compute CFL numbers for all
    GCFL = SourceAdvection(vx, sim_params)
    GCFL.compute_numbers_for_all_stages(sim_params, x, vx, t)

    # compute remap array f_S and object N_S that tells how many source particles reach the lower boundary

    # electron source
    f_Se, N_Se  = determine_source_effect_on_grid(
        x, vx, Se_data, GCFL, sim_params
        )

    # ion source
    f_Si, N_Si  = determine_source_effect_on_grid(
        x, vx, Si_data, GCFL, sim_params
        )

    store_global_vars(sim_params, f_Se, f_Si, N_Se, N_Si)


#==============================================================================#
# supporting routines
#==============================================================================#

def store_global_vars(sim_params, f_Se, f_Si, N_Se, N_Si):

    sim_params['sources'] = {}
    sim_params['sources']['N_Se'] = N_Se
    sim_params['sources']['N_Si'] = N_Si

    sim_params['sources']['f_Se'] = f_Se
    sim_params['sources']['f_Si'] = f_Si

    return None

def get_distribution(vx, source = 'maxwellian', **kwargs):
    """
    Returns the function object for a specified source
    """

    def source_function(vx):

        if source == 'constant':

            return np.where(vx < 0)

        if source == 'maxwellian':

            # check inputs, provide normalization factors as compared to electron species, e.g. mu, tau
            if 'mu' in kwargs.keys():
                mu = kwargs['mu']
                assert ( type(mu) != str)
            else:
                mu = 1.0 # dummy value

            if 'tau' in kwargs.keys():
                 tau = kwargs['tau']
                 assert ( type(tau) != str)

            else:
                tau = 1.0 # dummy value

            if 'vD' in kwargs.keys():
                 vD = kwargs['vD']
                 assert ( type(vD) != str)

            else:
                print "\nCourtesy notice: no value for the drift velocity vD has been specified, using vD = 0.0"

            return 1 / np.sqrt(2*np.pi * tau / mu) * np.exp(-1/2. * (vx - vD) ** 2 / (tau / mu))

    return source_function


def determine_source_effect_on_grid(x,
                                    vx,
                                    S,
                                    GCFL,
                                    sim_params
                                    ):
    """
    create objects f_S and N_S that are the objects that capture
    the effect of a bulk source flowing in from the right

    inputs:
    x  -- (instance)
    vx -- (instance)
    S  -- (ndarray, ndim=1) source function data for all vx < 0, len(S) < vx.N
          only negative values
    GCFL -- (instance) ghost CFL numbers and parts for source terms
    sim_params -- (dict) simulation parameters

    outputs:
    f_S -- (ndarray, ndim=3), shape = (stages+1, x.N, vx.N) density array that is
           the contribution of the source particles to every [i,j] for
           a velocity dependent distribution function

           Note: while the dimensions span vx.N columns, it is usual that
           half of these entries would correspond to positive velocities
           and thus are empty (source particles with postive velocity cannot
           flow inward [left] to the on-grid domain). The array
           is dimensioned of equal size to the arrays fe and fi in DECSKS
           simulations so that straightforward addition can accomplish
           the remapping on top of the evolved distributions fe, fi

    N_S -- (ndarray, ndim = 1), shape = (stages+1,) the number of source
            particles that reach the lower boundary

    see notebook s28 for details
    """

    num_stages = sim_params['splitting']['number_of_stages']['a']

    # create an array is indexed as (s,i,j) where s labels the stage (>= 1)
    f_S = np.zeros((num_stages + 1, x.N, vx.N))

    # create a 1D array indexed by s, the stage (s >= 1; s = 0 entry is empty)
    N_S = np.zeros(num_stages + 1, dtype = int)

    # send to Cython routine
    f_S, N_S = DECSKS.lib.source_advection.create_postpoint_objects(S,
                                GCFL.int,
                                GCFL.frac,
                                num_stages,
                                x.N,
                                vx.N,
                                x.width,
                                vx.width,
                                sim_params
                                )

    return f_S, N_S

class SourceAdvection:

    def __init__(self, vx, sim_params):
        # incoming source from right-side, i.e. lefgoing velocity

        dim1 = (sim_params['splitting']['number_of_stages']['a'] + 1,) # one-tuple
        dim2 = (len(vx.gridvalues[vx.gridvalues < 0]),) # leftgoing only, one-tuple
        dims = dim1 + dim2

        self.numbers = np.zeros(dims)
        self.int     = np.zeros(dims, dtype = int)
        self.frac    = np.zeros(dims)

        # dims
        self.nrows    = dims[0] # access 0th entry of the one-tuple
        self.ncols    = dims[1] # access 1st entry of the one-tuple

        return None

    def compute_numbers_for_all_stages(self, sim_params, x, vx, t):
        """
        Essentially a copy of the routine from lib.domain
        the difference is the return is not 2D, but 1D with size
        corresponding to the velocity grid domain
        """

        splitting = sim_params['splitting']

        # recall stages start numbering at 1, s = 1, 2, ..., not zero
        # when we are accessing the information as stored in
        # the dict sim_params['splitting']
        for s in range(1, splitting['number_of_stages']['a'] + 1):
            split_coeff = splitting['a'][s]
            self.compute_numbers_for_stage_s(x, vx, split_coeff*t.width, s)

        return None

    def compute_numbers_for_stage_s(self, x, vx, dt, s):

        self.numbers[s,:] = vx.gridvalues[:self.ncols] * dt / x.width
        self.int[s,:]     = np.ceil(self.numbers[s,:]) # round towards zero
        self.frac[s,:]    = self.numbers[s,:] - self.int[s,:]

        return None
