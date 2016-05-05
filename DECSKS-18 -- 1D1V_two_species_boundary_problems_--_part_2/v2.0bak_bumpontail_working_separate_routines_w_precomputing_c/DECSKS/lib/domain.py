import numpy as np

class Setup:
    """Returns (t, x, y, z, vx, vy, vz) instances

    inputs:
    sim_params -- (dict) simulation parameters from read_input()
    var -- (str) v, {x,y,z}, or t
    dim -- (str) if var = 'v', dim = x, y, z must be supplied
                 else dim = None by keyword

    outputs:
    self -- (instance) time or phase space variable
    """
    def __init__(self, sim_params, var, dim = None):
        if var[0] == 'a':
            self.prepointvaluemesh = np.zeros(sim_params['active_dims'])
            self.Ngridpoints = sim_params['N' + dim] # Nx, Ny, or Nz
            self.str = var + dim # string label ax, ay, or az

            if sim_params['BC']['f'][dim]['type'] == 'periodic':
                # if the instance eval('a' + dim), e.g. ax, ay, az has a periodic boundary condition on the
                # underlying affiliated dimension x, y, or z
                self.N = self.Ngridpoints - 1
            else:
                self.N = self.Ngridpoints

            self.prepoints = np.arange(self.N)
            self.prepointmesh = np.array(np.outer( self.prepoints, np.ones([1, sim_params['active_dims'][1]]) ), dtype = int)

        elif dim is not None and var.lower() != 't':
            # var = 'v', dim = 'x', 'y' or 'z'
            # instantiates vx, vy, or vz

            # set up number of total gridpoints (self.Ngridpoints)
            # vs. number of active gridpoints self.N (can be equal)
            self.Ngridpoints = sim_params['N' + var + dim]
            if sim_params['BC']['f'][var + dim]['type'] == 'periodic':
                self.N = self.Ngridpoints - 1
            else:
                self.N = self.Ngridpoints

            self.a = float(sim_params['a' + var + dim])
            self.b = float(sim_params['b' + var + dim])
            self.L = float(self.b - self.a)
            self.width = self.L / (self.Ngridpoints - 1)
            self.str = var + dim

            self.prepoints = np.arange(self.N)
            self.prepointvalues = self.generate_Eulerian_mesh(self.N)

            postpointmesh_dims = [2]
            for dim in sim_params['active_dims']:
                postpointmesh_dims.append(dim)
            # container, to be filled at each timestep with postpoint index map
            self.postpointmesh = np.zeros(postpointmesh_dims, dtype = int)

            self.prepointvaluemesh = np.outer( np.ones([sim_params['active_dims'][0], 1]), self.prepointvalues)
            self.prepointmesh = np.array(np.outer( np.ones([sim_params['active_dims'][0], 1]), self.prepoints), dtype = int)

            # for plots
            self.gridpoints = np.array(range(self.Ngridpoints))
            self.gridvalues = self.generate_Eulerian_mesh(self.Ngridpoints)

            # sub-instantiation of CFL parameter class
            self.CFL = CourantNumberVelocity(self)

        elif var.lower() != 't':
            # var = 'x','y', or 'z'; dim = None
            # instantiates x, y or z

            # set up number of total gridpoints (self.Ngridpoints)
            # vs. number of active gridpoints self.N (can be equal)
            self.Ngridpoints = sim_params['N' + var]
            if sim_params['BC']['f'][var]['type'] == 'periodic':
                self.N = self.Ngridpoints - 1
            else:
                self.N = self.Ngridpoints

            self.a = float(sim_params['a' + var])
            self.b = float(sim_params['b' + var])
            self.L = float(self.b - self.a)
            self.width = self.L / (self.Ngridpoints - 1)
            self.str = var

            self.prepoints = np.arange(self.N)
            self.prepointvalues = self.generate_Eulerian_mesh(self.N)

            postpointmesh_dims = [2]
            for dim in sim_params['active_dims']:
                postpointmesh_dims.append(dim)
            # container, to be filled at each timestep with postpoint index map for each [i,j]
            self.postpointmesh = np.zeros(postpointmesh_dims, dtype = int)

            self.prepointvaluemesh = np.outer( self.prepointvalues, np.ones([1, sim_params['active_dims'][1]] ) )
            self.prepointmesh = np.array(np.outer( self.prepoints, np.ones([1, sim_params['active_dims'][1]]) ), dtype = int)

            # for plots
            self.gridpoints = np.array(range(self.Ngridpoints))
            self.gridvalues = self.generate_Eulerian_mesh(self.Ngridpoints)

            # sub-instantiation of CFL parameter class
            self.CFL = CourantNumberConfiguration(self, sim_params)

        else:
            # var = 't', dim = None
            # instantiates t

            self.N = sim_params['N' + var] # number of timesteps
            self.Ngridpoints = self.N + 1 # total number of time grid points

            # the following are defined for us in the method self.generate_Eulerian_mesh
            self.a = 0.0
            self.b = float(sim_params['T'])

            self.T = float(self.b - self.a)
            self.width = self.T / self.N
            self.stepnumbers = np.array(range(1,self.Ngridpoints))
            self.times = self.generate_Eulerian_mesh(self.Ngridpoints)
            self.str = var

    def generate_Eulerian_mesh(self, Num):
        """Mesh generator for domain cells, w = {x,y,z,vx,vy,vz,t}

        Num is included in general since t contains t.N + 1 cells
        whereas x, ..., vx, ... contain x.N, y.N, etc."""
        w = np.zeros(Num)
        for i in range(Num):
            w[i] = self.a + i*self.width
        return w

class CourantNumberConfiguration:
    """Creates an instance affiliated with a configuration variable z
    that describes the advection (CFL numbers)
    """
    def __init__(self, z, sim_params):
        """
        Inits the three 3D arrays characterizing the CFL numbers for a phase
        space variable z whose dimensions are (S+1,) + z.prepointmesh.shape
        here, S is the number of stages. The zeroeth stage is not used as
        we start numbering from 1, i.e. stages s = 1, 2, 3, ... S

            z.CFL.numbers -- (ndarray, ndim=3) raw CFL numbers at each [i,j]
                             where i labels one of [x,y,z] and j labels the
                             respective velocity in [vx, vy, vz] for
                             stage s of a split scheme, e.g. if the
                             configuration variable is chosen (arbitrarily)
                             as the coefficient stored as 'a', and the split
                             scheme has the interleaving a1, b1, a2, b2, a3, ...
                             the stage refers to 1 when a1 is applied,
                             then 2 when a2 is applied, then 3 when a3 is applied,
                             and so on .... Note that the other coefficient ('b')
                             would be associated with the velocity advection by
                             our already choosing 'a' to represent a configuration
                             step.

            z.CFL.int -- (ndarray, ndim=3) integer part of the above CFL numbers
            z.CFL.frac -- (ndarray, ndim=3) fractional part of the above CFL numbers
        """
        # stage 0 has no meaning, we label stages s = 1, 2, 3, ...
        dim1 = (sim_params['splitting']['number_of_stages']['a'] + 1,)
        dims = dim1 + z.prepointmesh.shape

        self.numbers = np.zeros(dims)
        self.int = np.zeros(dims, dtype = int)
        self.frac = np.zeros(dims)

        return None

    def compute_numbers_for_stage_s(self, z, vz, dt, s):
        """
        Returns a 2D array that is stored in the first dimension at index
        s of the above listed 3D arrays,

        Inputs:
        self -- (instance) z.CFL
        z --(instance) configuration variable z
        vz -- (instance) physical velocity vz for the configuration variable z
        dt -- (float) split_coeff * t.width, the time step at this stage
        s -- (int) stage number: 1, 2, 3, ... (does not start at 0)

        updated attributes:
        z.CFL.numbers[s,:,:] -- (ndarray, ndim=2) raw CFL numbers at each [i,j]
        z.CFL.int[s,:,:] -- (ndarray, ndim=2) integer part of CFL numbers at each [i,j]
        z.CFL.frac[s,:,:] -- (ndarray,  ndim=2)  fractional part of CFL numbers at each [i,j]

        Outputs:
        None -- updates attriutes listed above, no return needed
        """


        self.numbers[s,:,:] = vz.prepointvaluemesh * dt / z.width

        # if >= 0 , self.int = np.floor(self.numbers), else np.ceil(self.numbers)
        # i.e. the sign of CFL.frac agrees with the sign of vz.prepointvaluemesh
        self.int[s,:,:] = np.where(self.numbers[s,:,:] < 0, np.ceil(self.numbers[s,:,:]),
                                   np.floor(self.numbers[s,:,:])).astype(int)

        # remaining portion is the fractional CFL number
        self.frac[s,:,:] = self.numbers[s,:,:] - self.int[s,:,:]

        return None

    def compute_numbers_for_all_stages(self, sim_params, z, vz, t):
        """
        NOTICE: we choose (arbitrarily) that 'a' will be a configuration variable,
        hence we use the split coefficients a1, a2, a3, ... for the stored
        split scheme whose coefficients are stored in the dict object
        sim_params['splitting'] which is a subdict object with a key 'a'
        that contains a list of the coefficients a1, a2, a3, i.e.

        sim_params['splitting'] is a subdict that has the key:value pair

            splitting['a'] : [None, a1, a2, a3, ...]

        where a1, a2, a3 are floats. Thus, we access a1, a2, a3 as

            sim_params['splitting']['a'][1], sim_params['splitting']['a'][2],
                sim_params['splitting']['a'][3], ...

        note that since this is an attribute method of the subinstance z.CFL,
        and python (though users may disagree with the specific language used here)
        passes attributes of instances by reference, the changes are already
        applied, hence we return None at the conclusion of this method.
        """

        splitting = sim_params['splitting']

        # recall stages start numbering at 1, s = 1, 2, ..., not zero
        # when we are accessing the information as stored in
        # the dict sim_params['splitting']
        for s in range(1, splitting['number_of_stages']['a'] + 1):
            split_coeff = splitting['a'][s]
            self.compute_numbers_for_stage_s(z, vz, split_coeff*t.width, s)

        return None


class CourantNumberVelocity:
    """Returns a CFL number instance of the phasespace variable z

    inputs:
    z -- (instance) phase space variable from class lib.domain.Setup

    outputs:
    self -- (instance) CFL number ascribed to variable z convection

    Note: broadcasting ensures self.numbers is the same shape as z.MCs
    """
    def __init__(self, z):

        self.numbers = np.zeros(z.prepointmesh.shape)
        self.frac = np.zeros(z.prepointmesh.shape)
        self.int = np.zeros(z.prepointmesh.shape)

    def compute_numbers(self, z, vz, dt):
        """Calculates the CFL numbers and corresponding integer and fractional
        parts for each [i,j] of a velocity instance z

        In the advection step inside lib.convect_configuration.advection_step,
        We implement the indicial displacement of each gridpoint according
        to the velocity values in each column by shifting from prepoints

            (i,j) --> (i,j+ CFL.numbers[j])

        where the index i is not explicitly referenced given it is obvious.

        inputs:
        self -- (lib.domain.CourantNumber instance) CFL instance with attribute containers
                containers CFL.numbers, CFL.int, CFL.frac.

                NOTE: the method is an attribute belonging to the subinstance z.CFL
                hence, self refers to z.CFL, not z.

        z -- (lib.domain.Setup instance) phasespace instance being advected
        vz -- (lib.domain.Setup instance) velocity for self, e.g. vx, ..., ax, ..
        dt -- (float) width of time step, can be a fraction of t.width for split
              methods

        outputs:
        None -- updates attributes
        """
        self.numbers = vz.prepointvaluemesh * dt / z.width

        # if >= 0 , self.int = floor(self.numbers), else ceil(self.numbers)
        # i.e. the sign of CFL.frac agrees with the sign of vz
        self.int = np.where(self.numbers < 0, np.ceil(self.numbers),
                            np.floor(self.numbers)).astype(int)

        # remaining portion is the fractional CFL number
        self.frac = self.numbers - self.int

        # format dtype as int
        self.int = np.array(self.int, dtype = int)

        return None

def velocity_advection_prep(f_initial, z, vz):
    """
    When advecting physical velocity variables, the implementation
    requires several transpositions. This method performs those
    operations. Note, the instances z and vz will are changed
    by reference, the ndarrays f_final and f_initial need to be
    referred (assigned) back to the caller.

    inputs:
    f_initial -- (ndarray, ndim=2) shape = (x.N, vx.N)
    z -- (instance) phase space variable being evolved

        z.prepointmesh -- (ndarray, ndim=2)
        z.postpointmesh -- (ndarray, ndim=3), shape = (2, x.N, vx.N)

        z.CFL.frac -- (ndarray, ndim=2) shape = (x.N, vx.N)
        z.CFL.int -- (ndarray, ndim=2) shape = (x.N, vx.N)

    vz -- (instance) generalized velocity, here vz = acceleration since z = vel.

        z.prepointmesh -- (ndarray, ndim=2) shape = (x.N, vx.N)

    outputs:

    f_initial -- (ndarray, ndim=2) shape = (vx.N, x.N)
    z -- (instance) phase space variable being evolved

        z.prepointmesh -- (ndarray, ndim=2) shape = (vx.N, x.N)
        z.postpointmesh -- (ndarray, ndim=3), shape = (2, vx.N, x.N)

        z.CFL.frac -- (ndarray, ndim=2) shape = (vx.N, x.N)
        z.CFL.int -- (ndarray, ndim=2) shape = (vx.N, x.N)

    vz -- (instance) generalized velocity, here vz = acceleration since z = vel.

        z.prepointmesh -- (ndarray, ndim=2) shape = (vx.N, x.N)

    NOTE: technically, we do not need to return the instances z and vz, as the
    changes persist outside the function scope. We do so for clarity.

    """

    f_initial = np.transpose(f_initial)

    z.prepointmesh = np.transpose(z.prepointmesh)
    z.postpointmesh = np.transpose(z.postpointmesh, (0,2,1))

    z.CFL.frac = np.transpose(z.CFL.frac)
    z.CFL.int = np.transpose(z.CFL.int)

    vz.prepointmesh = np.transpose(vz.prepointmesh)

    return f_initial, z, vz

def velocity_advection_postproc(f_remapped, z, vz):
    """
    This function undoes the transpositions of velocity_advection_prep.
    It is written as a new method just for clarity. It is not passed
    to the other function to not include needless function call overhead.

    By this point, the initial density f_initial (shape = (vx.N, x.N))
    has been remapped (after all advection) to f_remapped (shape = (vx.N, x.N)).
    This function returns f_remapped (shape = (x.N, vx.N)) so that it can be
    stuffed into the final density container f_final (shape = (x.N, vx.N))

    It also transposes the previously transposed ndarray attributes of
    the instances z and vz.

    inputs:
    f_remapped -- (ndarray, ndim=2) shape = (vx.N, x.N)
    z -- (instance) phase space variable being evolved

        z.prepointmesh -- (ndarray, ndim=2)
        z.postpointmesh -- (ndarray, ndim=3), shape = (2, vx.N, x.N)

        z.CFL.frac -- (ndarray, ndim=2) shape = (vx.N, x.N)
        z.CFL.int -- (ndarray, ndim=2) shape = (vx.N, x.N)

    vz -- (instance) generalized velocity, here vz = acceleration since z = vel.

        z.prepointmesh -- (ndarray, ndim=2) shape = (vx.N, x.N)

    outputs:

    f_remapped -- (ndarray, ndim=2) shape = (x.N, vx.N)
    z -- (instance) phase space variable being evolved

        z.prepointmesh -- (ndarray, ndim=2) shape = (x.N, vx.N)
        z.postpointmesh -- (ndarray, ndim=3), shape = (2, x.N, vx.N)

        z.CFL.frac -- (ndarray, ndim=2) shape = (x.N, vx.N)
        z.CFL.int -- (ndarray, ndim=2) shape = (x.N, vx.N)

    vz -- (instance) generalized velocity, here vz = acceleration since z = vel.

        z.prepointmesh -- (ndarray, ndim=2) shape = (x.N, vx.N)

    NOTE: technically, we do not need to return the instances z and vz, as the
    changes persist outside the function scope. We do so for clarity.

    """

    f_remapped = np.transpose(f_remapped)

    z.prepointmesh = np.transpose(z.prepointmesh)
    z.postpointmesh = np.transpose(z.postpointmesh, (0,2,1))

    z.CFL.frac = np.transpose(z.CFL.frac)
    z.CFL.int = np.transpose(z.CFL.int)

    vz.prepointmesh = np.transpose(vz.prepointmesh)

    return f_remapped, z, vz

def extract_active_grid(f_total_grid, z, sim_params):
    """We evolve the density from the previous time step, f_old
    only on the gridpoints that are 'active' (cf. DECSKS-09 notebook)
    We distinguish, e.g. in 1D, the two attributes of a phase space
    variable z:

        z.Ngridpoints -- (int) total number of gridpoints
        z.N           -- (int) total number of 'active' gridpoints

    The total grid indices  : 0, 1, ... , z.Ngridpoints - 1
    The active grid indices : 0, 1, ... , z.N - 1

    For all but periodic boundary conditions (PBCs), these are the same.
    That is, for periodic boundary conditions (PBCs):

        z.N = z.Ngridpoints - 1

    so we evolve f_old[:z.N] -> f_new[:z.N]

    and we complete by the density by periodicity:

        f_new[z.Ngridpoints - 1] = f_new[0]

    for all other BCs: z.N = z.Ngridpoints and this function has no
    effect.

    inputs:
    f_total_grid -- (ndarray, ndim=2) 2D density constituting total grid

    outputs:
    f_active_grid -- (ndarray, ndim=2) 2D density containing all
                     active gridpoints

    """

    f_active_grid = f_total_grid[0:sim_params['active_dims'][0], 0:sim_params['active_dims'][1]]

    return f_active_grid
