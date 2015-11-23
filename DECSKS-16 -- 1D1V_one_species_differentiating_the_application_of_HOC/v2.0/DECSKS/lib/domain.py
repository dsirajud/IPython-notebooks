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

            if ((sim_params['BC'][dim]['lower'] == 'periodic') and sim_params['BC'][dim]['upper'] == 'periodic'):
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
            if (sim_params['BC'][var + dim]['lower'] == 'periodic' and sim_params['BC'][var + dim]['upper']) == 'periodic':
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
            self.CFL = CourantNumber(self)

        elif var.lower() != 't':
            # var = 'x','y', or 'z'; dim = None
            # instantiates x, y or z

            # set up number of total gridpoints (self.Ngridpoints)
            # vs. number of active gridpoints self.N (can be equal)
            self.Ngridpoints = sim_params['N' + var]
            if (sim_params['BC'][var]['lower'] == 'periodic' and sim_params['BC'][var]['upper']) == 'periodic':
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
            self.CFL = CourantNumber(self)

        else:
            # var = 't', dim = None
            # instantiates t

            self.N = sim_params['N' + var] # number of timesteps
            self.Ngridpoints = self.N + 1 # total number of time grid points

            # the following are defined for us in the method self.generate_Eulerian_mesh
            self.a = 0.0
            self.b = float(sim_params['T'])

            self.T = float(sim_params['T'])
            self.width = self.T / self.N
            self.stepnumbers = np.arange(1,self.Ngridpoints)
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

class CourantNumber:
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
        parts for each col of z.prepointmesh and stores in the 2D stack

            z.CFL.compute_numbers(z,vz,dt)

        note that each number corresponds to advection in 1D for each 1D problem
        whose subdomain is the column, and whose column space constitutes the entire
        grid.

        Hence, we implement the indicial displacement of each gridpoint according
        to the velocity values in each column by shifting from prepoints

            (i,j) --> (i,j+ CFL.numbers[j])

        where the index i is not explicitly referenced given it is obvious.n

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
        self.int = np.where(self.numbers >=0, np.floor(self.numbers),
                            np.ceil(self.numbers)).astype(int)

        # remaining portion is the fractional CFL number
        self.frac = self.numbers - self.int

        # format dtype as int
        self.int = np.array(self.int, dtype = int)


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
