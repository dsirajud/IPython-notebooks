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
                # if the instance eval('a' + dim), e.g. ax, ay, az has a periodic boundary condition on the dimension
                self.N = self.Ngridpoints - 1
            else:
                self.N = self.Ngridpoints

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
            self.postpointmesh = np.zeros(sim_params['active_dims']) # container, to be filled at each timestep

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
            self.postpointmesh = np.zeros(sim_params['active_dims']) # container, to be filled at each timestep

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
            self.a = float(sim_params['a' + var])
            self.b = float(sim_params['b' + var])
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

class CourantNumber:
    """Returns a CFL number instance of the phasespace variable z

    inputs:
    z -- (instance) phase space variable from class lib.domain.Setup

    outputs:
    self -- (instance) CFL number ascribed to variable z convection

    Note: broadcasting ensures self.numbers is the same shape as z.MCs
    """
    def __init__(self, z):

        self.numbers = np.zeros_like(z.prepointmesh)
        self.frac = np.zeros_like(z.prepointmesh)
        self.int = np.zeros_like(z.prepointmesh)


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
                            np.ceil(self.numbers))

        # remaining portion is the fractional CFL number
        self.frac = self.numbers - self.int

        # format dtype as int
        self.int = np.array(self.int, dtype = int)
