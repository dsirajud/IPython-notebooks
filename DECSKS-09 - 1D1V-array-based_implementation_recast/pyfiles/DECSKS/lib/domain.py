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
        if dim is not None and var.lower() != 't':
            # var = 'v', dim = 'x', 'y' or 'z'
            # instantiates vx, vy, or vz

            # set up number of total gridpoints (self.Ngridpoints)
            # vs. number of active gridpoints self.N (can be equal)
            self.Ngridpoints = sim_params['N' + var + dim]
            if sim_params['BC'].lower() == 'periodic':
                self.N = self.Ngridpoints - 1
            else:
                self.N = self.Ngridpoints

            self.a = float(sim_params['a' + var + dim])
            self.b = float(sim_params['b' + var + dim])
            self.L = float(self.b - self.a)
            self.width = self.L / (self.Ngridpoints - 1)
            self.str = var + dim

            self.prepoints = np.array(range(self.N))
            self.postpoints = np.zeros(self.N) # container, to be filled at each timestep

            numdims = len(sim_params['phasespace_vars'])
            if numdims == 1:
                self.MCs = np.zeros(self.N) # container, to be filled at each timestep

            elif numdims == 2:
                N1_str = 'N' + sim_params['phasespace_vars'][0]
                N2_str = 'N' + sim_params['phasespace_vars'][1]

                N1, N2 = sim_params[N1_str], sim_params[N2_str]
                if sim_params['BC'].lower() == 'periodic':
                    self.MCs = np.zeros([N1 - 1,N2 - 1])
                else:
                    self.MCs = np.zeros([N1,N2])

            # for convecting cells
            self.prepointvalues = np.zeros(self.N)
            self.prepointvalues = self.generate_Eulerian_mesh(self.N)

            # for plots
            self.gridpoints = np.array(range(self.Ngridpoints))
            self.gridvalues = self.generate_Eulerian_mesh(self.Ngridpoints)

        elif var.lower() != 't':
            # var = 'x','y', or 'z'; dim = None
            # instantiates x, y or z

            # set up number of total gridpoints (self.Ngridpoints)
            # vs. number of active gridpoints self.N (can be equal)
            self.Ngridpoints = sim_params['N' + var]
            if sim_params['BC'].lower() == 'periodic':
                self.N = self.Ngridpoints - 1
            else:
                self.N = self.Ngridpoints

            self.a = float(sim_params['a' + var])
            self.b = float(sim_params['b' + var])
            self.L = float(self.b - self.a)
            self.width = self.L / (self.Ngridpoints - 1)
            self.str = var

            self.prepoints = np.array(range(self.N))
            self.postpoints = np.zeros(self.N) # container, to be filled at each timestep

            numdims = len(sim_params['phasespace_vars'])
            if numdims == 1:
                self.MCs = np.zeros(self.N) # container, to be filled at each timestep

            elif numdims == 2:
                N1_str = 'N' + sim_params['phasespace_vars'][0]
                N2_str = 'N' + sim_params['phasespace_vars'][1]

                N1, N2 = sim_params[N1_str], sim_params[N2_str]
                if sim_params['BC'].lower() == 'periodic':
                    self.MCs = np.zeros([N1 - 1,N2 - 1])
                else:
                    self.MCs = np.zeros([N1,N2])

            # for convecting cells
            self.prepointvalues = np.zeros(self.N)
            self.prepointvalues = self.generate_Eulerian_mesh(self.N)

            # for plots
            self.gridpoints = np.array(range(self.Ngridpoints))
            self.gridvalues = self.generate_Eulerian_mesh(self.Ngridpoints)

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
        """Mesh generator for domain cells, w = {x,y,z,vx,vy,vz,t}"""
        w = np.zeros(Num)
        for i in range(Num):
            w[i] = self.a + i*self.width
        return w

    def generate_Lagrangian_mesh(self, z, vz, dt):
        """Advects all MCs one time step for a phase space variable z.

        inputs:
        z -- (ndarray, dim=1) phase space variable mesh
        vz -- (ndarray, dim=1) velocity vector whose entries pair with each MC
        dt -- (float) time step taken

        outputs:
        z_MC -- (ndarray, dim=1) array of raw positions of each MC after push
        """

        z_MC = z + vz*dt

        return z_MC
