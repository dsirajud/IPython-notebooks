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
            self.v_str = 'a' + dim

            self.prepoints = np.array(range(self.N))
            self.prepointvalues = self.generate_Eulerian_mesh(self.N)

            if sim_params['numdims'] == 1:
                self.MCs = np.zeros(self.N) # container, to be filled at each timestep
                self.prepointvaluemesh = self.prepointvalues
                self.prepointmesh = self.prepoints
                self.postpointmesh = np.zeros(self.N) # container, to be filled at each timestep

            elif sim_params['numdims'] == 2:
                N1_str = 'N' + sim_params['phasespace_vars'][0] # Nx, Ny, or Nz
                N2_str = 'N' + sim_params['phasespace_vars'][1] # Nvx, Nvy, or Nvz

                N1, N2 = sim_params[N1_str], sim_params[N2_str]
                if sim_params['BC'].lower() == 'periodic':
                    N1 -= 1
                    N2 -= 1
                self.MCs = np.zeros([N1,N2])
                self.postpointmesh = np.zeros([N1, N2]) # container, to be filled at each timestep

                # prepointvaluemesh is of dimensions (N1, N2), self = x, y, or z, self.N = N1
                self.prepointvaluemesh = np.outer( np.ones([N1, 1]), self.prepointvalues)
                self.prepointmesh = np.outer( np.ones([N1, 1]), self.prepoints)

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
            self.v_str = 'v' + var

            self.prepoints = np.array(range(self.N))
            self.prepointvalues = self.generate_Eulerian_mesh(self.N)

            if sim_params['numdims'] == 1:
                self.MCs = np.zeros(self.N) # container, to be filled at each timestep
                self.prepointvaluemesh = self.prepointvalues
                self.prepointmesh = self.prepoints
                self.postpointmesh = np.zeros(self.N) # container, to be filled at each timestep

            elif sim_params['numdims'] == 2:
                N1_str = 'N' + sim_params['phasespace_vars'][0] # Nx, Ny, or Nz on params.dat
                N2_str = 'N' + sim_params['phasespace_vars'][1] # Nvx, Nvy, or Nvz on params.dat

                N1, N2 = sim_params[N1_str], sim_params[N2_str]
                if sim_params['BC'].lower() == 'periodic':
                    N1 -= 1
                    N2 -= 1
                self.MCs = np.zeros([N1,N2]) # container, to be filled at each timestep

                # prepointvaluemesh is of dimensions (N1, N2), self = x, y, or z, self.N = N1
                self.prepointvaluemesh = np.outer( self.prepointvalues, np.ones([1, N2]))
                self.prepointmesh = np.outer( self.prepoints, np.ones([1, N2]) )
                self.postpointmesh = np.zeros([N1, N2]) # container, to be filled at each timestep

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

    def generate_Lagrangian_mesh(self, v, dt):
        """Advects all MCs one time step for a phase space variable z (self)
        according for all velocities vz in the grid and timestep dt.

        inputs:
        self -- (instance) phase space variable being convected
        v -- (instance or ndarray[, ndim = 1]) velocity instance
               or vector whose entries pair with each MC
        dt -- (float) time step taken

        generates:
        self.MCs --

                if 1D
                (ndarray, ndim=1) 1D array of raw positions of each MC after
                acvection by v*dt for v = constant

                if 2D
                (ndarray, dim=2) 2D array of raw positions of each MC after
                advection by vz*dt for all vz in mesh.

                for 1D1V, the shape is (x.N, v.N) always

                if self.str = 'x', then z_MC[:,j] gives the postpoints of
                z_MCs after advection by vz[j]*dt

                if self.str = 'vx', then z_MC[i,:] gives the postpoints of
                z_MCS after advection by a[i]*dt

        outputs:
        None -- the purpose of this function is to store values in the attribute
                self.MCs
        """

        self.MCs = self.prepointvaluemesh + v*dt

        return None
