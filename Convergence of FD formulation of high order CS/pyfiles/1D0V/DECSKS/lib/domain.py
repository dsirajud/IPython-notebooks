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
    def __init__(self, sim_params, var):
        if var.lower() != 't':
            # var = 'x','y', or 'z'; dim = None
            self.N = sim_params['N' + var]
            self.a = float(sim_params['a' + var])
            self.b = float(sim_params['b' + var])
            self.L = float(self.b - self.a)
            self.width = self.L / self.N
            self.str = var
            # self.gridpoints and self.prepoints are
            # identical, but used in different loops in
            # the implementation depending on which
            # label is most suitable for a given context
            self.gridpoints = np.array(range(self.N))
            self.prepoints = np.array(range(self.N))
            self.postpoints = np.zeros(self.N)

            self.cells = self.generate_Eulerian_mesh(self.N)
            self.MCs = np.zeros(self.N)

        else:
            # var = 't', dim = None
            self.N = sim_params['N' + var]
            self.a = float(sim_params['a' + var])
            self.b = float(sim_params['b' + var])
            self.T = float(self.b - self.a)
            self.width = self.T / self.N
            self.stepnumbers = np.array(range(1,self.N+1))
            self.times = self.generate_Eulerian_mesh(self.N+1)
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
        vvz -- (ndarray, dim=1) velocity vector whose entries pair with each MC
        dt -- (float) time step taken

        outputs:
        z_MC -- (ndarray, dim=1) array of raw positions of each MC after push
        """

        z_MC = z + vz*dt

        return z_MC

def Lagrangian_mesh(z, vz, dt):
    """Advects all MCs one time step for a phase space variable z.

    inputs:
    z -- (ndarray, dim=1) phase space variable mesh
    vvz -- (ndarray, dim=1) velocity vector whose entries pair with each MC
    dt -- (float) time step taken

    outputs:
    z_MC -- (ndarray, dim=1) array of raw positions of each MC after push
    """

    z_MC = z + vz*dt

    return z_MC
