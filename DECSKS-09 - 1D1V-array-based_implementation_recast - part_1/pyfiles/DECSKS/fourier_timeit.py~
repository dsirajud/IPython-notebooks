import _mypath
import DECSKS

def setup():

    class phasespace_var:
        def __init__(self, Nz, az, bz, name = 'x'):
            self.N = Nz
            self.az = az
            self.bz = bz
            self.L = float(bz - az)
            self.width = self.L / (self.N - 1)
            self.gridvalues = np.linspace(az, bz, Nz)

            self.str = name

    import numpy as np

    # constants
    epsi = 0.04
    k = 0.3
    vT1 = 2.
    mu1 = 0.
    vT2 = 1/4.
    mu2 = 4.5
    A = 1 / np.sqrt(2*np.pi)
    B = 0.9
    C = 0.2

    N = 21

    Nx, Nv = 768, 1536
    ax, bx = -2*np.pi / k, 2*np.pi / k
    av, bv = -8., 8.

    x = phasespace_var(Nx, ax, bx, name = 'x')
    vx = phasespace_var(Nv, av, bv, name = 'v')
    X,V = np.meshgrid(x.gridvalues,vx.gridvalues)

    f = np.zeros([x.N,vx.N])
    f = A*(1 + epsi*np.cos(k*X))*(B*np.exp( -(V - mu1) ** 2 / vT1) + C*np.exp( -(V - mu2) ** 2 / vT2))

    f = f.T

    sim_params = DECSKS.lib.read.inputfile('./etc/params.dat')

    if sim_params['HOC'] == 'FD':
        sim_params['W'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrices(sim_params,x,vx)
        sim_params['W_dn1'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrix_single_derivative(sim_params,x,dn = 1, LTE = 6)

    return f, sim_params

def main(f, sim_params):
    d = sim_params['W']['x'].dot(f)
    return d

    # timeit setup

    # import timeit
    #stmt = "main(f, sim_params)"
    #setup = "from derivative_timeit import setup, main; import _mypath; import DECSKS; f, sim_params = setup()"
    #print(timeit.timeit(stmt, setup, number = 1))



