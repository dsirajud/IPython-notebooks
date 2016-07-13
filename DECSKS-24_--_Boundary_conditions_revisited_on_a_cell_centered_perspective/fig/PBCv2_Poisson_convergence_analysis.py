import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as LA
import Poisson_solvers

if __name__ == '__main__':

    Nx = [12, 24, 48, 96, 192, 384, 768, 1536]
    num_grids = len(Nx)

    error_norm = np.zeros(num_grids)
    orders = np.zeros(num_grids)

    plt.figure()

    for grid in range(num_grids):

        # grid dependent parameters
        a,b = 0, 2
        L = float(b - a)
        dx = L / Nx[grid]

        x = np.zeros(Nx[grid])
        x[0] = a + dx / 2.
        for i in range(Nx[grid]):
            x[i] = x[0] + i*dx

        n = np.sin(np.pi * x)

        phi = Poisson_solvers.PBCv2(n, ax = 0., bx = 2., Nx = Nx[grid], dx = dx, LBC = None, RBC = None)
        phi_exact = - 1 / np.pi ** 2 * np.sin(np.pi * x)

        error_norm[grid] = LA.norm(phi_exact - phi,2) * np.sqrt(dx / L)

        if grid == 0:
            print "Nx%d        error = %g       ----" % (Nx[grid], error_norm[grid])
        else:
            orders[grid] = np.log2(error_norm[grid-1] / error_norm[grid])
            print "Nx%d        error = %g       order = %g" % (Nx[grid], error_norm[grid], orders[grid])

        plt.plot(x,phi, 'o', lw = 2, label = 'Nx = %d' % Nx[grid])
        plt.plot(x,phi_exact, '--', lw = 2, label = 'exact')
        plt.grid()
        plt.legend(loc = 'upper left')
        grid_str = 'PBC_Nx_' + '%d' % Nx[grid] + '.png'
        plt.savefig(grid_str)
        plt.clf()

    print '\n'