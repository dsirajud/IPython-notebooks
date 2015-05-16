import _mypath
import pyfiles
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as LA

def function(_x, a = 3/4., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_c = -0.25):

    return a*np.exp(-((_x + x_a) / w_a) ** 2) + np.exp(-((_x/w_b)**2)) + c*np.exp(-((_x + x_c)/w_c) ** 2)

def df1(_x, a = 3/4., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_c = -0.25):

    return a*(-2./w_a**2)*(_x + x_a)*np.exp(-((_x + x_a) / w_a)**2) \
      + (-2./w_b**2)*_x*np.exp(-(_x/w_b)**2) \
      + c*(-2./w_c**2)*(_x + x_c)*np.exp(-((_x + x_c)/w_c)**2)

def df2(_x, a = 3/4., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_c = -0.25):

    return a*(-2./w_a**2) * (1 - 2./w_a**2 * (_x + x_a)**2) * np.exp(-( (_x + x_a) / w_a)**2) \
      + (-2./w_b**2) * (1 - 2./w_b**2 * _x ** 2) * np.exp(-(_x / w_b)**2) \
      + c*(-2./w_c**2) * (1 - 2/w_c**2 * (_x + x_c)**2) * np.exp(-( (_x + x_c) / w_c) ** 2)

def df3(_x, a = 3/4., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_c = -0.25):

    return a * (-2. / w_a**2) ** 2 * ( 3*(_x + x_a) - 2./w_a**2 * (_x + x_a) ** 3) * np.exp(-( (_x + x_a) / w_a) ** 2) \
      + (-2./w_b**2)**2 * (3*_x - 2./w_b**2 * _x**3) * np.exp(-(_x / w_b)**2) \
      + c*(-2./w_c**2) ** 2 * (3*(_x + x_c) - 2./w_c**2 * (_x + x_c)**3) * np.exp(-( (_x + x_c) / w_c) ** 2)

def df4(_x, a = 3/4., b = 1., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_b = 0., x_c = -0.25):

    return a * (-2./w_a ** 2) ** 2 * ( 3*(1 + 2*(-2./w_a**2)*(_x + x_a) ** 2) + (-2. / w_a**2)**2 * (_x + x_a) ** 4) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
            + b * (-2./w_b ** 2) ** 2 * ( 3*(1 + 2*(-2./w_b**2)*(_x + x_b) ** 2) + (-2. / w_b**2)**2 * (_x + x_b) ** 4) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
            + c * (-2./w_c ** 2) ** 2 * ( 3*(1 + 2*(-2./w_c**2)*(_x + x_c) ** 2) + (-2. / w_c**2)**2 * (_x + x_c) ** 4) * np.exp(-( (_x + x_c) / w_c) ** 2 ) \

def df5(_x, a = 3/4., b = 1., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_b = 0., x_c = -0.25):

    return a * (-2./w_a ** 2) ** 2 * ( 3*(1 + 2*(-2./w_a**2)*(_x + x_a) ** 2) + (-2. / w_a**2)**2 * (_x + x_a) ** 4) * (-2./w_a**2) * (_x + x_a) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
      + a * (-2./w_a ** 2) ** 2 * ( 3*(0 + 2*2*(-2./w_a**2)*(_x + x_a)) + 4* (-2. / w_a**2)**2 * (_x + x_a) ** 3) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
      + b * (-2./w_b ** 2) ** 2 * ( 3*(1 + 2*(-2./w_b**2)*(_x + x_b) ** 2) + (-2. / w_b**2)**2 * (_x + x_b) ** 4) * (-2./w_b**2) * (_x + x_b) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
      + b * (-2./w_b ** 2) ** 2 * ( 3*(0 + 2*2*(-2./w_b**2)*(_x + x_b)) + 4* (-2. / w_b**2)**2 * (_x + x_b) ** 3) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
      + c * (-2./w_c ** 2) ** 2 * ( 3*(1 + 2*(-2./w_c**2)*(_x + x_c) ** 2) + (-2. / w_c**2)**2 * (_x + x_c) ** 4) * (-2./w_c**2) * (_x + x_c) * np.exp(-( (_x + x_c) / w_c) ** 2 ) \
      + c * (-2./w_c ** 2) ** 2 * ( 3*(0 + 2*2*(-2./w_c**2)*(_x + x_c)) + 4* (-2. / w_c**2)**2 * (_x + x_c) ** 3) * np.exp(-( (_x + x_c) / w_c) ** 2 )

def domain(_a = -0.5, _b = 1.5, _Nx = 21):

    _L = float(_b - _a)
    _dx = _L / _Nx

    _x = np.zeros(_Nx)
    for i in range(_Nx):
        _x[i] = _a + i*_dx

    return _x, _dx, _L


def FD_derivative_matrix_formulation(_dn = 1, _p = 1, _Nx = 128):

    x, dx, L = domain(_a = -0.5, _b = 1.5, _Nx = _Nx)

    f = function(x)
    df_str = 'df' + str(_dn)
    df = eval(df_str)(x)

    # compute derivative approximations
    # extract subdictionary pertaining to scheme of derivative order dn
    dn_key = 'dn' + str(_dn)
    FD_schemes = pyfiles.lib.make_FD_schemes_dict.store(_dn, _p)
    FD_schemes = FD_schemes[dn_key] # dictionary containing the family of
                                    # FD schemes or order dn

    imax = len(x) - 1
    stencil_size = _p + _dn
    stencil_center = stencil_size // 2

    W = np.zeros([_Nx, _Nx])
    for i in range(_Nx):
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

        FD_scheme = FD_schemes[handedness][asymmetry]
        w = FD_scheme['w']
        stencil = FD_scheme['stencil']

        W[i, i + np.array(stencil)] = w

    df_approx = W.dot(f)
    df_approx /= dx ** _dn
    error = df_approx - df
    L2_error_norm = np.sqrt(dx/L)* LA.norm(error, 2) # L2 norm of error
    Linf_error = np.max(error) # L-infinity norm of error

    return L2_error_norm, df_approx

def convergence_routine_nonperiodic(NumGrids = 10, Nx = 18, LTE = 3, dn = 1, plots = 'yes'):
    """executes FD_derivative_matrix_formulation(*args) on progressively finer grids.
    The first grid contains Nx points (must be at least as large as the stencil size
    = O(LTE) + order of derivative), the subsequent (NumGrids-1)
    grids double the number of points from the previous grid

    inputs:
    NumGrids -- (int) number of grids the derivative is to be evaluated on
    Nx -- (int) number of gridpoints on coarsest grid
    LTE -- (int) local truncation error of the derivatives being used
    dn -- (int) order of derivative to be evaluated

    NOTE: derivative tables must be generated for the required dn at chosen LTE
    by executing

   $ python generate_tables_of_finite_difference_schemes_for_a_given_LTE.py LTE dn

    from ./bin/ which stores the table in ./etc/
    """
    error_norm = np.zeros(NumGrids)
    grids = np.zeros(NumGrids)
    orders = np.zeros(NumGrids)

    # calculate first order then loop over remaining
    q = 0
    error_norm[q], df_approx = FD_derivative_matrix_formulation(dn, LTE, Nx)
    grids[q] = Nx

    for q in range(1, NumGrids):
        Nx *= 2
        error_norm[q], df_approx = FD_derivative_matrix_formulation(dn, LTE, Nx)
        orders[q] = np.log2(error_norm[q-1] / error_norm[q])
        grids[q] = Nx

    order = -np.polyfit(np.log2(grids), np.log2(error_norm), 1)[0]
    print "slope of a linear fit = %g \n" % order

    print "order calculations at each refinement step:"

    for n in range(NumGrids):
        if n == 0:
            print "Nx%d        error = %g       ----" % (grids[n],error_norm[n])
        else:
            print "Nx%d        error = %g       order = %g" % (grids[n],error_norm[n],orders[n])

    print '\n'

    if plots.lower() == 'yes':
        fig, ax = plt.subplots()
        ax.loglog(grids, error_norm, '--o', markersize = 8, linewidth = 2)
        ax.set_xscale('log', basex=2)
        ax.set_yscale('log', basey=2)
        plt.xlabel(r'Number of grid points $N_x$ over $x\in [-1,1]$', fontsize = 16)
        plt.ylabel(r'$L^2$ error of $f_{\Delta x} - f_{exact}$', fontsize = 16)
        plt.text(np.max(grids), np.max(error_norm),'slope = -%g' % order, horizontalalignment = 'right', verticalalignment = 'top', fontsize = 16)
        plt.title(r'Error in $d^{%i}f/dx^{%i}$ using FD schemes at LTE = $O(\Delta x^{%i})$' % (dn,dn,LTE))
        plt.grid()

    return error_norm, grids


def convergence_for_several_derivatives_at_const_LTE(NumGrids = 10, _LTE = 3, _dn_min = 1, _dn_max = 4):
    """runs a convergence routine on a chosen number of grids (NumGrids) that whose uniform
    spacing between grid points is halved in each subsequent run, for the range of derivatives
    dn_min <= dn <= dn_max, dn_min must be greater than or equal to 1"""

    error_histories = np.zeros([NumGrids, _dn_max+1]) # no values stored in [:,0] entry
    grid_histories = np.zeros([NumGrids, _dn_max+1])  # no values stored in [:,0] entry

    fig, ax = plt.subplots()

    for d in range(_dn_min, _dn_max + 1):
        print "derivative dn = %d: \n" % d

        error_histories[:,d], grid_histories[:,d] = convergence_routine_nonperiodic(NumGrids = 10,
                                                    Nx = 21, LTE = _LTE, dn = d, plots = 'no')

        plot_label = '$dn = $' + str(d)
        ax.loglog(grid_histories[:,d], error_histories[:,d], '--o', markersize = 8, linewidth = 2, label = plot_label)
        ax.hold('on')
    ax.set_xscale('log', basex=2)
    ax.set_yscale('log', basey=2)
    plt.xlabel(r'Number of grid points $N_x$ over $x\in [-1,1]$', fontsize = 16)
    plt.ylabel('$L^2$ error', fontsize = 16)
    plt.legend(loc = 'best')
    plt.grid()

if __name__ == '__main__':
    """
    run from terminal as

        $ python main.py NumGrids Nx LTE dn

    sys.argv inputs:
    NumGrids -- (int) number of grids the derivative is to be evaluated on
    Nx -- (int) number of gridpoints on coarsest grid
    LTE -- (int) local truncation error of the derivatives being used
    dn -- (int) order of derivative to be evaluated

    NOTE: derivative tables must be generated for the required dn at chosen LTE
    by executing

   $ python generate_tables_of_finite_difference_schemes_for_a_given_LTE.py LTE dn

    from ./bin/ which stores the table in ./etc/
    """
    import sys

    NumGrids = int(sys.argv[1])
    Nx = int(sys.argv[2])
    LTE = int(sys.argv[3])
    dn = int(sys.argv[4])

    error_norm = np.zeros(NumGrids)
    grids = np.zeros(NumGrids)
    orders = np.zeros(NumGrids)
    p = LTE # relabeling is due to author's personal collocation with this symbol

    # calculate first order then loop over remaining
    q = 0
    error_norm[q] = FD_derivative_matrix_formulation(dn, p, Nx)
    grids[q] = Nx

    for q in range(1, NumGrids):
        Nx *= 2
        error_norm[q] = FD_derivative_matrix_formulation(dn, p, Nx)
        orders[q] = np.log2(error_norm[q-1] / error_norm[q])
        grids[q] = Nx

    order = -np.polyfit(np.log2(grids), np.log2(error_norm), 1)[0]
    print "slope of a linear fit = %g" % order

    print "order calculations at each refinement step:"

    for n in range(NumGrids):
        if n == 0:
            print "Nx%d        error = %g       ----" % (grids[n],error_norm[n])
        else:
            print "Nx%d        error = %g       order = %g" % (grids[n],error_norm[n],orders[n])

    fig, ax = plt.subplots()
    ax.loglog(grids, error_norm, '--o')
    #ax.hold('on')
    #ax.loglog(abscissa, order_line, '-b')
    ax.set_xscale('log', basex=2)
    ax.set_yscale('log', basey=2)
    #ax.set_xlim(1, 1e4)
    #ax.set_ylim(1e-15, 1)
    plt.ylabel('$L^2$ norm of $f_{num} - f_{exact}$', fontsize = 14)
    plt.xlabel('number of gridpoints $N_x$', fontsize = 14)
    plt.text(np.max(grids), np.max(error_norm),'slope = -%g' % order, horizontalalignment = 'right', verticalalignment = 'top', fontsize = 16)
    plt.grid()
    plt.show()
