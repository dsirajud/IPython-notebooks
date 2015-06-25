import numpy as np
import matplotlib.pyplot as plt

def dat2array(filepath):
    """reads in a dat file formatted where the delimiter for each
    datum is a newline \n character and returns a numpy array
    of data

    input:
    filepath -- (str) relative location of the dat file

    output:
    data -- (ndarray, ndim = 1) numpy array of data read in from file
    """
    infile = open(filepath)
    lines = infile.readlines()

    data = []
    for line in lines:
        data.append(eval(line))

    return np.array(data)

def fractional_error(data):
    """computes and returns a list of fractional error of data (list)
    with respect to its value at time zero (data[0])

    input:
    data -- (ndarray, ndim = 1) numpy array of data

    output:
    fractional_error_data -- (ndarray, ndim = 1) numpy array containing the fractional
                             deviation of each entry of data with
                             respect to its value at data[0]
    """

    fractional_error_data = np.zeros(len(data))
    fractional_error_data[0] = 0.    # no deviation at data[0] from itself by definition
    fractional_error_data[1:] = (data[1:] - data[0]) / data[0]

    return fractional_error_data

def main(scheme = 'fourier', N = 21):
    """main routine for error history calculations in Landau
    simulations, here we have conducted a 'coarse' and
    'fine' simulation and the errors for each invariant
    are plotted on the same axis for comparison in separate
    figures.

    input:
    scheme -- (str) 'fourier' or 'FD'
    N -- (int) global error

    output:
    display plots in jupyter notebook inline, but return None

    note: schemes employed so far:

    F21
    FD7
    """
    # directory strings
    outdir_parent = './../etc/outputs/'
    outdir_coarse_str = 'Landau_--_Nx8Nv256'
    outdir_fine_str = 'Landau_--_Nx16Nv512'

    if scheme.lower() == 'fourier':
        scheme_str = 'F' + str(N)
    else:
        scheme_str = 'FD' + str(N)

    outdir_coarse_time_splitter_stem = '_O6-4/'
    outdir_fine_time_splitter_stem = '_O11-6/'

    outdir_coarse = outdir_parent + outdir_coarse_str + scheme_str + outdir_coarse_time_splitter_stem
    outdir_fine = outdir_parent + outdir_fine_str + scheme_str + outdir_fine_time_splitter_stem

    # common filenames
    filename_I1 = 'out_I1'
    filename_I2 = 'out_I2'
    filename_IW = 'out_IW'
    filename_IS = 'out_S'

    # coarse simulation
    filepath_I1_coarse = outdir_coarse + filename_I1
    filepath_I2_coarse = outdir_coarse + filename_I2
    filepath_IW_coarse = outdir_coarse + filename_IW
    filepath_IS_coarse = outdir_coarse + filename_IS

    # fine simulation
    filepath_I1_fine = outdir_fine + filename_I1
    filepath_I2_fine = outdir_fine + filename_I2
    filepath_IW_fine = outdir_fine + filename_IW
    filepath_IS_fine = outdir_fine + filename_IS

    # store arrays from dat files

    # coarse simulation
    I1_coarse = dat2array(filepath_I1_coarse)
    I2_coarse = dat2array(filepath_I2_coarse)
    IW_coarse = dat2array(filepath_IW_coarse)
    IS_coarse = dat2array(filepath_IS_coarse)

    # fine simulation
    I1_fine = dat2array(filepath_I1_fine)
    I2_fine = dat2array(filepath_I2_fine)
    IW_fine = dat2array(filepath_IW_fine)
    IS_fine = dat2array(filepath_IS_fine)

    # compute fractional error arrays

    # coarse simulation
    I1_coarse_error = fractional_error(I1_coarse)
    I2_coarse_error = fractional_error(I2_coarse)
    IW_coarse_error = fractional_error(IW_coarse)
    IS_coarse_error = fractional_error(IS_coarse)

    # fine simulation
    I1_fine_error = fractional_error(I1_fine)
    I2_fine_error = fractional_error(I2_fine)
    IW_fine_error = fractional_error(IW_fine)
    IS_fine_error = fractional_error(IS_fine)

    # plot specs

    T = 60.

    Nt_coarse = 60
    dt_coarse = T / Nt_coarse

    Nt_fine = 600
    dt_fine = T / Nt_fine

    # abscissas
    t_coarse = np.zeros(Nt_coarse + 1)
    t_fine = np.zeros(Nt_fine + 1)

    for it in range(Nt_coarse + 1):
        t_coarse[it] = 0 + it * dt_coarse

    for it in range(Nt_fine + 1):
        t_fine[it] = 0 + it * dt_fine

    # reset dictionary values for legend
    params = {'legend.fontsize': 12,
              'legend.linewidth': 2}
    plt.rcParams.update(params)

    plt.figure(0)
    plt.semilogy(t_coarse, np.abs(I1_coarse_error), 'or',
                 label = r'$GE = \mathcal{O}(\Delta x^{21},\,\Delta v^{21}, \, \Delta t^4)$, $N_x = 8,\, N_v = 256$')
    plt.semilogy(t_fine, np.abs(I1_fine_error), 'b',
                 label = r'$GE = \mathcal{O}(\Delta x^{21},\,\Delta v^{21}, \, \Delta t^6)$, $N_x = 16,\, N_v = 512$')

    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Abs. value of fractional error in the $L^1$ norm', fontsize = 12)
    plt.title('Fractional error in the $L^1$ norm invariant (mass conservation)', fontsize = 14)
    plt.legend(loc = 'best')
    plt.axis([0, 60, 1e-16, 1e-6])

    plt.figure(1)
    plt.semilogy(t_coarse, np.abs(I2_coarse_error), 'or',
                 label = r'$GE = \mathcal{O}(\Delta x^{21},\,\Delta v^{21}, \, \Delta t^4)$, $N_x = 8,\, N_v = 256$')
    plt.semilogy(t_fine, np.abs(I2_fine_error), 'b',
                 label = r'$GE = \mathcal{O}(\Delta x^{21},\,\Delta v^{21}, \, \Delta t^6)$, $N_x = 16,\, N_v = 512$')

    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Abs. value of fractional error in the $L^2$ norm', fontsize = 12)
    plt.title('Fractional error in the $L^2$ norm invariant', fontsize = 14)
    plt.legend(loc = 'best')
    plt.axis([0, 60, 1e-16, 1e-6])


    plt.figure(2)
    plt.semilogy(t_coarse, np.abs(IW_coarse_error), 'or',
                 label = r'$GE = \mathcal{O}(\Delta x^{21},\,\Delta v^{21}, \, \Delta t^4)$, $N_x = 8,\, N_v = 256$')
    plt.semilogy(t_fine, np.abs(IW_fine_error), 'b',
                 label = r'$GE = \mathcal{O}(\Delta x^{21},\,\Delta v^{21}, \, \Delta t^6)$, $N_x = 16,\, N_v = 512$')

    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Abs. value of fractional error in total energy, $W$', fontsize = 12)
    plt.title('Fractional error in the total energy invariant $W$', fontsize = 14)
    plt.legend(loc = 'best')
    plt.axis([0, 60, 1e-16, 1e-6])

    plt.figure(3)
    plt.semilogy(t_coarse, np.abs(IS_coarse_error), 'or',
                 label = r'$GE = \mathcal{O}(\Delta x^{21},\,\Delta v^{21}, \, \Delta t^4)$, $N_x = 8,\, N_v = 256$')
    plt.semilogy(t_fine, np.abs(IS_fine_error), 'b',
                 label = r'$GE = \mathcal{O}(\Delta x^{21},\,\Delta v^{21}, \, \Delta t^6)$, $N_x = 16,\, N_v = 512$')

    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Abs. value of fractional error in the entropy $S$ invariant', fontsize = 12)
    plt.title('Fractional error in the entropy invariant $S$', fontsize = 14)
    plt.legend(loc = 'best')
    plt.axis([0, 60, 1e-16, 1e-6])

    return None
