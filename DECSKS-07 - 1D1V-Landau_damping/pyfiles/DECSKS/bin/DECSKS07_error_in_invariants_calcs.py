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

def main(scheme = 'fourier'):
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

    if scheme.lower() == 'fourier':
        outdir_coarse = outdir_parent + 's7-01/'
        outdir_fine = outdir_parent + 's7-02/'
    else:
        outdir_coarse= outdir_parent + 's7-03/'
        outdir_fine = outdir_parent + 's7-04/'
        outdir_Nx768 = outdir_parent + 's7-05/'
        outdir_s7_07 = outdir_parent + 's7-07/'

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

    # Nx768 FD simulation, Nv = 256, LF2, Nt = 60

    filepath_I1_Nx768 = outdir_Nx768 + filename_I1
    filepath_I2_Nx768 = outdir_Nx768 + filename_I2
    filepath_IW_Nx768 = outdir_Nx768 + filename_IW
    filepath_IS_Nx768 = outdir_Nx768 + filename_IS

    # s7_07 FD simulation, Nv = 128, O6-4, Nt = 60
    filepath_I1_s7_07 = outdir_s7_07 + filename_I1
    filepath_I2_s7_07 = outdir_s7_07 + filename_I2
    filepath_IW_s7_07 = outdir_s7_07 + filename_IW
    filepath_IS_s7_07 = outdir_s7_07 + filename_IS

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

    # Nx768 FD simulation, Nv = 256, LF2, Nt = 60
    I1_Nx768 = dat2array(filepath_I1_Nx768)
    I2_Nx768 = dat2array(filepath_I2_Nx768)
    IW_Nx768 = dat2array(filepath_IW_Nx768)
    IS_Nx768 = dat2array(filepath_IS_Nx768)

    # s7_07 FD simulation, Nv = 128, O6-4, Nt = 60
    I1_s7_07 = dat2array(filepath_I1_s7_07)
    I2_s7_07 = dat2array(filepath_I2_s7_07)
    IW_s7_07 = dat2array(filepath_IW_s7_07)
    IS_s7_07 = dat2array(filepath_IS_s7_07)

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

    # Nx768 FD simulation, Nv = 256, LF2, Nt = 60
    I1_Nx768_error = fractional_error(I1_Nx768)
    I2_Nx768_error = fractional_error(I2_Nx768)
    IW_Nx768_error = fractional_error(IW_Nx768)
    IS_Nx768_error = fractional_error(IS_Nx768)

    # s7_07 FD simulation, Nv = 128, O6-4, Nt = 60
    I1_s7_07_error = fractional_error(I1_s7_07)
    I2_s7_07_error = fractional_error(I2_s7_07)
    IW_s7_07_error = fractional_error(IW_s7_07)
    IS_s7_07_error = fractional_error(IS_s7_07)
    
    # plot specs

    T = 60.

    if scheme.lower() == 'fourier':
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

    else:
        Nt_coarse = 60
        dt_coarse = T / Nt_coarse

        Nt_fine = 200
        dt_fine = T / Nt_fine

        # abscissas
        t_coarse = np.zeros(Nt_coarse + 1)
        t_fine = np.zeros(Nt_fine + 1)

        for it in range(Nt_coarse + 1):
            t_coarse[it] = 0 + it * dt_coarse

        for it in range(Nt_fine + 1):
            t_fine[it] = 0 + it * dt_fine

        t_Nx768 = t_coarse
        t_s7_07 = t_coarse


    # plot tableau

    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
    for i in range(len(tableau20)):
        r, g, b = tableau20[i]
        tableau20[i] = (r / 255., g / 255., b / 255.)

    # reset dictionary values for legend
    params = {'legend.fontsize': 10,
              'legend.linewidth': 2}
    plt.rcParams.update(params)
    plt.figure(0)

    if scheme.lower() == 'fourier':
        plt.semilogy(t_coarse, np.abs(I1_coarse_error), marker = 'o', color = tableau20[0],
                 label = r'Fourier/O6-4 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^5)$, $N_x = 8,\, N_v = 256$')
        plt.semilogy(t_fine, np.abs(I1_fine_error), color = tableau20[5],
                 label = r'Fourier/O11-6 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^7)$, $N_x = 16,\, N_v = 512$')
    else:
        plt.semilogy(t_coarse, np.abs(I1_coarse_error), 'o', color = tableau20[0],
                 label = r'FD7/LF2 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^3)$, $N_x = 384,\, N_v = 256$')
        plt.semilogy(t_fine, np.abs(I1_fine_error), color = tableau20[5],
                 label = r'FD7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 384,\, N_v = 256$')
        plt.semilogy(t_Nx768, np.abs(I1_Nx768_error), marker = 'H', color = tableau20[2],
                 label = r'FD7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 768,\, N_v = 256$')
        #        plt.semilogy(t_s7_07, np.abs(I1_s7_07_error), '>g',
        #                 label = r'FD7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 1536,\, N_v = 128$')

    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Abs. value of fractional error in the $L^1$ norm', fontsize = 12)
    plt.title('Fractional error in the $L^1$ norm invariant (mass conservation)', fontsize = 14)
    plt.legend(loc = 'best')
    plt.axis([0, 60, 1e-16, 1e-6])
    if scheme.lower() == 'fourier':
        plt.savefig('../../../fig/I1_-_F21.png')
    else:
        plt.savefig('../../../fig/I1_-_FD7.png')
    plt.figure(1)

    if scheme.lower() == 'fourier':
        plt.semilogy(t_coarse, np.abs(I2_coarse_error), marker = 'o', color = tableau20[0],
                 label = r'Fourier/O6-4 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^5)$, $N_x = 8,\, N_v = 256$')
        plt.semilogy(t_fine, np.abs(I2_fine_error), color = tableau20[5],
                 label = r'Fourier/O11-6 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^7)$, $N_x = 16,\, N_v = 512$')

    else:
        plt.semilogy(t_coarse, np.abs(I2_coarse_error), marker = 'o', color = tableau20[0],
                 label = r'FD7/LF2 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^3)$, $N_x = 384,\, N_v = 256$')
        plt.semilogy(t_fine, np.abs(I2_fine_error), marker = 'o', color = tableau20[5],
                 label = r'FD7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 384,\, N_v = 256$')
        plt.semilogy(t_Nx768, np.abs(I2_Nx768_error), marker = 'H', color = tableau20[2],
                 label = r'FD7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 768,\, N_v = 256$')
        #        plt.semilogy(t_s7_07, np.abs(I2_s7_07_error), '>g',
        #                 label = r'FD7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 1536,\, N_v = 128$')


    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Abs. value of fractional error in the $L^2$ norm', fontsize = 12)
    plt.title('Fractional error in the $L^2$ norm invariant', fontsize = 14)
    plt.legend(loc = 'best')
    if scheme.lower() == 'fourier':
        plt.axis([0, 60, 1e-16, 1e-6])
        plt.savefig('../../../fig/I2_-_F21.png')
    else:
        plt.axis([0,60,1e-16, 1e-1])
        plt.savefig('../../../fig/I2_-_FD7.png')

    plt.figure(2)
    if scheme.lower() == 'fourier':
        plt.semilogy(t_coarse, np.abs(IW_coarse_error), marker = 'o', color = tableau20[0],
                 label = r'Fourier/O6-4 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^5)$, $N_x = 8,\, N_v = 256$')
        plt.semilogy(t_fine, np.abs(IW_fine_error), color = tableau20[5],
             label = r'Fourier/O11-6 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^7)$, $N_x = 16,\, N_v = 512$')

    else:
        plt.semilogy(t_coarse, np.abs(IW_coarse_error), marker = 'o', color = tableau20[0],
                 label = r'FD7/LF2 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^3)$, $N_x = 384,\, N_v = 256$')
        plt.semilogy(t_fine, np.abs(IW_fine_error), marker = 'o', color = tableau20[5],
                 label = r'FD7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 384,\, N_v = 256$')
        plt.semilogy(t_Nx768, np.abs(IW_Nx768_error), marker = 'H', color = tableau20[2],
                 label = r'FD7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 768,\, N_v = 256$')
        #        plt.semilogy(t_s7_07, np.abs(IW_s7_07_error), '>g',
        #                 label = r'FD7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 1536,\, N_v = 128$')


    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Abs. value of fractional error in total energy, $W$', fontsize = 12)
    plt.title('Fractional error in the total energy invariant $W$', fontsize = 14)
    plt.legend(loc = 'best')
    if scheme.lower() == 'fourier':
        plt.axis([0, 60, 1e-16, 1e-6])
        plt.savefig('../../../fig/IW_-_F21.png')
    else:
        plt.axis([0, 60, 1e-16, 1e-1])
        plt.savefig('../../../fig/IW_-_FD7.png')


    plt.figure(3)
    if scheme.lower() == 'fourier':
        plt.semilogy(t_coarse, np.abs(IS_coarse_error), marker = 'o', color = tableau20[0],
                 label = r'Fourier/O6-4 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^5)$, $N_x = 8,\, N_v = 256$')
        plt.semilogy(t_fine, np.abs(IS_fine_error), color = tableau20[5],
                 label = r'Fourier/O11-6 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^7)$, $N_x = 16,\, N_v = 512$')

    else:
        plt.semilogy(t_coarse, np.abs(IS_coarse_error), marker = 'o', color = tableau20[0],
                 label = r'FD7/LF2 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^3)$, $N_x = 384,\, N_v = 256$')
        plt.semilogy(t_fine, np.abs(IS_fine_error), marker = 'o', color = tableau20[5],
                 label = r'FD7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 384,\, N_v = 256$')
        plt.semilogy(t_Nx768, np.abs(IS_Nx768_error), marker = 'H', color = tableau20[2],
                 label = r'FD7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 768,\, N_v = 256$')
        #        plt.semilogy(t_s7_07, np.abs(IS_s7_07_error), '>g',
        #                 label = r'FD7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 1536,\, N_v = 128$')


    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Abs. value of fractional error in the entropy $S$ invariant', fontsize = 12)
    plt.title('Fractional error in the entropy invariant $S$', fontsize = 14)
    plt.legend(loc = 'best')
    if scheme.lower() == 'fourier':
        plt.axis([0, 60, 1e-16, 1e-6])
        plt.savefig('../../../fig/IS_-_F21.png')
    else:
        plt.axis([0,60, 1e-16, 1e-1])
        plt.savefig('../../../fig/IS_-_FD7.png')
    return None
