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
    """main routine for electrostatic energy  history calculations in Landau
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
        outdir_Nx768_F7 = outdir_parent + 's7-06/'

    # common filenames
    filename_WE = 'out_WE'

    filepath_WE_coarse = outdir_coarse + filename_WE
    filepath_WE_fine = outdir_fine + filename_WE
    filepath_WE_Nx768 = outdir_Nx768 + filename_WE
    filepath_WE_Nx768_F7 = outdir_Nx768_F7 + filename_WE

    # store arrays from dat files
    WE_coarse = dat2array(filepath_WE_coarse)
    WE_fine = dat2array(filepath_WE_fine)
    WE_Nx768 = dat2array(filepath_WE_Nx768)
    WE_Nx768_F7 = dat2array(filepath_WE_Nx768_F7)

    # normalize each electrostatic energy by the time zero value

    WE_coarse /= WE_coarse[0]
    WE_fine /= WE_fine[0]
    WE_Nx768 /= WE_Nx768[0]
    WE_Nx768_F7 /= WE_Nx768_F7[0]

    # plot specs
    T = 60.

    Nt_coarse = 60
    dt_coarse = T / Nt_coarse

    if scheme.lower() == 'fourier':
        Nt_fine = 600
        dt_fine = T / Nt_fine

    else:
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
    t_Nx768_F7 = t_coarse

    # plots
    # reset dictionary values for legend
    params = {'legend.fontsize': 11,
              'legend.linewidth': 2}
    plt.rcParams.update(params)

    # Linear Landau damping prediction:

    gamma = -0.153359 # from linear theory for this distribution function
    WE_linear = np.exp(2*gamma*t_coarse) # Linear Landau theory prediction for normalized W_E / W_E0


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

    if scheme.lower() == 'fourier':
        plt.semilogy(t_coarse, WE_coarse, 'o',
                     label = r'Fourier/O6-4 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^5)$, $N_x = 8,\, N_v = 256$')
        plt.semilogy(t_fine, WE_fine, '-g',
                     label = r'Fourier/O11-6 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^7)$, $N_x = 16,\, N_v = 512$')

    else:
        #        plt.semilogy(t_coarse, WE_coarse, marker = 'o', color = tableau20[1],
        #                     label = r'FD7/LF2 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^2)$, $N_x = 384,\, N_v = 256$')
        #        plt.semilogy(t_fine, WE_fine, marker = 'D', color = tableau20[0], linestyle = 'None',
        #                     label = r'FD7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 384,\, N_v = 512$')
        plt.semilogy(t_Nx768, WE_Nx768, marker = 'H', color = tableau20[2],
                     label = r'FD7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 768,\, N_v = 256$')
        plt.semilogy(t_Nx768_F7, WE_Nx768_F7, marker = '<', color = tableau20[3],
                     label = r'F7/O6-4 $LTE = \mathcal{O}(\Delta x^{8},\,\Delta v^{8}, \, \Delta t^5)$, $N_x = 16,\, N_v = 256$')

    plt.semilogy(t_coarse, WE_linear, '--k', linewidth = 2, label = 'Linear theory')

    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Normalized electrostatic energy, $W_E/W_{E0}$', fontsize = 12)
    plt.title('Normalized electrostatic energy, $W_E / W_{E0}$')
    plt.legend(loc = 'best')


    if scheme.lower() == 'fourier':
        plt.savefig('../../../fig/WE_-_F21.png')
    else:
        plt.savefig('../../../fig/WE_-_FD7.png')
    return None
