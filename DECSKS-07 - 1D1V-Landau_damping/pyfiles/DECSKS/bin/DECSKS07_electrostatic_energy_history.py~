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
    filename_WE = 'out_WE'

    filepath_WE_coarse = outdir_coarse + filename_WE
    filepath_WE_fine = outdir_fine + filename_WE

    # store arrays from dat files
    WE_coarse = dat2array(filepath_WE_coarse)
    WE_fine = dat2array(filepath_WE_fine)

    # normalize each electrostatic energy by the time zero value

    WE_coarse /= WE_coarse[0]
    WE_fine /= WE_fine[0]

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


    # plots
    # reset dictionary values for legend
    params = {'legend.fontsize': 14,
              'legend.linewidth': 2}
    plt.rcParams.update(params)

    # Linear Landau damping prediction:

    gamma = -0.153359 # from linear theory for this distribution function
    WE_linear = np.exp(2*gamma*t_coarse) # Linear Landau theory prediction for normalized W_E / W_E0


    fig, ax = plt.subplots()
    ax.semilogy(t_coarse, WE_coarse, 'o',
                 label = r'$GE = \mathcal{O}(\Delta x^{21},\,\Delta v^{21}, \, \Delta t^4)$, $N_x = 8,\, N_v = 256$')
    ax.semilogy(t_fine, WE_fine, '-g',
                 label = r'$GE = \mathcal{O}(\Delta x^{21},\,\Delta v^{21}, \, \Delta t^6)$, $N_x = 16,\, N_v = 512$')
    ax.semilogy(t_coarse, WE_linear, '--m', linewidth = 2, label = 'Linear theory')

    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Normalized electrostatic energy, $W_E/W_{E0}$', fontsize = 12)
    plt.title('Normalized electrostatic energy, $W_E / W_{E0}$')

    # shrink frame of the plot to make room for a legend on the right side
    frame = ax.get_position() # position of plot center
    ax.set_position([frame.x0, frame.y0,
                     frame.width * 0.8, frame.height])

    # Place legend to the right of the shrunken frame
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
              fancybox=True, ncol=1)



    return None
