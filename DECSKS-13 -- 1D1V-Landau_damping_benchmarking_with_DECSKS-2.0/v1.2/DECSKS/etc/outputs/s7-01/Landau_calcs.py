import numpy as np
import pylab as plt

    # Plot the electrostatic energy of the system

    # WE
    # COARSE SOLUTION
    infile = open('out_WE','r')
    lines = infile.readlines()

    WE_coarse = []
    for line in lines:
        WE_coarse.append(eval(line))

    WE_coarse = np.array(WE_coarse)
    WE_coarse = WE_coarse / WE_coarse[0] # normalize

    Nt_coarse = 60
    dt_coarse = 1.0
    T_coarse = 60.

    t_coarse = np.zeros(Nt_coarse+1)

    for it in range(Nt_coarse+1):
        t_coarse[it] = 0 + it * dt_coarse

    # CONVERGED SOLUTION

    infile = open('./../Landau_--_Nx16Nv512F21_O11-6/out_WE','r')
    lines = infile.readlines()

    WE_converged = []
    for line in lines:
        WE_converged.append(eval(line))

    WE_converged = np.array(WE_converged)
    WE_converged = WE_converged / WE_converged[0] # normalize

    Nt_converged = 600
    dt_converged = .1
    T_converged = 60.

    t_converged = np.zeros(Nt_converged+1)

    for it in range(Nt_converged+1):
        t_converged[it] = 0 + it * dt_converged

    # LINEAR THEORY, WE = WE[0] exp(-2*gamma*t), if
    #                E = E[0] exp(jw_r t) exp(-gamma t)
    #                normalize: WE / WE[0] = exp(-gamma*t)

    gamma = 0.153359 # damping constant

    WE_linear = np.exp(-2*gamma*t_coarse)

    # plotting

    params = {'legend.fontsize': 10,
              'legend.linewidth': 2}
    plt.rcParams.update(params)

    plt.semilogy(t_coarse, WE_coarse, 'o',
                 label = r'Fourier/O6-4 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^5)$, $N_x = 8,\, N_v = 256$')
    plt.semilogy(t_converged, WE_converged, '-g',
                 label = r'Fourier/O11-6 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^7)$, $N_x = 16,\, N_v = 512$')
    plt.semilogy(t_coarse, WE_linear, '--m', linewidth = 2, label = 'Linear theory')
    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Normalized electrostatic energy, $W_E/W_{E0}$', fontsize = 12)
    plt.legend(loc = 'lower left')
    #plt.savefig('Landau_WE.png')
    plt.clf()


    # IW
    # COARSE SOLUTION
    infile = open('out_IW','r')
    lines = infile.readlines()

    IW_coarse = []
    for line in lines:
        IW_coarse.append(eval(line))

    IW_coarse = np.array(IW_coarse)
    IW_coarse_error = np.zeros(len(IW_coarse))

    for i in range(len(IW_coarse)):
        IW_coarse_error[i] = (IW_coarse[i] - IW_coarse[0]) / IW_coarse[0]
        IW_coarse_error[i] = np.where(np.sign(IW_coarse_error[i]) == -1, -IW_coarse_error[i], IW_coarse_error[i])

    # CONVERGED SOLUTION
    infile = open('./../Landau_--_Nx16Nv512F21_O11-6/out_IW','r')
    lines = infile.readlines()

    IW_converged = []
    for line in lines:
        IW_converged.append(eval(line))

    IW_converged = np.array(IW_converged)
    IW_converged_error = np.zeros(len(IW_converged))

    for i in range(len(IW_converged)):
        IW_converged_error[i] = (IW_converged[i] - IW_converged[0]) / IW_converged[0]
        IW_converged_error[i] = np.where(np.sign(IW_converged_error[i]) == -1, -IW_converged_error[i], IW_converged_error[i])


    params = {'legend.fontsize': 12,
              'legend.linewidth': 1}
    plt.rcParams.update(params)

    plt.semilogy(t_coarse, IW_coarse_error, 'or',
                 label = r'Fourier/O6-4 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^5)$, $N_x = 8,\, N_v = 256$')
    plt.semilogy(t_converged, IW_converged_error, 'b',
                 label = r'Fourier/O11-6 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^7)$, $N_x = 16,\, N_v = 512$')

    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Relative error in total energy, $W$', fontsize = 12)
    plt.legend(loc = 'best')
    plt.axis([0, 60, 1e-16, 1e-6])
    #plt.savefig('Landau_IW.png')
    plt.clf()

    # S
    # COARSE SOLUTION
    infile = open('out_S','r')
    lines = infile.readlines()

    S_coarse = []
    for line in lines:
        S_coarse.append(eval(line))

    S_coarse = np.array(S_coarse)
    S_coarse_error = np.zeros(len(S_coarse))

    for i in range(len(S_coarse)):
        S_coarse_error[i] = (S_coarse[i] - S_coarse[0]) / S_coarse[0]

    # CONVERGED SOLUTION
    infile = open('./../Landau_--_Nx16Nv512F21_O11-6/out_S','r')
    lines = infile.readlines()

    S_converged = []
    for line in lines:
        S_converged.append(eval(line))

    S_converged = np.array(S_converged)
    S_converged_error = np.zeros(len(S_converged))

    for i in range(len(S_converged)):
        S_converged_error[i] = (S_converged[i] - S_converged[0]) / S_converged[0]
        S_converged_error[i] = np.where(np.sign(S_converged_error[i]) == -1, -S_converged_error[i], S_converged_error[i])

    params = {'legend.fontsize': 12,
              'legend.linewidth': 1}
    plt.rcParams.update(params)

    plt.semilogy(t_coarse, np.abs(S_coarse_error), 'or',
                 label = r'Fourier/O6-4 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^5)$, $N_x = 8,\, N_v = 256$')
    plt.semilogy(t_converged, np.abs(S_converged_error), 'b',
                 label = r'Fourier/O11-6 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^7)$, $N_x = 16,\, N_v = 512$')

    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Relative error in total entropy $S$', fontsize = 12)
    plt.axis([0, 60, 1e-16, 1e-6])
    plt.legend(loc = 'upper left')
    #plt.savefig('Landau_S.png')
    plt.clf()

    # I2
    # COARSE SOLUTION
    infile = open('out_I2','r')
    lines = infile.readlines()

    I2_coarse = []
    for line in lines:
        I2_coarse.append(eval(line))

    I2_coarse = np.array(I2_coarse)
    I2_coarse = np.sqrt(I2_coarse)
    I2_coarse_error = np.zeros(len(I2_coarse))

    for i in range(len(I2_coarse)):
        I2_coarse_error[i] = (I2_coarse[i] - I2_coarse[0]) / I2_coarse[0]

    # CONVERGED I2OLUTION
    infile = open('./../Landau_--_Nx16Nv512F21_O11-6/out_I2','r')
    lines = infile.readlines()

    I2_converged = []
    for line in lines:
        I2_converged.append(eval(line))

    I2_converged = np.array(I2_converged)
    I2_converged_error = np.zeros(len(I2_converged))

    for i in range(len(I2_converged)):
        I2_converged_error[i] = (I2_converged[i] - I2_converged[0]) / I2_converged[0]
        I2_converged_error[i] = np.where(np.sign(I2_converged_error[i]) == -1, -I2_converged_error[i], I2_converged_error[i])

    params = {'legend.fontsize': 12,
              'legend.linewidth': 1}
    plt.rcParams.update(params)

    plt.semilogy(t_coarse, np.abs(I2_coarse_error), 'or',
                 label = r'Fourier/O6-4 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^5)$, $N_x = 8,\, N_v = 256$')
    plt.semilogy(t_converged, np.abs(I2_converged_error), 'b',
                 label = r'Fourier/O11-6 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^7)$, $N_x = 16,\, N_v = 512$')

    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Relative error in $L^2$ norm', fontsize = 12)
    plt.axis([0, 60, 1e-16, 1e-6])
    plt.legend(loc = 'best')
    #plt.savefig('Landau_I2.png')
    plt.clf()

    # I1
    # COARSE SOLUTION
    infile = open('out_I1','r')
    lines = infile.readlines()

    I1_coarse = []
    for line in lines:
        I1_coarse.append(eval(line))

    I1_coarse = np.array(I1_coarse)
    I1_coarse_error = np.zeros(len(I1_coarse))

    for i in range(len(I1_coarse)):
        I1_coarse_error[i] = (I1_coarse[i] - I1_coarse[0]) / I1_coarse[0]

    # CONVERGED I1OLUTION
    infile = open('./../Landau_--_Nx16Nv512F21_O11-6/out_I1','r')
    lines = infile.readlines()

    I1_converged = []
    for line in lines:
        I1_converged.append(eval(line))

    I1_converged = np.array(I1_converged)
    I1_converged_error = np.zeros(len(I1_converged))

    for i in range(len(I1_converged)):
        I1_converged_error[i] = (I1_converged[i] - I1_converged[0]) / I1_converged[0]
        I1_converged_error[i] = np.where(np.sign(I1_converged_error[i]) == -1, -I1_converged_error[i], I1_converged_error[i])

    params = {'legend.fontsize': 12,
              'legend.linewidth': 1}
    plt.rcParams.update(params)

    plt.semilogy(t_coarse, np.abs(I1_coarse_error), 'or',
                 label = r'Fourier/O6-4 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^5)$, $N_x = 8,\, N_v = 256$')
    plt.semilogy(t_converged, np.abs(I1_converged_error), 'b',
                 label = r'Fourier/O11-6 $LTE = \mathcal{O}(\Delta x^{22},\,\Delta v^{22}, \, \Delta t^7)$, $N_x = 16,\, N_v = 512$')

    plt.grid()
    plt.xlabel('time', fontsize = 14)
    plt.ylabel('Relative error in $L^1$ norm', fontsize = 12)
    plt.legend(loc = 'best')
    plt.axis([0, 60, 1e-16, 1e-6])
    #plt.savefig('Landau_I1.png')
    plt.clf()

