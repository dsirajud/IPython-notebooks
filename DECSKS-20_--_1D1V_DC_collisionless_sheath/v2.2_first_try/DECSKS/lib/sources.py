import numpy as np

def maxwellian(v, tau = 1/30., mu = 1836.15267389, vD = 0., species = 'electron'):

    if species == 'electron':

        S = np.where(v.gridvalues < 0, 1 / np.sqrt(2*np.pi) * np.exp(-1/2. * (v.gridvalues - vD) ** 2), 0)

    elif species == 'ion':
        # vD for ions -- (float) is measured in multiples of the sound speed here, which in these normalized units
        #  are multiples of 1 / sqrt(mu), see DECSKS-20 for the details

        vD *= 1 / np.sqrt(mu)

        # multiples of the sound speed in these normalized units
        S = np.where(v.gridvalues < 0, 1 / np.sqrt(2*np.pi*tau/mu) * np.exp(-(v.gridvalues - vD) ** 2 / (2 * tau / mu)), 0)
    
    return S

