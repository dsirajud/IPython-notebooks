import numpy as np
import pylab

class Plots:
    """Plots are generated in ./DECSKS/plots/
    simulation will ask if you want to remove ALL
    plot files in this folder. Move any files from
    this folder that wish to be kept"""
    def __init__(self, t, x, v, it, sim_params):
        self.divider = '_-_'
        self.underscore = '_'

        self.Nx_str = 'Nx%d' % x.N
        self.Nt_str = 'Nt%d' % sim_params['Nt']

        if sim_params['HOC'] == 'FOURIER':
            self.N_str  = 'F%d'  % sim_params['N']
        elif sim_params['HOC'] == 'FD':
            self.N_str  = 'FD%d' % sim_params['N']
        else:
            self.N_str = 'Classic_CS'

        self.it_str = 'it%05d' % it

        self.t = t
        self.x = x
        self.v = v
        self.it = it
        self.Nt = sim_params['Nt']

        self.filetype = 'plot' # as opposed to movie
        self.fileformat = '.png'
        self.path = './plots/'

class PlotSetup(Plots):
    def __init__(self, f, it, t, x, v, sim_params):
        Plots.__init__(self, t, x, v, it, sim_params)

        plot_params = sim_params['plot_params']
        self.xmin = plot_params['xmin']
        self.xmax = plot_params['xmax']
        self.ymin = plot_params['ymin']
        self.ymax = plot_params['ymax']

        self.dimensionality = '1D'
        self.filename = self.filetype + self.divider \
          + self.dimensionality + self.underscore \
          + self.N_str + self.Nx_str + self.Nt_str \
          + self.underscore + self.it_str + self.fileformat

        self.f = f

    def __call__(self, n):
        # f = f[x], 1 dim in phase space
        ft = self.f[n,:]
        pylab.plot(self.x.cells,ft,'ob')
        pylab.grid()
        pylab.axis([self.xmin, self.xmax, self.ymin, self.ymax])
        pylab.xlabel('$x$', fontsize = 18)
        pylab.ylabel('$f(x)$', fontsize = 18)
        pylab.savefig(self.path + self.filename)
        pylab.clf()
        return None
