import numpy as np
import pylab
from mpl_toolkits.axes_grid1 import make_axes_locatable

class Plots:
    """Plots are generated in ./DECSKS/plots/
    simulation will ask if you want to remove ALL
    plot files in this folder. Move any files from
    this folder that wish to be kept"""
    def __init__(self, t, x, v, it, sim_params):
        self.divider = '_-_'
        self.underscore = '_'

        self.Nx_str = 'Nx%d' % x.Ngridpoints
        self.Nv_str = 'Nv%d' % v.Ngridpoints
        self.Nt_str = 'Nt%d' % t.N

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
    def __init__(self, f, it, t, x, v, sim_params, species = None):
        Plots.__init__(self, t, x, v, it, sim_params)

        plot_params = sim_params['plot_params']
        self.xmin = plot_params['xmin']
        self.xmax = plot_params['xmax']
        self.ymin = plot_params['ymin']
        self.ymax = plot_params['ymax']

        if len(f.shape) == 2: # f = f[x,vx], 2 dim in phase space
            self.dimensionality = '1D1V'
            self.splitscheme = sim_params['split_scheme']
            self.filename = self.filetype + self.divider \
              + self.dimensionality + self.underscore \
              + species + self.underscore \
              + self.splitscheme + self.underscore \
              + self.Nx_str + self.Nv_str + self.Nt_str \
              + self.underscore + self.it_str + self.fileformat

            self.X, self.V = np.meshgrid(self.x.gridvalues,self.v.gridvalues)
            self.f = f

        elif len(f.shape) == 1:
            self.dimensionality = '1D'
            self.splitscheme = ''
            self.filename = self.filetype + self.divider \
              + self.dimensionality + self.underscore \
              + self.splitscheme + self.underscore \
              + self.Nx_str + self.Nv_str + self.Nt_str \
              + self.underscore + self.it_str + self.fileformat

            self.f = f

    def __call__(self, n):
        if len(self.f.shape) == 2:
            # f = f[x,vx], 2 dim in phase space
            pylab.pcolormesh(self.X, self.V, self.f.T, cmap = 'jet')
            pylab.colorbar()

            #            pylab.clim(0.0, 0.45)
            # clim for various test cases
            pylab.clim(0,0.38) # for bump on tail test cases (s17-01, 02, 03, 04)
                # pylab.clim(-0.004, 0.004) # for linear landau test cases (s07-01, 02, 03, 04, 05)
                # pylab.clim(0.0, 0.45) # for two stream instability (s18-18)
                # pylab.clim(0.0, 0.6) # for strong landau test cases (s18-17)


            pylab.grid()
            pylab.axis([self.xmin, self.xmax, self.ymin, self.ymax])
            pylab.xlabel('$x$', fontsize = 18)
            pylab.ylabel('$v_x$', fontsize = 18)
            pylab.title('sXX-XX DECSKS-2.3: $N_x$ = %d, $N_v$ = %d, $t^n$ = %2.3f, n = %03d' % (self.x.Ngridpoints, self.v.Ngridpoints, self.it*self.t.width, n))
            pylab.savefig(self.path + 'sXX-XX_' + self.filename)
            pylab.clf()
            return None

        if len(self.f.shape) == 1:
            # f = f[x], 1 dim in phase space
            pylab.plot(self.x.gridvalues,self.f,'ob')
            pylab.grid()
            pylab.axis([self.xmin, self.xmax, self.ymin, self.ymax])
            pylab.xlabel('$x$', fontsize = 18)
            pylab.ylabel('$f(x)$', fontsize = 18)
            pylab.savefig(self.path + self.filename)
            return None
