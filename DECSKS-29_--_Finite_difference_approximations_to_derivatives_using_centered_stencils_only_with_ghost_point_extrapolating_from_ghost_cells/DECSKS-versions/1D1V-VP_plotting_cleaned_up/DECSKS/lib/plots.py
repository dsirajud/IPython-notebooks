import numpy as np
import pylab
from mpl_toolkits.axes_grid1 import make_axes_locatable

class Plots:
    """Plots are generated in ./DECSKS/plots/
    simulation will ask if you want to remove ALL
    plot files in this folder. Move any files from
    this folder that wish to be kept"""
    def __init__(self, t, x, v, sim_name, sim_params):

        # create common attributes for all plots

        # for filenames
        self.divider = '_--_'
        self.underscore = '_'
        self.sim_name = sim_name

        # to access gridvalues for plotting and x.N, v.N, t.width for labels
        self.t = t
        self.x = x
        self.v = v
        self.Nt = sim_params['Nt']

        self.fileformat = '.png'
        self.path = './plots/'

class PlotSetup(Plots):
    def __init__(self, f, t, x, v = None, sim_name = 'unspecified_simulation_number', sim_params = None, species = 'unspecified_species', quantity = 'unspecified_quantity'):
        Plots.__init__(self, t, x, v, sim_name, sim_params) # inherit superclass constructor attributes

        plot_params = sim_params['plot_params'] # dictionary of plot window specs
        self.xmin = plot_params['xmin']
        self.xmax = plot_params['xmax']
        self.ymin = plot_params['ymin']
        self.ymax = plot_params['ymax']

        self.species = species   # to label 2D plots title and filename
        self.quantity = quantity # to label 1D plots title and filename

        if f.ndim == 2: # f = f[x,vx], 2 dim in phase space
            self.filename_stem = self.sim_name + \
              self.divider + self.species + self.underscore \
              + quantity + self.underscore

            self.X, self.V = np.meshgrid(self.x.gridvalues,self.v.gridvalues)

        elif f.ndim == 1:
            self.filename_stem = self.sim_name + \
              self.divider + quantity + self.underscore

    def __call__(self, f, n, sim_name = ''):
        if f.ndim == 2:

            timestep_str = '%06d' % n

            filename = self.filename_stem + timestep_str

            plot_title = '[' + sim_name + '] ' + self.species + 's' \
              + ': $N_x$ = %d, $N_v$ = %d, $t^{%d}$ = %2.3f' \
              % (self.x.Ngridpoints, self.v.Ngridpoints, n, n*self.t.width)

            pylab.pcolormesh(self.X, self.V, f.T, cmap = 'jet')
            pylab.colorbar()

            #------------------------------------------------------------------------------------------
            # clim for various test cases
            #------------------------------------------------------------------------------------------
            # pylab.clim(0,0.38)              # for bump on tail test cases (s17-01, 02, 03, 04)
            # pylab.clim(-0.004, 0.004)       # for linear landau test cases (s07-01, 02, 03, 04, 05)
            # pylab.clim(0.0, 0.45)           # for two stream instability (s18-18)
            # pylab.clim(0.0, 0.6)            # for strong landau test cases (s18-17)
            #------------------------------------------------------------------------------------------

            pylab.clim(0,0.38)
            pylab.grid()
            pylab.axis([self.xmin, self.xmax, self.ymin, self.ymax])
            pylab.xlabel('$x$', fontsize = 18)
            pylab.ylabel('$v_x$', fontsize = 18)
            pylab.title(plot_title)
            pylab.savefig(self.path + filename +  self.fileformat)
            pylab.clf()
            return None

        if f.ndim == 1:

            timestep_str = '%06d' % n
            if self.quantity.lower() == 'potential':
                self.quantity = r'Potential $\phi$'

            filename = self.filename_stem + timestep_str

            plot_title = '[' + sim_name + '] ' + self.quantity \
              + ': $t^{%d}$ = %2.3f' % (n, n*self.t.width)

            # f = f[x], 1 dim in phase space
            pylab.plot(self.x.gridvalues, f, lw = 2, color = 'blue')
            pylab.grid()
            pylab.xlim(self.xmin, self.xmax)
            #            pylab.axis([self.xmin, self.xmax, self.ymin, self.ymax])
            pylab.xlabel('$x$', fontsize = 18)
            pylab.xlabel(r'position $x$', fontsize = 18)
            pylab.ylabel(self.quantity, fontsize = 18)
            pylab.title(plot_title)
            pylab.savefig(self.path + filename)
            pylab.clf()
            return None
