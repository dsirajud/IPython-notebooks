import matplotlib
import numpy as np


def densities(fe, fi, x, vx, t, n, sim_name, sim_params):
    # plotting densities
    ni = np.sum(fi, axis = 1)*vx.width
    ne = np.sum(fe, axis = 1)*vx.width

    # plot ni
    matplotlib.pyplot.plot(x.gridvalues, ni, linewidth = 2, color = 'red')
    matplotlib.pyplot.grid()
    #    matplotlib.pyplot.axis([x.gridvalues[0], x.gridvalues[-1], -5.1, 5.1])
    matplotlib.pyplot.xlabel(r'position $x$', fontsize = 18)
    matplotlib.pyplot.ylabel(r'$n_i (t^n,x)$', fontsize = 18)
    matplotlib.pyplot.title(sim_name + r' ion density $n_i (t^n,x)$: $N_x$ = %d, $N_v$ = %d, $t^n$ = %2.3f, n = %04d' % (sim_params['total_dims'][0], sim_params['total_dims'][1], n*t.width, n))
    it_str = 'it%05d' % n
    matplotlib.pyplot.savefig('./plots/' + sim_name + '_--_ni_' + it_str)
    matplotlib.pyplot.clf()

    # plot ne
    matplotlib.pyplot.plot(x.gridvalues, ne, linewidth = 2, color = 'blue')
    matplotlib.pyplot.grid()
    #    matplotlib.pyplot.axis([x.gridvalues[0], x.gridvalues[-1], -5.1, 5.1])
    matplotlib.pyplot.xlabel(r'position $x$', fontsize = 18)
    matplotlib.pyplot.ylabel(r'$n_i (t^n,x)$', fontsize = 18)
    matplotlib.pyplot.title(sim_name + r' ion density $n_i (t^n,x)$: $N_x$ = %d, $N_v$ = %d, $t^n$ = %2.3f, n = %04d' % (sim_params['total_dims'][0], sim_params['total_dims'][1],n*t.width, n))
    it_str = 'it%05d' % n
    matplotlib.pyplot.savefig('./plots/' + sim_name + '_--_ne_' + it_str)
    matplotlib.pyplot.clf()
