import math
import numpy
import pylab


def func(x):
    # get midpoint value and store
    x_crit = 0.315128259098501215428873837159675538342649935961030535963
    if x < x_crit:
        return 0.2*numpy.sin(2*numpy.pi * x / 1.) + .25
    else:
        b = (0.2*numpy.sin(2*numpy.pi * x_crit / 1.) + .25) + x_crit * .5 # last term is (x - x0)*slope, where x - x0 is the distance to x = 0
        return  b - .5 *x


    #def func(x):
    #    # get midpoint value and store
    #    if x < 0.5:
    #        return 0.2*numpy.sin(2*numpy.pi * x / 1.) + .25
    #    else:
    #        return .25 + .25/2 - .25 *x

def simplegrid():

    xmin = 0.0
    xmax = 1.0

    Nx = 50
    x = numpy.linspace(xmin, xmax, Nx)

    f = numpy.zeros_like(x)
    for i in range(Nx):
        f[i] = func(x[i])

    pylab.plot(x,f, lw = 2, color = 'cornflowerblue')
    nzones = 7


    dx = (xmax - xmin)/float(nzones)

    # plot DBC boundary conditions
    LDBC = f[0]
    RDBC = f[-1]
    pylab.plot([xmin-3/2.*dx,xmin], [LDBC,LDBC], color="b", lw=2)
    pylab.plot([xmax,xmax + 3/2.*dx], [RDBC,RDBC], color="b", lw=2)


    xl = numpy.arange(nzones)*dx
    xr = (numpy.arange(nzones)+1)*dx

    xc = 0.5*(xl + xr)

    pylab.plot([xmin,xmax], [0,0], color="k", lw=2)
    pylab.plot([xmin - 2*dx,xmin], [0,0], color="k", linestyle = '--', lw=2)
    pylab.plot([xmax, xmax + 2*dx], [0,0], color="k", linestyle = '--', lw=2)

    # plot ghost points and ghost boundaries

    # left ghost point
    pylab.plot(xc[0] - dx, 0, 'D', markeredgecolor = 'k', markerfacecolor = 'white')

    # left ghost boundary
    pylab.plot([xl[0] - dx, xl[0] - dx], [-0.05, 0.05], color="g")

    # left ghost boundary label
    pylab.text(xl[0] - dx, -.075, r"$-3/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'g')


    # right ghost point
    pylab.plot(xc[-1] + dx, 0, 'D', markeredgecolor = 'k', markerfacecolor = 'white')

    # right ghost boundary
    pylab.plot([xr[-1] + dx, xr[-1] + dx], [-0.05, 0.05], color="g")

    # left ghost boundary label
    pylab.text(xr[-1] + dx, -.075, r"$N+1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'g')




    # left wall
    pylab.plot([xl[0], xl[0]], [0.0, 0.42], color = 'k', lw = 2)


    # right wall
    pylab.plot([xr[-1], xr[-1]], [0, 0.42], color = 'k', lw = 2)

    n = 0
    while (n < nzones):


        # draw center marker
        pylab.plot(xc[n], 0, 'D', color="k")

        # draw edge marker
        if (n == 0):
            pylab.plot([xl[0], xl[0]], [-0.05, 0.05], color="g")

        pylab.plot([xr[n], xr[n]], [-0.05, 0.05], color="g")

        n += 1

    # left boundary
    pylab.text(xl[0], -0.075, r"$-1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'g')

    # left boundary
    pylab.text(xl[0], 0.55, r"$\phi (x_{-1/2}) = g_{\ell}$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium", color = 'b')


    # left ghost point
    pylab.text(xc[0] - dx, -0.1, r"$-1$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'k')

    # right ghost point
    pylab.text(xc[-1] + dx, -0.1, r"$N$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'k')

    # left boundary
    pylab.text(xl[0], -0.155, r"$x = a$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium", color = 'red')


    # right boundary
    pylab.text(xr[-1], 0.55, r"$\phi (x_{(N-1)+1/2}) = g_{r}$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium", color = 'b')


    # right boundary
    pylab.text(xr[-1], -0.2, r"$(N-1) + 1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'g')


    # right boundary
    pylab.text(xr[-1], -0.125, r"$x = b$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium", color = 'red')


    # label a few edges
    pylab.text(xl[1], -0.075, r"$1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'g')

    pylab.text(xl[2], -0.075, r"$3/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'g')

    pylab.text(xl[3], -0.075, r"$i-1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'g')

    pylab.text(xl[4], -0.075, r"$i + 1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'g')

    pylab.text(xl[5], -0.075, r"$\ldots$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'g')

    pylab.text(xl[6], -0.2, r"$(N-1) - 1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'g')


    # label a few cell-centers at the center, general i
    pylab.text(xc[0], -0.1, r"$0$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    pylab.text(xc[1], -0.1, r"$1$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")


    pylab.text(xc[2], -0.1, r"$\ldots$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    pylab.text(xc[3], -0.1, r"$i$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    pylab.text(xc[4], -0.1, r"$\ldots$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    pylab.text(xc[5], -0.1, r"$N-2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    pylab.text(xc[6], -0.1, r"$N-1$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")


    # plot cell values, note that x is measured in terms of data units
    # whereas y is measured in relative units (0 to 1)
    # where 0 = bottom, 0.5 = middle, 1 = top
    #
    # for example, at x = 0, we have the function f = func(x[0])
    # if we want to plot a rectangle from the data values
    # y = 0 to y = func(x[0]), where the axis is spans
    #
    #         y in [-.25, 0.6], length = .6 - (-.25) = .85
    #
    # note that y = 0 is .25 units from the bottom, and
    # .25 / .85 relative units from the bottom
    #
    # note that y = .25 is .5 units from the bottom, or
    # 0.5 / .85 units of the whole axis
    #
    # in general, we at each f(x), we have the relative
    # y height to be (distance from the bottom to the origin) + f(x)
    # i.e. at f[xc[i]], the y value in axvspan should be
    # ymax = (f[xc[i]] + .25) / .85

    pylab.axis([xmin-3/2.*dx,xmax+3/2.*dx, -0.25, 0.6])
    pylab.axis("off")

    # plot rectangles at all function values at each center xc[i]
    for i in range(len(xc)):
        pylab.axvspan(xl[i], xl[i] + dx, 0.25 / .85, (func(xc[i]) + .25) / .85, facecolor='b', alpha=0.2)

    # boundary ghost cells
    pylab.axvspan(xl[0] - dx, xl[0], 0.25 / .85, (func(xl[0]) + .25) / .85, facecolor='cyan', alpha=0.2, hatch = '//')
    pylab.axvspan(xr[-1], xr[-1] + dx, 0.25 / .85, (func(xr[-1]) + .25) / .85, facecolor='cyan', alpha=0.2, hatch = '\\')

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(9.0,3.0)

    pylab.savefig("cc_grid_DBCs_with_cell_values.png")

if __name__== "__main__":
    simplegrid()
