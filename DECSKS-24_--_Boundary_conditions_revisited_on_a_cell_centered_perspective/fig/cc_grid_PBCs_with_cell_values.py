import math
import numpy
import pylab


def func(x):
    # get midpoint value and store, the function will be piecewise but continuous at his critical point
    return 0.2*numpy.sin(2*numpy.pi * x / 1.) + .25

    #def func(x):
    #    # get midpoint value and store
    #    if x < 0.5:
    #        return 0.2*numpy.sin(2*numpy.pi * x / 1.) + .25
    #    else:
    #        return .25 + .25/2 - .25 *x

def simplegrid():

    xmin = 0.0
    xmax = 1.0
    nzones = 7

    Nx = 50
    dx = (xmax - xmin) / float(nzones)
    x = numpy.linspace(xmin-dx, xmax+dx, Nx)

    f = numpy.zeros_like(x)
    for i in range(Nx):
        f[i] = func(x[i])

    pylab.plot(x,f, lw = 2, color = 'cornflowerblue')

    dx = (xmax - xmin)/float(nzones)

    # plot DBC boundary conditions
    LDBC = f[0]
    RDBC = f[-1]
    #    pylab.plot([xmin-3/2.*dx,xmin], [LDBC,LDBC], color="b", lw=2)
    #    pylab.plot([xmax,xmax + 3/2.*dx], [RDBC,RDBC], color="b", lw=2)


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

    # plot rectangles at the ghost values, left
    pylab.axvspan(xl[0] - dx, xl[0], 0.25 / .85, (func(xc[-1]) + .25) / .85, facecolor='purple', alpha=0.2)
    # plot rectangles at the ghost values, right
    pylab.axvspan(xr[-1], xr[-1] + dx, 0.25 / .85, (func(xc[0]) + .25) / .85, facecolor='purple', alpha=0.2)

    # line connecting periodic pairs, i = 0 and i = N
    pylab.plot([xc[0], xc[-1]+dx], [func(xc[0]), func(xc[0])], color = 'k', linestyle = ':', lw = 2)
    pylab.plot(xc[0], func(xc[0]), marker = '<', color = 'k', lw = 2)
    pylab.plot(xc[-1]+dx, func(xc[0]), marker = '>', color = 'k', lw = 2)
    # text annotation of BC
    pylab.text((xc[-1] +dx - xc[0]) / 2. + .2, func(xc[0]) + .08, r"$\mathrm{Periodic\, BC:\,\, } \phi_{0} = \phi_{N}$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")


    # line connecting periodic pairs, i = -1 and i = N-1
    pylab.plot([xc[0]-dx, xc[-1]], [func(xc[-1]), func(xc[-1])], color = 'k', linestyle = ':', lw = 2)
    pylab.plot(xc[0]-dx, func(xc[-1]), marker = '<', color = 'k', lw = 2)
    pylab.plot(xc[-1], func(xc[-1]), marker = '>', color = 'k', lw = 2)
    # text annotation of BC
    pylab.text((xc[-1] - xc[0] - dx) / 2., func(xc[-1]) + .08, r"$\mathrm{Periodic \,BC:\,\, } \phi_{-1} = \phi_{N-1}$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")
    


    # boundary ghost cells

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(9.0,3.0)

    pylab.savefig("cc_grid_PBCs_with_cell_values.png")
    pylab.clf()

if __name__== "__main__":
    simplegrid()
