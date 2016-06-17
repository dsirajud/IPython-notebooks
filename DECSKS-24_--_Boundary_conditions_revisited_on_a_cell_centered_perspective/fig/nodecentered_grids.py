import math
import numpy
import pylab

def simplegrid():

    xmin = 0.0
    xmax = 1.0

    nzones = 7

    dx = (xmax - xmin)/float(nzones)

    xl = numpy.arange(nzones)*dx
    xr = (numpy.arange(nzones)+1)*dx

    xc = 0.5*(xl + xr)

    pylab.plot([xmin - .5*dx,xmax + .5*dx], [0,0], color="k", lw=2)

    n = 0
    while (n < nzones):


        # draw center marker
        pylab.plot([xc[n], xc[n]], [-0.05, 0.05], color="g")

        # draw edge marker
        if (n == 0):
            pylab.plot(xl[0], 0, 'D', color="k")

        pylab.plot(xr[n], 0, 'D', color="k")

        n += 1

    # draw outermost edges
    pylab.plot([xl[0] - .5*dx, xl[0] - .5*dx], [-0.05, 0.05], lw = 2, color="g")

    # draw outermost edges
    pylab.plot([xr[-1] + .5*dx, xr[-1] + .5*dx], [-0.05, 0.05], lw = 2, color="g")

    # label outermost edges
    pylab.text(xl[0] - .5*dx, -0.1, r"$-1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'green')


    # label outermost edges
    pylab.text(xr[-1] + .5*dx, -0.1, r"$(N-1)$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'green')

    pylab.text(xr[-1] + .5*dx, -0.175, r"$+1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'green')


    # left boundary
    pylab.text(xl[0], -0.075, r"$0$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    # left boundary
    pylab.text(xl[0], -0.155, r"$x = a$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium", color = 'red')

    # right boundary
    pylab.text(xr[-1], -0.075, r"$N-1$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    # right boundary
    pylab.text(xr[-1], -0.125, r"$x = b$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium", color = 'red')


    # label a few edges
    pylab.text(xl[1], -0.075, r"$1$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    pylab.text(xl[2], -0.075, r"$2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    pylab.text(xl[3], -0.075, r"$\ldots$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    pylab.text(xl[4], -0.075, r"$i$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    pylab.text(xl[5], -0.075, r"$\ldots$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    pylab.text(xl[6], -0.075, r"$N-2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    pylab.text(xr[6], -0.075, r"$N-1$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")


    # label a few cell-centers at the center, general i
    pylab.text(xc[0], -0.1, r"$1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'green')

    pylab.text(xc[1], -0.1, r"$3/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'green')


    pylab.text(xc[3], -0.1, r"$i - 1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'green')

    pylab.text(xc[4], -0.1, r"$i+1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'green')

    pylab.text(xr[-1] - .5*dx, -0.1, r"$(N-1)$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'green')

    pylab.text(xr[-1] - .5*dx, -0.175, r"$-1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'green')




    pylab.axis([xmin-0.5*dx,xmax+0.5*dx, -0.25, 0.25])
    pylab.axis("off")

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(9.0,2.0)

    pylab.savefig("node_centered_grid.png")
    pylab.clf()
def simplegrid2():

    xmin = 0.0
    xmax = 1.0

    nzones = 7

    dx = (xmax - xmin)/float(nzones)

    xl = numpy.arange(nzones)*dx
    xr = (numpy.arange(nzones)+1)*dx

    xc = 0.5*(xl + xr)

    pylab.plot([xmin,xmax], [0,0], color="k", lw=2)
    pylab.plot([xmin-.5*dx,xmin], [0,0], linestyle = '--', color="k", lw=2)
    pylab.plot([xmax,xmax+.5*dx], [0,0], linestyle = '--',color="k", lw=2)

    n = 0
    while (n < nzones):

        # draw center marker
        pylab.plot([xl[n], xl[n]], [-0.05, 0.05], color="purple")

        # draw edge marker
        if (n == 0):
            pylab.plot(xl[0], 0, 'D', color="k")
            pylab.plot([xl[0], xl[0]], [0, 0.15], color="purple", lw = 2)
            pylab.text(xl[0], 0.22, r"$\mathrm{Periodic \, boundary}$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium", color = 'purple')

        pylab.plot(xr[n], 0, 'D', color="k")

        n += 1

    pylab.plot([xl[n-1] + dx, xl[n-1] + dx], [-0.05, 0.05], color="purple")
    pylab.plot([xl[n-1] + dx, xl[n-1] + dx], [0, 0.15], color="purple", lw= 2)
    pylab.text(xl[n-1] + dx, 0.22, r"$\mathrm{Periodic \, boundary}$",
       horizontalalignment='center', verticalalignment='top',
       fontsize="medium", color = 'purple')


    # left boundary
    pylab.text(xl[0], -0.075, r"$0$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")

    # left boundary
    pylab.text(xl[0], -0.155, r"$x = a$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium", color = 'red')

    # right boundary
    pylab.text(xr[-1], -0.075, r"$N-1$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")

    # right boundary
    pylab.text(xr[-1], -0.125, r"$x = b$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium", color = 'red')

    pylab.text(xr[-1], -0.165, r"$f_{N-1} = f_0$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium", color = 'blue')

    # label a few edges
    pylab.text(xl[1], -0.075, r"$1$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")

    pylab.text(xl[2], -0.075, r"$2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")

    pylab.text(xl[2], -0.175, r"$k$", color = 'blue',
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")

    pylab.plot([xl[0], xl[2]], [.125, .125], lw = 2, linestyle= ':',color = 'k')
    pylab.plot(xl[2], .125, marker = '>',
               lw = 2, linestyle= ':',color = 'k')


    pylab.text(xl[3], -0.075, r"$\ldots$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")

    pylab.text(xl[4], -0.075, r"$i$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")

    pylab.plot([xl[4], xr[6]], [.125, .125], lw = 2, linestyle= ':',color = 'k')
    pylab.plot(xr[5], .125, marker = '>',
               lw = 2, linestyle= ':',color = 'k')

    pylab.plot(xl[4], .125, marker = 'x',
               lw = 2, linestyle= ':',color = 'k')

    pylab.text(xl[5], -0.075, r"$\ldots$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")

    pylab.text(xl[6], -0.075, r"$N-2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")

    pylab.text(xr[6], -0.075, r"$N-1$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")



    pylab.axis([xmin-0.5*dx,xmax+0.5*dx, -0.25, 0.25])
    pylab.axis("off")

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(9.0,2.0)

    pylab.savefig("node_centered_grid_2.png")


if __name__== "__main__":
    simplegrid()
    simplegrid2()
