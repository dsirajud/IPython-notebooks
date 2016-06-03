import math
import numpy
import pylab

def simplegrid():

    xmin = 0.0
    xmax = 2.0

    nzones = 7

    dx = (xmax - xmin)/float(nzones)

    xl = numpy.arange(nzones)*dx
    xr = (numpy.arange(nzones)+1)*dx

    xc = 0.5*(xl + xr)

    pylab.plot([xmin,xmax], [0,0], color="k", lw=2)
    pylab.plot([xmin - 2*dx,xmin], [0,0], color="k", linestyle = '--', lw=2)
    pylab.plot([xmax, xmax + 2*dx], [0,0], color="k", linestyle = '--', lw=2)

    # plot ghost points and ghost boundaries

    # left ghost point
    pylab.plot(xc[0] - dx, 0, 'D', markeredgecolor = 'k', markerfacecolor = 'white')

    # left ghost boundary, vertical
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
    pylab.plot([xl[0], xl[0]], [0, 0.85], color = 'k', lw = 2)

    # right wall
    pylab.plot([xr[-1], xr[-1]], [0, 0.85], color = 'k', lw = 2)

    # label regions

    # left absorber region
    pylab.text(xl[0] - dx, .85, r"$\mathrm{absorber}$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="large", color = 'b')


    # right absorber region
    pylab.text(xr[-1] + dx, .85, r"$\mathrm{absorber}$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="large", color = 'b')




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
    pylab.text(xl[0], 1., r"$x = a$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium", color = 'red')


    # right boundary
    pylab.text(xr[-1], -0.2, r"$(N-1) + 1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'g')


    # right boundary
    pylab.text(xr[-1], 1.05, r"$x = b$",
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

    pylab.text(xc[4], -0.1, r"$i+1$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    pylab.text(xc[5], -0.1, r"$N-2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")

    pylab.text(xc[6], -0.1, r"$N-1$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")



    # label density packet at prepoint
    pylab.text(xc[3] + .2*dx, .5, r"$f^n_{i,j}$",
               color = 'k',
               horizontalalignment='center', verticalalignment='top',
               fontsize="large")


    # label on characteristic from i to (-1, -1/2)
    pylab.text((xc[3] - xc[0] + .2*dx) / 2 + .1*dx, .45, r"$\mathcal{C}_i = S^n_i + U^n_i$",
               color = 'purple',
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")

    # line draw, characteristic from i to (-1, -1/2)
    pylab.plot([xc[3], xc[0]- .20*dx], [.25, .25], color="purple")
    # start marker
    pylab.plot(xc[3], 0.5, '|', markeredgecolor = 'magenta', markerfacecolor = 'm', markersize = 16)

    # line draw integral part of characteristic from i to i = 0, S
    pylab.plot([xc[3], xc[0]], [.5, .5], color="magenta")
    # end marker
    pylab.plot(xc[0]+.05*dx, 0.5, '<', markeredgecolor = 'magenta', markerfacecolor = 'm', markersize = 6)
    # label on line
    pylab.text((xc[3] - xc[0] + .2*dx) / 2 + .1*dx, .7, r"$S_i$",
               color = 'm',
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")


    # line draw fractional part of characteristic from i to i = 0, U
    pylab.plot([xc[0], xc[0] - .20*dx], [.5, .5], color="magenta")
    # end marker, <
    pylab.plot(xc[0] - .20*dx, 0.5, '<', markeredgecolor = 'magenta', markerfacecolor = 'm', markersize = 6)
    # label on line
    pylab.text(xc[0]- .1*dx, .73, r"$-U^n_i$",
               color = 'm',
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")

    # line draw (1-U) from true postpoint to k2 = -1
    pylab.plot([xc[0] - .2*dx, xc[0] -dx], [.5, .5], color="magenta")
    # start marker, >
    pylab.plot(xc[0] - .32*dx, .5, '>',color="magenta", markersize = 6, markeredgecolor = 'm')
    # end marker, <
    pylab.plot(xc[0] -.93*dx, .5, '<',color="magenta", markersize = 6, markeredgecolor = 'm')
    # label 1 - alpha
    pylab.text(xc[0] - .75*dx, .73, r"$1 + U^n_i$",
               color = 'm',
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")


    pylab.text(xc[0] - .27*dx, -.225, r"$k$",
               color = 'm',
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")


    # prepoint marker
    pylab.plot(xc[3], 0.25, '|', markeredgecolor = 'purple', markerfacecolor = 'purple', markersize = 16)
    # postpoint marker
    pylab.plot(xc[0] - .20*dx, 0.25, '<', markeredgecolor = 'purple', markerfacecolor = 'purple', markersize = 6)

    # vert. line at postpoint exact, magenta portion
    pylab.plot([xc[0] - .27*dx, xc[0] - .27*dx], [.35, .6], color="magenta")
    # vert. line at postpoint, purple portion
    pylab.plot([xc[0] - .27*dx, xc[0] - .27*dx], [0, .35], color="purple")

    # vert. line at k1 = 0, magenta portion
    pylab.plot([xc[0], xc[0]], [.35, .6], color="magenta")
    # vert. line at k1 = 0, purple portion
    pylab.plot([xc[0], xc[0]], [0, .35], color="purple")
    # label at k1
    pylab.text(xc[0] - dx, -.225, r"$k_2$",
               color = 'm',
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")

    # vert. line at k2 = -1
    pylab.plot([xc[0] - dx, xc[0] - dx], [0.05, .6], color="purple")
    # label at k2
    pylab.text(xc[0], -.225, r"$k_1$",
               color = 'm',
               horizontalalignment='center', verticalalignment='top',
               fontsize="medium")
    

    pylab.axis([xmin-3/2.*dx,xmax+3/2.*dx, -0.35, 1.2])
    pylab.axis("off")

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(10.0,2.0)

    pylab.savefig("remap_rule_absorbing.png")

if __name__== "__main__":
    simplegrid()
