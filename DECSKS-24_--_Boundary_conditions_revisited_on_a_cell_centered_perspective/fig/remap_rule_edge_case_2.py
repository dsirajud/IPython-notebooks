import math
import numpy
import pylab

def simplegrid():

    xmin = -5.0
    xmax = 5.0

    nzones = 2

    dx = (xmax - xmin)/float(nzones)

    xl = numpy.arange(nzones)*dx
    xr = (numpy.arange(nzones)+1)*dx

    xc = 0.5*(xl + xr)

    pylab.plot([xmin+2,xmax-2], [0,0], color="k", lw=2)
    pylab.plot([xmin + 2 - .125*dx,xmin + 2], [0,0], color="k", linestyle = ':', lw=2)
    pylab.plot([xmax - 2, xmax + .125*dx - 2], [0,0], color="k", linestyle = ':', lw=2)

    # plot ghost points and ghost boundaries

    # left ghost point
    pylab.plot(xc[0] - dx, 0, 'D', markeredgecolor = 'k', markerfacecolor = 'white')


    # left wall
    pylab.plot([xl[0], xl[0]], [0, 0.85], color = 'k', lw = 2)

    # draw center marker
    pylab.plot(xc[0], 0, 'D', color="k")

    # draw edge marker

    # left boundary
    pylab.text(xl[0], -0.075, r"$-1/2$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'g')

    # left ghost point
    pylab.text(xc[0] - dx, -0.1, r"$-1$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small", color = 'k')

    # left boundary
    pylab.text(xl[0], .95, r"$x = a$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="large", color = 'red')

    # label a few cell-centers at the center, general i
    pylab.text(xc[0], -0.1, r"$0$",
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")


    # line draw fractional part of characteristic from i = 0 to k_true, magenta
    pylab.plot([xc[0], xc[0] - .77*dx], [.5, .5], color="magenta")
    # start marker, <
    pylab.plot(xc[0] - .01*dx, 0.5,'>', color = 'm', markeredgecolor = 'magenta', markerfacecolor = 'm', markersize = 8)
    # end marker, <
    pylab.plot(xc[0] - .77*dx + .01*dx, 0.5, '<', color = 'm',markeredgecolor = 'magenta', markerfacecolor = 'm', markersize = 8)
    # label on line
    pylab.text(xc[0]- dx/2. + .12*dx, .6, r"$-U^n_i$",
               color = 'm',
               horizontalalignment='center', verticalalignment='top',
               fontsize="large")

    # line draw fractional part of characteristic from i = 0 to k_true, blue
    pylab.plot([xc[0]-dx/2., xc[0] - .77*dx], [.25, .25], color="blue")
    # start marker, <
    pylab.plot(xc[0] - dx/2. - .01*dx, 0.25,'>', color = 'blue', markeredgecolor = 'blue', markerfacecolor = 'blue', markersize = 8)
    # end marker, <
    pylab.plot(xc[0] - .77*dx + .01*dx, 0.25, '<', color = 'm',markeredgecolor = 'blue', markerfacecolor = 'blue', markersize = 8)
    # label on line
    pylab.text(xc[0]- dx + .12*dx, .35, r"$1+U^n_i$",
               color = 'blue',
               horizontalalignment='center', verticalalignment='top',
               fontsize="large")



    # line draw (1-U) from true postpoint to k2 = -1
    pylab.plot([xc[0] - .77*dx, xc[0] -dx], [.5, .5], color="magenta")
    # start marker, >
    pylab.plot(xc[0] - .77*dx - .01*dx, .5, '>',color="magenta", markersize = 8, markeredgecolor = 'm', markerfacecolor='m')
    # end marker, <
    pylab.plot(xc[0] - dx + .01*dx, .5, '<',color="magenta", markersize = 8, markeredgecolor = 'm', markerfacecolor='m')
    # label 1 - alpha
    pylab.text(xc[0] - dx + .12*dx, .6, r"$1 + U^n_i$",
               color = 'm',
               horizontalalignment='center', verticalalignment='top',
               fontsize="large")


    # line draw, from i = 0 to i = -1/2, total length = 1/2
    pylab.plot([xc[0], xc[0] - dx/2.], [.25, .25], '--', lw = 2, color="k")
    # <
    pylab.plot(xc[0] - dx/2. + .01*dx, .25, '<', lw = 2, color="k", markersize = 8)
    # >
    pylab.plot(xc[0] - .01*dx, .25, '>', lw = 2, color="k", markersize = 8)
    # label = 1/2
    pylab.text(xc[0] - 1/4.*dx,  .35, r"$1/2$",
               color = 'k',
               horizontalalignment='center', verticalalignment='top',
               fontsize="large")



    # line draw (1-U) from true postpoint to k2 = -1/2
    pylab.plot([xc[0]- .77*dx, xc[0] - dx], [.25, .25], color="blue")
    # start marker, >
    pylab.plot(xc[0] - .77*dx - .01*dx, .25, '>',color="blue", markersize = 8, markeredgecolor = 'blue', markerfacecolor='blue')
    # end marker, <
    pylab.plot(xc[0] - dx + .01*dx, .25, '<',color="blue", markersize = 8, markeredgecolor = 'blue', markerfacecolor='blue')
    # label 1 - alpha
    pylab.text(xc[0] - 1/2.*dx +  ((xc[0] - .77*dx) - (xc[0]-1/2.*dx))/2., .35, r"$-U_i^n - 1/2$",
               color = 'blue',
               horizontalalignment='center', verticalalignment='top',
               fontsize="large")

    # vert. line at true postpoint, purple portion
    pylab.plot([xc[0]-.77*dx, xc[0] - .77*dx], [0, .6], color="purple")

    # vert. line at k1 = 0, purple portion
    pylab.plot([xc[0], xc[0]], [0, .6], color="purple")
    # label at k1
    pylab.text(xc[0] - dx, -.16, r"$k_2, $",
               color = 'm',
               horizontalalignment='center', verticalalignment='top',
               fontsize="large")

    # true postpoint k
    pylab.text(xc[0] - .77*dx, -.16, r"$k$",
               color = 'purple',
               horizontalalignment='center', verticalalignment='top',
               fontsize="large")
    

    # vert. line at k2 = -1
    pylab.plot([xc[0] - dx, xc[0] - dx], [0.0, .6], color="purple")
    # label at k2
    pylab.text(xc[0], -.16, r"$k_1$",
               color = 'm',
               horizontalalignment='center', verticalalignment='top',
               fontsize="large")

    # label at k1 blue
    pylab.text(xc[0] - dx/2., -.16, r"$k_1$",
               color = 'blue',
               horizontalalignment='center', verticalalignment='top',
               fontsize="large")

    # label at k2 blue
    pylab.text(xc[0] - dx, -.16, r"$\qquad k_2$",
               color = 'blue',
               horizontalalignment='center', verticalalignment='top',
               fontsize="large")



    pylab.axis([xmin+.45*dx,xmax-.45*dx, -0.35, 1.2])
    pylab.axis("off")

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()

    pylab.savefig("remap_rule_absorbing_edge_case_2.png")
if __name__== "__main__":
    simplegrid()
