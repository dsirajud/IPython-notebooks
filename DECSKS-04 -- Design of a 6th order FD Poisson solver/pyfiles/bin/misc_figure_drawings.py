#!/usr/bin/env python
import matplotlib.pyplot as plt


# numberline_without_ghost_cell.png

# set up the figure
fig = plt.figure()
fig.set_size_inches(8, 1)
ax = fig.add_subplot(111)
ax.set_xlim(0,10)
ax.set_ylim(0,15)

# draw lines

# left and right boundaries
xmin = 1
xmax = 9

# number of cells
N = 10 #
L = float(xmax - xmin) # (number line drawing range)
dx =  L / (N-1)

y = 6.5
height = 1
plt.hlines(y, xmin, xmax)

x = []
for i in range(N):
    x.append(xmin + i*dx)
    plt.vlines(x[i], y - height / 2., y + height / 2.)
    plt.text(x[i], y - 3*height, '%d' % i, horizontalalignment='center')

    if i == 0:
        plt.text(x[i], y + 2*height, r'$0 = x_{%d}$' % i, horizontalalignment='center')

    elif i == N-1:
        plt.text(x[i], y + 2*height, r'$x_{%d} = 1$' % i, horizontalalignment='center')

    else:
        plt.text(x[i], y + 2*height, r'$x_{%d}$' % i, horizontalalignment='center')

plt.axis('off')

plt.text(x[0] - 2*dx, y - 3*height, 'grid index', horizontalalignment = 'left')
plt.text(x[0] - 2*dx, y + 2*height, 'abscissa', horizontalalignment = 'left')
plt.text(x[N-1] + dx, y - 3*height, r'$N_x = 10$', horizontalalignment = 'left')

# px = 4
#plt.plot(x[px],y, 'ro', ms = 15, mfc = 'r')

#plt.annotate('Price five days ago', (x[px],y), xytext = (px - 1, y + 1),
#              arrowprops=dict(facecolor='black', shrink=0.1),
#              horizontalalignment='right')


plt.arrow(x[2], y + 5*height, dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.arrow(x[3], y + 5*height, -dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.text(x[2] + dx/2., y + 7.25*height, r'$\Delta x_R = \frac{L}{N_x - 1} = \frac{1}{9}$', horizontalalignment = 'center')

plt.arrow(x[0], y - 4.5*height, (N-1)*dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.arrow(x[N-1], y - 4.5*height, -(N-1)*dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.text(x[4] + dx/2., y - 7.5*height, r'$L = 1$', horizontalalignment = 'center')
plt.savefig('./../fig/numberline_without_ghost_cell.png')




# numberline_with_ghost_cell.png

# set up the figure
fig = plt.figure()
fig.set_size_inches(8, 1)
ax = fig.add_subplot(111)
ax.set_xlim(0,10)
ax.set_ylim(0,15)

# draw lines

# left and right boundaries
xmin = 1
xmax = 9

# number of cells
N = 10 #
L = float(xmax - xmin) # (number line drawing range)
dx =  L / N

y = 6.5
height = 1

x = []
for i in range(N):
    x.append(xmin + i*dx)
    plt.vlines(x[i], y - height / 2., y + height / 2.)
    plt.text(x[i], y - 3*height, '%d' % i, horizontalalignment='center')

    if i == 0:
        plt.text(x[i], y + 2*height, r'$0 = x_{%d}$' % i, horizontalalignment='center')

    elif i == N-1:
        plt.text(x[i], y + 2*height, r'$x_{%d} \neq 1$' % i, horizontalalignment='center')

    else:
        plt.text(x[i], y + 2*height, r'$x_{%d}$' % i, horizontalalignment='center')

plt.hlines(y, x[0], x[N-1])

plt.hlines(y,x[N-1],x[N-1] + dx, linestyle = 'dotted')
plt.vlines(x[N-1] + dx, y - height / 2., y + height/2.)

plt.axis('off')
plt.text(x[0] - 2*dx, y - 3*height, 'grid index', horizontalalignment = 'left')
plt.text(x[0] - 2*dx, y + 2*height, 'abscissa', horizontalalignment = 'left')
plt.text(x[N-1] + 1.5*dx, y - 3*height, r'$N_x = 10$', horizontalalignment = 'left')

# px = 4
#plt.plot(x[px],y, 'ro', ms = 15, mfc = 'r')

#plt.annotate('Price five days ago', (x[px],y), xytext = (px - 1, y + 1),
#              arrowprops=dict(facecolor='black', shrink=0.1),
#              horizontalalignment='right')


plt.arrow(x[2], y + 5*height, dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.arrow(x[3], y + 5*height, -dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.text(x[2] + dx/2., y + 7.25*height, r'$\Delta x_{NR} = \frac{L}{N_x} = \frac{1}{10}$', horizontalalignment = 'center')

plt.hlines(y,x[N-1],x[N-1] + dx, linestyle = 'dotted')
plt.vlines(x[N-1] + dx, y - height / 2., y + height/2.)

plt.text(x[N-1] + dx, y + 2*height, r'$x_{10}$', horizontalalignment='center')
plt.text(x[N-1] + dx , y - 3*height, r'10', horizontalalignment='center')

plt.arrow(x[0], y - 4.5*height, (N-1)*dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.arrow(x[N-1], y - 4.5*height, -(N-1)*dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.text(x[4] + dx/2., y - 7.5*height, r'$L = 1$', horizontalalignment = 'center')
plt.savefig('./../fig/numberline_with_ghost_cell.png')

# both figures


# numberline_with_redundancy.png

# set up the figure
fig = plt.figure()
fig.set_size_inches(8, 3)
ax = fig.add_subplot(111)
ax.set_xlim(0,10)
ax.set_ylim(0,45)

# draw lines

# left and right boundaries
xmin = 1
xmax = 9

# number of cells
N = 10 #
L = float(xmax - xmin) # (number line drawing range)
dx =  L / (N-1)

y = 35
height = 1
plt.hlines(y, xmin, xmax)

x = []
for i in range(N):
    x.append(xmin + i*dx)
    plt.vlines(x[i], y - height / 2., y + height / 2.)
    plt.text(x[i], y - 3*height, '%d' % i, horizontalalignment='center')

    if i == 0:
        plt.text(x[i], y + 2*height, r'$0 = x_{%d}$' % i, horizontalalignment='center')

    elif i == N-1:
        plt.text(x[i], y + 2*height, r'$x_{%d} = 1$' % i, horizontalalignment='center')

    else:
        plt.text(x[i], y + 2*height, r'$x_{%d}$' % i, horizontalalignment='center')

plt.axis('off')

plt.text(x[0] - 2*dx, y - 3*height, 'grid index', horizontalalignment = 'left')
plt.text(x[0] - 2*dx, y + 2*height, 'abscissa', horizontalalignment = 'left')
plt.text(x[N-1] + dx, y - 3*height, r'$N_x = 10$', horizontalalignment = 'left')

# px = 4
#plt.plot(x[px],y, 'ro', ms = 15, mfc = 'r')

#plt.annotate('Price five days ago', (x[px],y), xytext = (px - 1, y + 1),
#              arrowprops=dict(facecolor='black', shrink=0.1),
#              horizontalalignment='right')


plt.arrow(x[2], y + 5*height, dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.arrow(x[3], y + 5*height, -dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.text(x[2] + dx/2., y + 7.25*height, r'$\Delta x_R = \frac{L}{N_x - 1} = \frac{1}{9}$', horizontalalignment = 'center')

plt.arrow(x[0], y - 4.5*height, (N-1)*dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.arrow(x[N-1], y - 4.5*height, -(N-1)*dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.text(x[4] + dx/2., y - 7.5*height, r'$L = 1$', horizontalalignment = 'center')

plt.text(-1.5, 47, 'Grid without ghost cell', horizontalalignment = 'left', color = 'blue' )

# number of cells
dx =  L / N

y = 6.5
height = 1

x = []
for i in range(N):
    x.append(xmin + i*dx)
    plt.vlines(x[i], y - height / 2., y + height / 2.)
    plt.text(x[i], y - 3*height, '%d' % i, horizontalalignment='center')

    if i == 0:
        plt.text(x[i], y + 2*height, r'$0 = x_{%d}$' % i, horizontalalignment='center')

    elif i == N-1:
        plt.text(x[i], y + 2*height, r'$x_{%d} \neq 1$' % i, horizontalalignment='center')

    else:
        plt.text(x[i], y + 2*height, r'$x_{%d}$' % i, horizontalalignment='center')

plt.hlines(y, x[0], x[N-1])

plt.hlines(y,x[N-1],x[N-1] + dx, linestyle = 'dotted')
plt.vlines(x[N-1] + dx, y - height / 2., y + height/2.)

plt.axis('off')
plt.text(x[0] - 2*dx, y - 3*height, 'grid index', horizontalalignment = 'left')
plt.text(x[0] - 2*dx, y + 2*height, 'abscissa', horizontalalignment = 'left')
plt.text(x[N-1] + 2*dx, y - 3*height, r'$N_x = 10$', horizontalalignment = 'left')

# px = 4
#plt.plot(x[px],y, 'ro', ms = 15, mfc = 'r')

#plt.annotate('Price five days ago', (x[px],y), xytext = (px - 1, y + 1),
#              arrowprops=dict(facecolor='black', shrink=0.1),
#              horizontalalignment='right')


plt.arrow(x[2], y + 5*height, dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.arrow(x[3], y + 5*height, -dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.text(x[2] + dx/2., y + 7.25*height, r'$\Delta x_{NR} = \frac{L}{N_x} = \frac{1}{10}$', horizontalalignment = 'center')

plt.arrow(x[0], y - 4.5*height, (N-1)*dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.arrow(x[N-1], y - 4.5*height, -(N-1)*dx, 0, head_width = 1, head_length = 0.1, fc = 'k', ec = 'k', shape = 'full')
plt.text(x[4] + dx/2., y - 7.5*height, r'$L = 1$', horizontalalignment = 'center')

plt.text(x[N-1] + dx, y + 2*height, r'$x_{10}$', horizontalalignment='center')
plt.text(x[N-1] + dx , y - 3*height, r'10', horizontalalignment='center')

plt.text(-1.5, 20, 'Grid with ghost cell', horizontalalignment = 'left', color = 'blue' )

plt.vlines(x[6], -15, 50, linestyle = 'dotted')
plt.text(x[6] + .1, 45, 'gridpoints do not overlap', horizontalalignment = 'left')

plt.savefig('./../fig/both_grids.png')
