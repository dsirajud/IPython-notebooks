#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# set up the figure
fig = plt.figure()
fig.set_size_inches(14,5)
ax = fig.add_subplot(111)
ax.set_xlim(0,10)
ax.set_ylim(0,10)

# color schemes
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

# ellipses module

el = Ellipse((2, -1), 0.5, 0.5)
ax.add_patch(el)


prepoint = 2
width = 1
halfwidth = width / 2.
left_edge = prepoint - halfwidth
right_edge = prepoint + halfwidth
height = 8.

plt.bar(left_edge, height, width = width, color = tableau20[1])
plt.text(prepoint, height + 0.5, r'$f(t^n, x_2)$', ha = 'center', fontsize = 16)
abscissa = range(11)
plt.xticks(abscissa)



postpoint = 7.25
postpoint_int = postpoint - (postpoint % 1)
left_edge = postpoint - halfwidth
right_edge = postpoint + halfwidth

# cell at i = 7
plt.bar(postpoint_int - halfwidth, height, width = width, color = [1,1,1], linestyle = 'dashed')

# cell at i = 8
plt.bar(postpoint_int + 1 - halfwidth, height, width = width, color = [1,1,1], linestyle = 'dashed')


# advected cell at x + vt
plt.bar(left_edge, height, width = width, color = tableau20[1])
plt.text(postpoint, height + 0.5, r'$f(t^n, x_2 + v\Delta t)$', ha = 'center', fontsize = 16)

# shade overlap fractions on cells at i = 7 and i = 8

# shading overlap of cell i = 7
plt.bar(left_edge, height, width = (7 + halfwidth) - left_edge, color = tableau20[13], linewidth = 0)
plt.bar(7 + halfwidth, height, width = right_edge - (7 + halfwidth), color = tableau20[9], linewidth = 0)

# label fractional parts of advected cell
plt.text(left_edge + (7 + halfwidth - left_edge) / 2., 4, r'$\zeta$', ha = 'center', fontsize = 16)
plt.text((7 + halfwidth) + ( right_edge - 7 - halfwidth) / 2., 4, r'$1 - \zeta$', ha = 'center', fontsize = 16)

# retrace cell i = 7 right-hand boundary (above overwrites it)
plt.vlines(postpoint_int + halfwidth, 0, height, linestyles = 'dashed')


arrow_length = 2
plt.arrow(3, 4, arrow_length, 0, width = 1, head_width = 2, head_length = 1, color = tableau20[3])
plt.text(3.75, 3.9, r'$\partial_t f + v\partial_x f = 0$', fontsize = 16)

plt.xlabel('grid indices, $i$', fontsize = 16)



# left and right boundaries
xmin = 1
xmax = 9

plt.savefig('./../fig/advection.png')
