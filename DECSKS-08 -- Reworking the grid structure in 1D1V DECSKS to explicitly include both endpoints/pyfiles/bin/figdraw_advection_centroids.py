#!/usr/bin/env python
import matplotlib.pyplot as plt

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

prepoint = 2
width = 1
halfwidth = width / 2.
left_edge = prepoint - halfwidth
right_edge = prepoint + halfwidth
height = 8.

plt.bar(left_edge, height, width = width, color = [1,1,1], linestyle = 'dashed')
plt.text(prepoint, height + 0.5, r'$f(t^n, x_2)$', ha = 'center', fontsize = 16)
abscissa = range(11)
plt.xticks(abscissa)

postpoint = 7.25
postpoint_int = postpoint - (postpoint % 1)
left_edge = postpoint - halfwidth
right_edge = postpoint + halfwidth

# cell at i = 7
plt.vlines(7, 0, height, linestyles = 'dotted')

# cell at i = 8
plt.vlines(8, 0, height, linestyles = 'dotted')

# cell at i'' 
plt.vlines(7.25, 0, height, linestyles = 'dotted')

# shading overlap of cell between i'' and i = 7
plt.bar(7, height, width = 7.25 - 7, color = tableau20[13], linestyle = 'dotted')
plt.bar(7.25, height, width = 8 - 7.25, color = tableau20[9], linestyle = 'dotted')

# label fractional parts of remapped cells
plt.text(7.125, 4, r'$1 - \alpha$', ha = 'center', fontsize = 16)
plt.text((8 - 7.25) / 2 + 7.25,4, r'$\alpha$', ha = 'center', fontsize = 16)

# label remapped cells according to location at t + dt

plt.text(7, height + 0.5, r'$f(t^{n+1},x_7)$', ha = 'center', fontsize = 16)
plt.text(8, height + 0.5, r'$f(t^{n+1},x_8)$', ha = 'center', fontsize = 16)

arrow_length = 2
plt.arrow(3, 4, arrow_length, 0, width = 1, head_width = 2, head_length = 1, color = tableau20[3])
plt.text(3.75, 3.9, r'$\partial_t f + v\partial_x f = 0$', fontsize = 16)

plt.xlabel('grid indices, $i$', fontsize = 16)


plt.savefig('./../fig/advection_centroids.png')
 