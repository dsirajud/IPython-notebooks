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
#plt.bar(postpoint_int - halfwidth, height, width = width, color = [1,1,1], linestyle = 'dashed')
plt.vlines(6.5, 0.75*height, height, linestyles = 'dashed')
plt.hlines(height, 6.5, 7.5, linestyles = 'dashed')

# cell at i = 8
#plt.bar(postpoint_int + 1 - halfwidth, height, width = width, color = [1,1,1], linestyle = 'dashed')

plt.vlines(8.5, .25*height, height, linestyles = 'dashed')
plt.hlines(height, 8.5, 7.5, linestyles = 'dashed')

# cell interface line at 7.5
plt.vlines(7.5, .75*height, height, linestyles = 'dashed')

# shade overlap fractions on cells at i = 7 and i = 8

# shading overlap of cell i = 7
plt.bar(7 - halfwidth, .75*height, width = width, color = tableau20[13])
plt.bar(7 + halfwidth, .25*height, width = width, color = tableau20[9])

# label fractional parts of remapped cells
plt.text(7, 4, r'$1 - \alpha$', ha = 'center', fontsize = 16)
plt.text(8,1, r'$\alpha$', ha = 'center', fontsize = 16)

# label remapped cells according to location at t + dt

plt.text(7, height + 0.5, r'$f(t^{n+1},x_7)$', ha = 'center', fontsize = 16)
plt.text(8, height + 0.5, r'$f(t^{n+1},x_8)$', ha = 'center', fontsize = 16)

arrow_length = 2
plt.arrow(3, 4, arrow_length, 0, width = 1, head_width = 2, head_length = 1, color = tableau20[3])
plt.text(3.75, 3.9, r'$\partial_t f + v\partial_x f = 0$', fontsize = 16)

plt.xlabel('grid indices, $i$', fontsize = 16)



# left and right boundaries
xmin = 1
xmax = 9

plt.savefig('./../fig/remap_centroids.png')
