import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-1,1,1000)
f =  np.sin(np.pi*x) + 1
nmax = 4 # maximum derivative desired to be plotted

fig, ax = plt.figure(), plt.subplot(111)

ax.plot(x,f, linewidth = 2, label = '$f(x)$')
ax.hold('on')

for n in range(1,nmax + 1):
    if np.mod(n,2) == 1: # n is odd
        dnf = (-1) ** (n + 1) * np.pi ** (2*n - 1) * np.cos(np.pi * x)
        dnf /= np.pi ** (2*n - 1) # i.e. we plot normalized derivatives

    else: # n is even
        dnf = (-1) ** n * np.pi ** (2*n) * np.sin(np.pi * x)
        dnf /= np.pi ** (2*n) # i.e. we plot normalized derivatives

    ax.plot(x,dnf, linewidth = 2, linestyle = '--', label = r'$\pi^{-%i}\partial_x^{%i} f$' % (n,n) )

ax.hold('off')
frame = ax.get_position() # position of plot center

# shrink frame of the plot to make room for a legend on the right side
ax.set_position([frame.x0, frame.y0, 
                 frame.width * 0.8, frame.height])

# Place legend to the right of the shrunken frame
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
          fancybox=True, shadow=True, ncol=1)

plt.show()
