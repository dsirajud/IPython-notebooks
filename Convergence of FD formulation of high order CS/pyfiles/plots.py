import numpy as np
import matplotlib.pyplot as plt
import convergence_routines

ax = -0.5
bx = 1.5
L = bx - ax
Nx = 1000
dx = L / Nx

x = np.zeros(Nx)

for i in range(Nx):
    x[i] = ax + i*dx

f =  convergence_routines.function(x)

df1 = convergence_routines.df1(x)
df2 = convergence_routines.df2(x)
df3 = convergence_routines.df3(x)
df4 = convergence_routines.df4(x)

fig, ax = plt.figure(), plt.subplot(111)

ax.plot(x,f, linewidth = 2, label = '$f(x)$')
ax.hold('on')

nmax = 4
for n in range(1,nmax + 1):
    df_str = 'df' + str(n)
    minimum = np.min(eval(df_str))
    maximum = np.max(eval(df_str))

    if np.abs(minimum) >= np.abs(maximum):
        absmax = np.abs(minimum)
    else:
        absmax = np.abs(maximum)
        
    ax.plot(x,eval(df_str) / absmax, linewidth = 2, linestyle = '--', label = r'$\partial_x^{%i} f$' % n )

ax.hold('off')
frame = ax.get_position() # position of plot center

# shrink frame of the plot to make room for a legend on the right side
ax.set_position([frame.x0, frame.y0,
                 frame.width * 0.8, frame.height])

# Place legend to the right of the shrunken frame
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
          fancybox=True, shadow=True, ncol=1)

plt.grid()
plt.show()
