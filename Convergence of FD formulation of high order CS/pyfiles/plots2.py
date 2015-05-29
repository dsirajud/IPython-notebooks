import numpy as np
import matplotlib.pyplot as plt
import convergence_routines

ax, bx = -0.5, 1.5
x = np.linspace(ax, bx,1000)
Dx = 0.02
for i in range(6):
    f = convergence_routines.function_edge(x, x_a = -1.025 - i*Dx, x_b = -1.275 - i*Dx, x_c = -1.525 - i*Dx, x_d = 0.475 - i*Dx)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x,f,'-b', linewidth = 2)
    ax.grid()
    plt.xlabel('$x$', fontsize = 14)
    plt.ylabel('$y$', fontsize = 14)
