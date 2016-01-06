from pylab import *
import matplotlib.pyplot as plt
import numpy as np

def density(x,vx, xc1 = 5., vc1= 5., xc2 = -5., vc2 = -5.):
    a = 7.
    r1 = np.zeros([len(x), len(vx)])
    for i in range(len(x)):
        for j in range(len(vx)):
            r1[i,j] = np.sqrt( (x[i] - xc1) ** 2
                              + (vx[j] - vc1) ** 2)

    r2 = np.zeros([len(x), len(vx)])
    for i in range(len(x)):
        for j in range(len(vx)):
            r2[i,j] = np.sqrt( (x[i] - xc2) ** 2
                              + (vx[j] - vc2) ** 2)


    f0 = np.zeros([len(x), len(vx)])
    for i in range(len(x)):
        for j in range(len(vx)):
            if r1[i,j] <= a:
                f0[i,j] = (np.cos(np.pi * r1[i,j] / (2.*a) )) ** 22
            elif r2[i,j] <= a:
                f0[i,j] = (np.cos(np.pi * r2[i,j] / (2.*a) )) ** 22
            else:
                f0[i,j] = 0

    return f0

x = np.linspace(-10,10, 512)
vx = np.linspace(-10, 10, 512)
X,V = np.meshgrid(x,vx)

colors = [('white')] + [(cm.jet(i)) for i in xrange(1,256)]
new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=256)

f_left, f_right = np.zeros([len(x), len(vx)]), np.zeros([len(x), len(vx)])
f = density(x,vx)
plt.pcolormesh(X,V,f.T, cmap = new_map)
plt.grid()
plt.title(r'$t = 0$', fontsize = 20)
plt.text(3, .5, r'$f(t = 0, x,v_x)$', fontsize = 18)
plt.text(-7, -1, r'$f^C(t = 0, x,v_x)$', fontsize = 18)
xlabel(r'$x$', fontsize = 18)
ylabel(r'$v_x$', fontsize = 18)
plt.show()
