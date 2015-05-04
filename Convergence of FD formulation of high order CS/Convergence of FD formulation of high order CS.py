# the function advected is below

    f = 0.5*np.exp(-((x + 0.25) / 0.03)**2) + 0.5*np.exp(-((x + 0.2) / 0.03)**2) + \
      np.exp(-(x / 0.06)**2) + 0.5*np.exp(-((x - 0.2) / 0.03)**2) + 0.5*np.exp(-((x - 0.25) / 0.03)**2)

    df1 = (-2/0.03 * (x + 0.25) / 0.03)* 0.5*np.exp(-((x + 0.25) / 0.03)**2) + \
      (-2 / 0.03 * (x + 0.2) / 0.03) * 0.5*np.exp(-((x + 0.2) / 0.03)**2) + \
      (-2 / 0.06 * (x / 0.06)) * np.exp(-(x / 0.06)**2) + \
      (-2 / 0.03 * (x - 0.2 / 0.03)) * 0.5*np.exp(-((x - 0.2) / 0.03)**2) + \
      (-2 / 0.03 * (x - 0.25) / 0.03) * 0.5*np.exp(-((x - 0.25) / 0.03)**2)



%matplotlib inline
import numpy as np

x = np.linspace(-0.5,1.5,1000)
f = 0.5*np.exp(-((x + 0.25) / 0.03)**2) + 0.5*np.exp(-((x + 0.2) / 0.03)**2) + np.exp(-(x / 0.06)**2) + 0.5*np.exp(-((x - 0.2) / 0.03)**2) + 0.5*np.exp(-((x - 0.25) / 0.03)**2)

plt.plot(x,f,'-b', linewidth = 2)
plt.grid()
plt.xlabel('$x$', fontsize = 14)
plt.ylabel('$y$', fontsize = 14)
