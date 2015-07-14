import numpy as np
import matplotlib.pyplot as plt


# plot tableau

tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
         (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
         (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
         (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
         (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

N = 6
simtimes = np.zeros(N + 1) # times in minutes

simtimes[1] = 1044.48 / 60
# s7-01 F22, Nx8Nv256Nt60, O6-4

plt.plot(1, simtimes[1], marker = 'o', color = tableau20[1], markersize = 15, linestyle = 'None',
         label = r'F21/O6-4 $N_x = 8, N_v = 256, N_t = 60$')

simtimes[2] = 48091.4 / 60
# s7-02, F22, Nx16Nv512Nt600, O11-6
plt.plot(2, simtimes[2], marker = 'o', color = tableau20[2], markersize = 15, linestyle = 'None',
         label = r'F21/O11-6 $N_x = 16, N_v = 512, N_t = 600$')

simtimes[3] = 52.6565
# s7-03, FD7, Nx384Nv256Nt60, LF2
plt.plot(3, simtimes[3], marker = 'o', color = tableau20[3], markersize = 15,linestyle = 'None',
         label = r'FD7/O6-4 $N_x = 384, N_v = 256, N_t = 60$')

simtimes[4] = 1466.77
# s7-04, FD7, Nx384Nv512Nt200, O6-4
plt.plot(4, simtimes[4], marker = 'o', color = tableau20[4], markersize = 15, linestyle = 'None',
         label = r'FD7/O6-4 $N_x = 384, N_v = 256, N_t = 200$')

simtimes[5] = 479.047
# s7-05, FD7, Nx768Nv256Nt60, O6-4
plt.plot(5, simtimes[5], marker = 'o', color = tableau20[5], markersize = 15, linestyle = 'None',
         label = r'FD7/O6-4 $N_x = 768, N_v = 256, N_t = 60$')

simtimes[6] = 7.96783
# s7-06, F7, Nx16Nv256Nt60, O6-4
plt.plot(6, simtimes[6], marker = 'o', color = tableau20[6], markersize = 15, linestyle = 'None',
         label = r'F7/O6-4 $N_x = 16, N_v = 256, N_t = 60$')


plt.legend(loc = 'best')
plt.grid()
plt.axis([0, N, 0, 4000])
plt.title('Simulation times for various test cases from notebook DECSKS-07', fontsize = 14)
plt.xlabel('simulation number, s7-#', fontsize = 14)
plt.ylabel('total processor time [minutes]', fontsize = 14)
plt.savefig('../fig/simtimes.png')
