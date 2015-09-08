import lib
import matplotlib.pyplot as plt

def main(N = 20, steps = 100):

    w = lib.sim.World(N, dt = .5)
    r_histories = lib.sim.evolve(w,steps) # r_histories is 3D: (steps, N, 2)

    # plot routine
    tableau = plot_getpallete()

    # plot all N trajectories
    for n in xrange(N):
        plt.plot(r_histories[:,n,0], r_histories[:,n,1], linewidth = 3, color = tableau[n])

    plt.grid('on')
    plt.show()

def plot_getpallete():
    # reset dictionary values for legend
    params = {'legend.fontsize': 11,
              'legend.linewidth': 2}
    plt.rcParams.update(params)

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

    return tableau20
