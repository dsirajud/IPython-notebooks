import os
import time
import multiprocessing as mp
import numpy as np
import lib

def test_fn(evolve_fn, name, steps=1000, dt=1e-3, bodies=101, threads=1):
    """
    evolve_fn -- (function)
    name      -- (str) name of function above for print calls
    steps     -- (int) number of steps in sim
    dt        -- (float) time width
    bodies    -- (int) number N of bodies in N-body sim
    threads   -- (int) number of threads for MP threading
    """
    print "\n"

    # Test the speed of the evolution function evolve_fn,
    w = lib.World(bodies, threads=threads, dt=dt)

    t0 = time.time()
    evolve_fn(w, steps)
    t1 = time.time()
    simtime = t1 - t0
    print "{0} ({1}): {2} steps/sec".format(
        name, threads, int(steps / simtime))

    print "for %s, simtime = %g seconds for %d steps" % (name, simtime, steps)
    # precision checking:

    # Compare the evolution function evolve_fn to the pure python version.
    # for case of 10 bodies

    w1 = lib.World(10, threads=threads, dt=dt)
    w2 = w1.copy()

    lib.evolve(w1, 1024) # Python version
    evolve_fn(w2, 1024)  # evolution function indicated in args above

    def f(name):
        wA = w1
        wB = w2
        dvmax = eval("np.absolute(wA.{0} - wB.{0}).max()".format(name))
        print("    max(delta {0}): {1:2.2}".format(name, dvmax))

    f("r")
    f("v")
    f("F")


if __name__ == "__main__":
    # Single CPU only tests.
    test_fn(lib.evolve, "Python", steps=512)
    test_fn(lib.evolve_c_simple1, "C Simple 1", steps=32000)
    test_fn(lib.evolve_c_simple2, "C Simple 2", steps = 32000)
    test_fn(lib.evolve_c_omp1, "C OMP 1", steps = 32000)


     # Multi-threaded tests. 
    threads = 0

    while True:

        threads += 1
        if threads > mp.cpu_count() + 2:
            break

        steps = threads * 32000

        test_fn(
            lib.evolve_c_omp1, "C OpenMP 1", steps=steps, threads=threads)

    print "\n"
