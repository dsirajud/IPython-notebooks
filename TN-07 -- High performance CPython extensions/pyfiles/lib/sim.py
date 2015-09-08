import numpy as np
import evolve1

class World(object):
    """World is a structure that holds the state of N bodies and
    additional variables.

    dt      : (float) The time-step.

    STATE OF THE WORLD:

    N : (int) The number of bodies in the simulation.
    m : (1D ndarray) The mass of each body.
    r : (2D ndarray) The position of each body.
    v : (2D ndarray) The velocity of each body.
    F : (2D ndarray) The force on each body.

    TEMPORARY VARIABLES:

    R  : (2D ndarray) The vectors from one body to all others.
    R3 : (1D ndarray) The norm of each s vector.
    """
    def __init__(self, N,
                 m_min=1, m_max=30.0, r_max=50.0, v_max=4.0, dt=1e-3):

        self.N  = N
        self.m  = np.random.uniform(m_min, m_max, N)
        self.r  = np.random.uniform(-r_max, r_max, (N, 2)) # sets up N duples: (x,y)
        self.v  = np.random.uniform(-v_max, v_max, (N, 2)) # sets up N duples: (vx,vy)
        self.F  = np.zeros_like(self.r) # force at each r(i) (sum over all j =/= i)
        self.R  = np.zeros_like(self.r) # vector of duples that are the distance from self to all others
        self.R3 = np.zeros_like(self.m) # (2-norm) ** 3, collapses above duples to collection of scalars

        self.dt = dt

def compute_F(w):
    """Compute the force on each body in the world, w."""
    for i in xrange(w.N):
        w.R[:] = w.r - w.r[i]
        w.R3[:] = (w.R[:,0]**2 + w.R[:,1]**2)**1.5
        w.R3[i] = 1.0 # This makes the self-force zero.
        w.F[i] = (w.m[i] * w.m[:,None] * w.R / w.R3[:,None]).sum(0)

def evolve(w, steps):
    """Evolve the world, w, through the given number of steps."""

    steps += 1 # so that iterating below has step 1 as 1, step 2 as 2, ...

    # initialize space curve history
    r_histories = np.zeros([steps, w.N, 2])
    r_histories[0,:,:] = w.r

    for _ in xrange(1,steps):
        compute_F(w)
        # velocities incremented according to forces w.F
        w.v += w.F * w.dt / w.m[:,None]

        # positions incremented according to velocities w.v
        w.r += w.v * w.dt

        # track history (space curves)
        r_histories[_, :, :] = w.r

    return r_histories

def evolve1_c_(w, steps):
    """Evolve the world using the evolve1 C module."""
    evolve1.evolve(w.dt, steps, w.N, w.m, w.r, w.v, w.F)

