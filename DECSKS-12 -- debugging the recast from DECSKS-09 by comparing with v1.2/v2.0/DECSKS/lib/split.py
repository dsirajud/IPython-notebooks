import DECSKS
import numpy as np
import time

def scheme(
        f,
        t,x,vx,ax,
        n,
        sim_params
        ):
    """Steps through 1D-1V Vlasov with chosen splitting scheme.

    inputs:
    f -- (ndarray, dim=3) f(t,x,v) with f(0,x,v) = f_0 initialized
    x -- (instance) x, space
    v -- (instance) v, velocity
    t -- (instance) time
    n -- (int) current time step index, t^n
    sim_params -- (dict) simulation parameters

    outputs:
    f -- (ndarray, dim=3) f(n+1,x,v)
    """
    # retrieve sub-dictionary containing splitting coefficients and composition order
    splitting = sim_params['splitting']
    coeff = splitting['order']['coeffs']
    stage = splitting['order']['stages']
    tic = time.time()
    for s in range(len(stage)):
        split_coeff = splitting[coeff[s]][int(stage[s])]
        if s == 0: # on first pass, the previous time step (n - 1) needs to be
            if coeff[s] == 'a': # advect x
                x.CFL.compute_numbers(x, vx, split_coeff*t.width)
                f[n,:,:] = DECSKS.lib.convect.scheme(
                        f[n-1,:,:],
                        n,
                        sim_params,
                        z = x,
                        vz = vx)

            elif coeff[s] == 'b': # advect v
                Ex = DECSKS.lib.fieldsolvers.Gauss1D1V(sim_params['ni'], f, x, vx, n-1, sim_params) # calculate accelerations at time zero (n-1)
                ax.prepointvaluemesh = -Ex
                vx.CFL.compute_numbers(vx, ax, split_coeff*t.width)
                f[n,:,:] = DECSKS.lib.convect.scheme(
                    f[n-1,:,:],
                    n,
                    sim_params,
                    z = vx,
                    vz = ax)

        else: # each subsequent steps overwrites the previous step, all at time n until all split steps complete
            if coeff[s] == 'a': # advect x
                x.CFL.compute_numbers(x, vx, split_coeff*t.width)
                f[n,:,:] = DECSKS.lib.convect.scheme(
                        f[n,:,:],
                        n,
                        sim_params,
                        z = x,
                        vz = vx)

            elif coeff[s] == 'b': # advect v
                Ex = DECSKS.lib.fieldsolvers.Gauss1D1V(sim_params['ni'], f, x, vx, n, sim_params)
                ax.prepointvaluemesh = -Ex
                vx.CFL.compute_numbers(vx, ax, split_coeff*t.width)
                f[n,:,:] = DECSKS.lib.convect.scheme(
                        f[n,:,:],
                        n,
                        sim_params,
                        z = vx,
                        vz = ax)

    toc = time.time()
    print "time step %d of %d completed in %g seconds" % (n,t.N, toc - tic)

    return f
