import DECSKS
import time

def scheme(
        fe, fi,
        t,
        x, vx, ax,
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

    # compute all CFL numbers for configuration variables beforehand
    x.CFL.compute_all_numbers(sim_params, x, vx, t)
    c = DECSKS.lib.HOC.correctors_on_configuration(sim_params, x, vx, t)

    # retrieve sub-dictionary containing splitting coefficients and composition order
    splitting = sim_params['splitting']
    coeff = splitting['order']['coeffs']
    stage = splitting['order']['stages']
    tic = time.time()
    for s in range(len(stage)):
        split_coeff = splitting[coeff[s]][int(stage[s])]
        if coeff[s] == 'a': # advect x
            # the number of cells a prepoint density packet at [i,j] travels and
            # which direction it travels in is given by the CFL numbers
            # themselves, x.CFL.numbers[stage[s],i,j]

            # when evolving a configuration variable, the distance travelled and
            # direction depends on the physical velocity and time step (e.g. "vx*dt").
            # Since the velocity is a grid quantity (fixed) and all time substeps
            # are known beforehand (uniform time step, and fractional [perhaps negative]
            # steps taken according to the chosen split scheme)
            # the same CFL numbers are reused in every full time step to evolve the
            # configuration; the only change is the density values at each [i,j]
            # from the previous time steps (advections along characteristics)

            # To save cost, we compute all CFL numbers above and reuse them rather than
            # compute the same set in each stage of the split scheme for every full time step.
            # We take further opportunity here to obtain the corresponding
            # high order corrected terms c (high order flux = c.dot(d), c ~ CFL.frac)
            # as well which also are the same from one full time step to the next.
            # hence, we pass the corrector set c below, which is combined inside the
            # convect_configuration.scheme routine with the derivative tensor d to
            # compute the high order flux Uf = c.dot(d). Note that the derivatives
            # must be calculated at each time step as the function, of course, changes
            # from one substep to the next.
            fe = DECSKS.lib.convect_configuration.scheme(
                    fe,
                    int(stage[s]), n,
                    sim_params, c,
                    z = x,
                    vz = vx,
                    split_coeff = split_coeff,
                    charge = -1)

            fi = DECSKS.lib.convect_configuration.scheme(
                    fi,
                    int(stage[s]),n,
                    sim_params, c,
                    z = x,
                    vz = vx,
                    split_coeff = split_coeff,
                    charge = 1)

        elif coeff[s] == 'b': # advect vx
            # calculate electric field at most recent positions of ions and electrons
            Ex = eval(sim_params['compute_electric_field_orchestrator_handle'])(fe, fi, x, vx, sim_params)

            # advect electron velocities
            ax.prepointvaluemesh = -Ex
            vx.CFL.compute_numbers(vx, ax, split_coeff*t.width)
            fe = DECSKS.lib.convect_velocity.scheme(
                fe,
                n,
                sim_params,
                z = vx,
                vz = ax,
                split_coeff = split_coeff,
                charge = -1)

            # advect ion velocities
            ax.prepointvaluemesh = 1. / sim_params['mu'] * Ex
            vx.CFL.compute_numbers(vx, ax, split_coeff*t.width)
            fi = DECSKS.lib.convect_velocity.scheme(
                fi,
                n,
                sim_params,
                z = vx,
                vz = ax,
                split_coeff = split_coeff,
                charge = 1)

    toc = time.time()
    print "time step %d of %d completed in %g seconds" % (n,t.N, toc - tic)

    return fe, fi
