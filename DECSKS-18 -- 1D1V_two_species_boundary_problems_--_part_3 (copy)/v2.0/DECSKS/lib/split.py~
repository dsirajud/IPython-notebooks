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

            # the advection along characteristics are calculated
            # a priori and stored in x.CFL.numbers[stage[s],:,:],
            # where stage[s] labels the (sub)stage of the full time step.
            # Hence, all the information needed to accomplish the
            # advection of each cell is communicated by passing
            # the stage[s] argument below along with x.CFL.numbers
            # which permits the advection distances to be accessed
            # as needed.

            fe = DECSKS.lib.convect_configuration.scheme(
                    fe,
                    int(stage[s]), n,
                    sim_params, c,
                    z = x,
                    vz = vx,
                    charge = -1)

            fi = DECSKS.lib.convect_configuration.scheme(
                    fi,
                    int(stage[s]),n,
                    sim_params, c,
                    z = x,
                    vz = vx,
                    charge = 1)

        elif coeff[s] == 'b': # advect vx
            # calculate electric field at most recent positions of ions and electrons
            #                Ex = DECSKS.lib.fieldsolvers.compute_electric_field_fourier(fe, fi, x, vx, n-1, sim_params)
            Ex = eval(sim_params['compute_electric_field_function_handle'])(fe, fi, x, vx, sim_params)

            # advect electron velocities
            ax.prepointvaluemesh = -Ex
            vx.CFL.compute_numbers(vx, ax, split_coeff*t.width)
            fe = DECSKS.lib.convect_velocity.scheme(
                fe,
                n,
                sim_params,
                z = vx,
                vz = ax,
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
                charge = 1)

    toc = time.time()
    print "time step %d of %d completed in %g seconds" % (n,t.N, toc - tic)

    return fe, fi
