import DECSKS
import numpy as np
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
    x.CFL.compute_numbers_for_all_stages(
        sim_params,
        x,
        vx,
        t
        )

    cx = DECSKS.lib.HOC.compute_all_correctors_on_a_configuration_variable(
        sim_params,
        x,
        vx,
        t
        )
    # retrieve sub-dictionary containing splitting coefficients and composition order
    splitting = sim_params['splitting']
    coeff = splitting['order']['coeffs']
    stage = splitting['order']['stages']

    tic = time.time()
    for s in range(len(stage)):
        split_coeff = splitting[coeff[s]][int(stage[s])]
        if coeff[s] == 'a': # advect x

            fe = DECSKS.lib.convect_configuration.scheme(
                    fe,
                    t,
                    n,
                    int(stage[s]),
                    sim_params,
                    cx[int(stage[s]),:,:], # pass the 2D array coresponding to current stage
                    z = x,
                    vz = vx,
                    dt = split_coeff*t.width,
                    charge = -1)

            fi = DECSKS.lib.convect_configuration.scheme(
                    fi,
                    t,
                    n,
                    int(stage[s]),
                    sim_params,
                    cx[int(stage[s]),:,:],
                    z = x,
                    vz = vx,
                    dt = split_coeff*t.width,
                    charge = 1)

    
            # FACTOR IN SOURCES
            if sim_params['BC']['f']['x']['upper'] == 'source':
                # add on its remap contribution for this time step
                fe += sim_params['sources']['f_Se'][int(stage[s]),:,:] # fe = fe[s,i,j] where s is the stage
                fi += sim_params['sources']['f_Si'][int(stage[s]),:,:] # fi = fi[s,i,j] where s is the stage

                # update wall charge due to incoming sources if lower boundary is a collector
                if sim_params['BC']['f']['x']['lower'] == 'collector': 
                    sim_params['sigma'][x.str]['lower'] += -1 * sim_params['sources']['N_Se'][int(stage[s])]
                    sim_params['sigma'][x.str]['lower'] +=  1 * sim_params['sources']['N_Si'][int(stage[s])]
                else:
                    pass

            else:
                pass

        elif coeff[s] == 'b': # advect v
            # calculate electric field at most recent positions of ions and electrons
            Ex = eval(sim_params['compute_electric_field_orchestrator_handle']['x'])(fe, fi, x, vx, n, sim_params)


            # advect electron velocities
            ax.prepointvaluemesh = -Ex
            vx.CFL.compute_numbers(vx, ax, split_coeff*t.width)

            fe = DECSKS.lib.convect_velocity.scheme(
                fe,
                t,
                n,
                sim_params,
                z = vx,
                vz = ax,
                dt = split_coeff*t.width,
                charge = -1)

            # advect ion velocities
            ax.prepointvaluemesh = 1. / sim_params['mu'] * Ex
            vx.CFL.compute_numbers(vx, ax, split_coeff*t.width)

            fi = DECSKS.lib.convect_velocity.scheme(
                fi,
                t,
                n,
                sim_params,
                z = vx,
                vz = ax,
                dt = split_coeff*t.width,
                charge = 1)

    toc = time.time()
    print "time step %d of %d completed in %g seconds" % (n,t.N, toc - tic)

    return fe, fi
