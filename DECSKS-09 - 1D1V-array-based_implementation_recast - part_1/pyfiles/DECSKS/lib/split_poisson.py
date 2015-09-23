import DECSKS
import time

def scheme(
        f,
        t,x,v,
        n,
        sim_params
        ):
    """Steps through 1D-1V Vlasov by leapfrog (Strang) splitting method.

    inputs:
    f -- (ndarray, dim=3) f(t,x,v) with f(0,x,v) = f_0 initialized
    x -- (instance) x, space
    v -- (instance) v, velocity
    t -- (instance) time
    n -- (int) current time step index, t^n
    sim_params -- (dict) simulation parameters

    outputs:
    f -- (ndarray, dim=3) f(it+1, x,v)
    """

    # retrieve sub-dictionary containing splitting coefficients and composition order
    splitting = sim_params['splitting']
    coeff = splitting['order']['coeffs']
    stage = splitting['order']['stages']
    tic = time.time()
    for s in range(len(stage)):
        split_coeff = splitting[coeff[s]][int(stage[s])]
        if s == 0: # on first pass, the previous time step (n - 1) needs to be
            if coeff[s] == 'a': # convect x
                for j in v.prepoints:
                    x.MCs   = x.generate_Lagrangian_mesh(x.prepointvalues, v.prepointvalues[j], split_coeff*t.width)
                    f[n,:,j] = DECSKS.lib.convect.scheme(
                        f[n-1,:,j],
                        x,n,
                        sim_params)

            elif coeff[s] == 'b': # convect v
                subtic = time.time()
                phi = DECSKS.lib.fieldsolvers.Poisson_PBC_6th(sim_params['ni'], f, x, v, n-1)
                dphi = 1 / x.width ** 1 * sim_params['W_dn1'].dot(phi) # currently W_dn1 is a 6th order LTE matrix of FD coeffs for first derivative

                a = dphi
                for i in x.prepoints:
                    v.MCs   = v.generate_Lagrangian_mesh(v.prepointvalues, a[i], split_coeff*t.width)
                    f[n,i,:] = DECSKS.lib.convect.scheme(
                        f[n-1,i,:],
                        v,n,
                        sim_params)
                subtoc = time.time()
                print "(b) substep %d complete in %g sec" % (s, subtoc - subtic)

        else: # each subsequent steps overwrites the previous step, all at time n until all split steps complete
            if coeff[s] == 'a': # convect x
                subtic = time.time()
                for j in v.prepoints:
                    x.MCs   = x.generate_Lagrangian_mesh(x.prepointvalues, v.prepointvalues[j], split_coeff*t.width)
                    f[n,:,j] = DECSKS.lib.convect.scheme(
                        f[n,:,j],
                        x,n,
                        sim_params)
                subtoc = time.time()
                print "(a) substep %d complete in %g sec" % (s, subtoc - subtic)

            elif coeff[s] == 'b': # convect v
                subtic = time.time()
                phi = DECSKS.lib.fieldsolvers.Poisson_PBC_6th(sim_params['ni'], f, x, v, n)
                dphi =  1 / x.width ** 1 * sim_params['W_dn1'].dot(phi) # currently W_dn1 is a 6th order LTE matrix of FD coeffs for first derivative

                a = dphi
                for i in x.prepoints:
                    v.MCs   = v.generate_Lagrangian_mesh(v.prepointvalues, a[i], split_coeff*t.width)
                    f[n,i,:] = DECSKS.lib.convect.scheme(
                        f[n,i,:],
                        v,n,
                        sim_params)
                subtoc = time.time()
                print "(b) substep %d complete in %g sec" % (s, subtoc - subtic)

    toc = time.time()
    print "time step %d of %d completed in %g seconds" % (n,t.N, toc - tic)

    return f
